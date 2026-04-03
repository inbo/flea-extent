#' Convert a single point location to a grid cell polygon
#'
#' @param xy a `SpatVector` with geometry type points
#' @param cell_width_m cell width in meter, default 500
#' @param point_position default centre of grid cell
#' @param crs default EPSG code 31370
#'
#' @return a `SpatVector` with geometry type polygon
#' @export
#'
#' @examples
point_to_gridcell <- function(
    xy,
    cell_width_m = 500,
    point_position =
      c("center", "lowerleft", "upperleft", "lowerright", "upperright"),
    crs = 31370) {
  point_position <- match.arg(point_position)

  if (point_position != "center") stop(point_position, " not yet implemented")

  stopifnot(inherits(xy, "SpatVector"))
  # convert to sf
  xy <- sf::st_as_sf(xy)
  xy_df <- sf::st_drop_geometry(xy)
  xy <- sf::st_geometry(xy)

  # buffer with 1 point per quandrant
  halflength <- cell_width_m / 2
  xy_buffer <- sf::st_buffer(
    x = xy,
    dist = sqrt(2 * halflength^2),
    nQuadSegs = 1
  )

  # rotate 45 degrees around centroid
  rot <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
  pl <- (xy_buffer - xy) * rot(pi / 4) + xy
  pl <- sf::st_sf(data.frame(xy_df, pl), crs = crs)
  # convert to SpatVector
  pl <- terra::vect(pl)
  return(pl)
}

#' Extract a sample of points from a raster
#'
#' This function extracts a specified number of points from a raster,
#' selecting the lowest values after sorting.
#'
#' @param rast A raster object (`terra` `SpatRaster`) containing GRTS rank
#' numbers
#' @param n The number of points to extract
#'
#' @return `SpatVector` object containing the extracted points with their values,
#'         cell numbers, coordinates, and assigned stratum name
#'
#' @importFrom terra extract as.points
#' @importFrom sf st_as_sf
#'
#' @examples
#' # Assuming 'my_raster' is a terra SpatRaster object
#' # result <- extract_sample_helper(my_raster, 100)
#'
#' @export
extract_sample_helper <- function(
    rast,
    n) {

  # Extract values, exclude NA
  extracted <- terra::extract(
    x = rast,
    y = as.points(rast, na.rm = TRUE),
    cells = TRUE,
    xy = TRUE,
    ID = FALSE
  )

  # Sort and select the lowest n
  sorted_indices <- order(extracted[[1]])[1:n]
  selected <- extracted[sorted_indices, ]
  selected$stratum_name <- names(selected)[1]
  names(selected)[1] <- "grts_rank"

  # Convert to SpatVector
  selected <- selected |>
    vect(geom = c("x", "y"), crs = "EPSG:31370")

  return(selected)
}


#' Create raster with layers corresponding to temporal strata and values
#' corresponding to GRTS rankings from the master GRTS sample
#'
#' @param stratum_raster A raster object representing the stratification
#' @param fleagrts A raster object containing GRTS rankings
#' @param stratum_name Character string specifying the name of the stratum
#' @param min_stratum_size Integer specifying the minimum stratum size in terms
#' of number of raster cells
#'
#' @return A `SpatRaster`
#'
#' @importFrom terra `activeCat<-` classify segregate crop
#' @importFrom dplyr %>% mutate
#' @importFrom tibble as_tibble
#'
#' @export
separate_grts_strata <- function(
    stratum_raster,
    fleagrts,
    stratum_name) {

  # assertions
  assertthat::assert_that(inherits(stratum_raster, "SpatRaster"))
  assertthat::assert_that(inherits(fleagrts, "SpatRaster"))
  assertthat::assert_that(terra::compareGeom(stratum_raster, fleagrts))
  assertthat::assert_that(assertthat::is.string(stratum_name))

  # Set the active category to the desired factor
  ts1 <- `activeCat<-`(stratum_raster, stratum_name)
  # Adjust the number of levels
  # select second column which holds the gain, loss, ... cats
  levelvec <- levels(ts1)[[1]][[2]]
  # set stable absence to NA
  levelvec[levelvec == "Stable absence"] <- NA
  # reclassify original levels (all transitions) to just the four change cats
  ts2 <- classify(
    ts1,
    cbind(
      levels(ts1)[[1]][[1]],
      as.numeric(as.factor(levelvec))
    )
  )
  names(ts2) <- "changecat"
  levels(ts2) <- data.frame(
    changecat = seq_along(levels(as.factor(levelvec))),
    label = levels(as.factor(levelvec))
  )

  # make separate layer for each change cat
  ts2 <- ts2 %>% segregate(classes = NULL, other = NA)
  names(ts2) <- levels(as.factor(levelvec))

  # make separate layers containing grts rankings for each change cat
  fleagrts_ts2 <- crop(fleagrts, ts2, mask = TRUE)
  names(fleagrts_ts2) <- names(ts2)
  return(fleagrts_ts2)
}



#' Extract stratified sample from the cropped set of GRTS
#' (Generalized Random Tessellation Stratified) rankings
#'
#' @param ntot Integer specifying the total number of samples to extract
#' @param nmin Integer specifying the minimum number of samples per stratum
#'
#' @return A `SpatVector`
#'
extract_sample <- function(
    separate_grts,
    stratum_name,
    ntot,
    nmin,
    min_stratum_size) {
  assertthat::assert_that(inherits(separate_grts, "SpatRaster"))
  assertthat::assert_that(is.numeric(ntot), ntot > 0)
  assertthat::assert_that(is.numeric(nmin), nmin < ntot)
  assertthat::assert_that(
    is.numeric(min_stratum_size),
    min_stratum_size > 0
  )

  # determine stratum population sizes
  popsize <- global(separate_grts, fun = "notNA") |>
    as_tibble(rownames = "layername") |>
    rename(count = notNA)

  # check for strata that are too small
  remove_me <- popsize$layername[popsize$count < min_stratum_size]
  if (length(remove_me) > 0) {
    popsize <- popsize[!(popsize$layername %in% remove_me), ]
    if (nrow(popsize) > 0) {
      separate_grts <- terra::subset(separate_grts, remove_me, negate = TRUE)
    } else {
      return(terra::vect())
    }
  }

  # determine sample size allocation
  # first distribute nmin to each stratum,
  # remaining allocate proportional to stratum size
  # give more weight to changes than to stable?
  ntot <- ntot * nrow(popsize) / 4 # reduce ntot in case less than 4 changeclass

  allocation <- popsize |>
    mutate(
      stable = layername == "Stable presence",
      n_changeclasses = n(),
      ntot_stable = round(ntot * 1 / n_changeclasses),
      ntot_changed = ntot - ntot_stable,
      n_h = ifelse(
        stable,
        ntot_stable,
        nmin +
          round((ntot_changed - nmin * (n() - 1)) *
                  (count / sum(count[!stable]))) # nolint
      )
    )


  sample_ts2 <- vector(mode = "list", length = nlyr(separate_grts))
  sample_ts2 <- setNames(sample_ts2, names(separate_grts))
  for (i in names(sample_ts2)) {
    sample_ts2[[i]] <- extract_sample_helper(
      separate_grts[[i]],
      allocation$n_h[allocation$layername == i]
    )
    sample_ts2[[i]]$changecat <- i
  }

  sample_ts2 <- vect(sample_ts2)
  sample_ts2$stratum_name <- stratum_name

  return(sample_ts2)
}



#' Get data from the web feature service for GRB
#'
#' @param layer A string. Should be one of the `WFS` layers available in the
#' service
#' @param bbox A `SpatExtent` or an object from which a `SpatExtent` can be
#' determined. The `bbox` values should be in CRS 31370.
#'
#' @return A `SpatVector` containing `GRB` layer objects that intersect `bbox`.
#' @export
#'
#' @examples
get_grb <- function(layer, bbox) {
  wfs_grb <- "https://geo.api.vlaanderen.be/GRB/wfs"

  assertthat::assert_that(assertthat::is.string(layer))

  bbox <- terra::ext(bbox)

  bbox_sf <- sf::st_as_sfc(sf::st_bbox(bbox, crs = sf::st_crs(31370)))

  grb <- sf::read_sf(
    paste0("WFS:", wfs_grb),
    layer = layer,
    wkt_filter = sf::st_as_text(bbox_sf)
  )

  grb <- grb |>
    sf::st_cast("GEOMETRYCOLLECTION") %>%
    sf::st_collection_extract("LINESTRING") %>%
    sf::st_cast("POLYGON")

  grb <- terra::vect(grb)

  # crop to bbox
  grb <- terra::crop(grb, bbox)

  return(grb)
}

read_layernames <- function(x) {
  assertthat::assert_that(is.character(x))
  return(x)
}

get_grb_by_row <- function(layer, polygons) {
  assertthat::assert_that(inherits(polygons, "SpatVector"))

  # check if polygon has data
  if (nrow(polygons) == 0) {
    return(polygons)
  }

  out <- vector("list", length = nrow(polygons))
  namesvec <- polygons$grts_rank
  out <- setNames(out, nm = namesvec)
  for (i in seq_along(out)) {
    bbox <- polygons[i, ]
    grb <- get_grb(layer = layer, bbox = bbox)
    grb$grts_rank <-  namesvec[i]
    out[[as.character(namesvec[i])]] <- grb
  }
  out <- terra::vect(out)
  out$layer <- layer
  # convert date(time) fields to ISO-8601 format
  time <- terra::datatype(out) == "time"
  datetimecols <- terra::names(out)[time]
  if (length(datetimecols) > 0) {
    for (i in datetimecols) {
      out[[i]] <-
        format(out[[i]], format = "%Y-%m-%dT%H:%M:%S.000Z", tz = "UTC")
    }
  }
  # catch case empty records
  # can be removed when terra or geotargets deals with this natively
  # https://github.com/ropensci/geotargets/issues/187
  if (nrow(out) == 0) {
    out <- .create_empty_geom(out, type = terra::geomtype(polygons))
    return(out)
  }

  return(out)
}

process_water_wtz <- function(grb_wtz) {
  grb <- vect(grb_wtz)
  grb <- sf::st_as_sf(grb)
  grb <- grb |>
    dplyr::mutate(
      jaar = pmin(lubridate::year(BEGINDATUM),
                  lubridate::year(OPNDATUM), na.rm = TRUE),
      layer = "GRB:WTZ",
      value = ifelse(VHAG %in% c(-9, -8), 900, 800)
    )

  # select columns
  grb <- grb |>
    select(
      gml_id,
      grts_rank,
      jaar,
      layer,
      value
    )

  grb <- vect(grb)

  return(grb)
}

process_parcels <- function(grb) {
  grb <- vect(grb)
  grb <- grb[, c("gml_id", "grts_rank", "BEGINDATUM", "FISCDATUM")]
  grb <- sf::st_as_sf(grb) |>
    dplyr::mutate(
      jaar = pmin(lubridate::year(BEGINDATUM),
                  lubridate::year(FISCDATUM), na.rm = TRUE),
      layer = "GRB:ADP",
      value = NA
    ) |>
    dplyr::select(gml_id, grts_rank, jaar, layer, value)
  # cast to lines
  grb <- terra::vect(grb) |> terra::as.lines()
  return(grb)
}

process_settlement <- function(grb) {
  grb <- vect(grb) # this should combine grb settlement layers
  grb <- sf::st_as_sf(grb)

  # enkel TRN met bodembedekking verhard
  grb <- grb |>
    dplyr::filter(is.na(LBLBDMBD) | LBLBDMBD == "verhard") |>
    tidyr::unite(
      col = lbl,
      c(LBLTYPE, LBLFNCT),
      na.rm = TRUE,
      remove = FALSE) |>
    dplyr::mutate(
      jaar = pmin(lubridate::year(BEGINDATUM),
                  lubridate::year(OPNDATUM), na.rm = TRUE)
    )

  # assign values
  # 101 1.1 Settlements - buildings
  # 102 1.2 Settlements - sealed soil
  # 105 1.5 Settlements - water
  # 106 1.6 Settlements - unknown land cover
  grb <- grb |>
    dplyr::mutate(
      value = case_when(
        layer %in% c("GRB:GBG", "GRB:GBA", "GRB:KNW") ~ 101,
        layer %in% c(
          "GRB:WBN",
          "GRB:SBN") |
          (layer == "GRB:TRN" & LBLBDMBD == "verhard") ~ 102,
        TRUE ~ NA
      )
    )

  # select columns
  grb <- grb |>
    select(
      gml_id,
      grts_rank,
      layer,
      jaar,
      lbl,
      value
    )

  # deal with overlapping polygons
  # within same grts_rank
  grb <- grb |>
    arrange(value, factor(
      layer,
      levels = c(
        "GRB:GBG", "GRB:GBA", "GRB:KNW", "GRB:WBN", "GRB:SBN", "GRB:TRN"
        )
      )
    )
  grb <- vect(grb)
  grb <- unique(grb)

  gbg <- subset(grb, layer == "GRB:GBG", NSE = TRUE) |>
    aggregate(
      by = c(
        "gml_id",
        "grts_rank",
        "layer",
        "jaar",
        "lbl",
        "value"
      ),
      count = FALSE
    )
  gba <- subset(grb, layer == "GRB:GBA", NSE = TRUE) |>
    aggregate(
      by = c(
        "gml_id",
        "grts_rank",
        "layer",
        "jaar",
        "lbl",
        "value"
      ),
      count = FALSE
    )
  knw <- subset(grb, layer == "GRB:KNW", NSE = TRUE) |>
    aggregate(
      by = c(
        "gml_id",
        "grts_rank",
        "layer",
        "jaar",
        "lbl",
        "value"
      ),
      count = FALSE
    )
  wbn <- subset(grb, layer == "GRB:WBN", NSE = TRUE) |>
    aggregate(
      by = c(
        "gml_id",
        "grts_rank",
        "layer",
        "jaar",
        "lbl",
        "value"
      ),
      count = FALSE
    )
  sbn <- subset(grb, layer == "GRB:SBN", NSE = TRUE) |>
    aggregate(
      by = c(
        "gml_id",
        "grts_rank",
        "layer",
        "jaar",
        "lbl",
        "value"
      ),
      count = FALSE
    )
  trn <- subset(grb, layer == "GRB:TRN", NSE = TRUE) |>
    aggregate(
      by = c(
        "gml_id",
        "grts_rank",
        "layer",
        "jaar",
        "lbl",
        "value"
      ),
      count = FALSE
    )

  c1 <- cover(gba, gbg)
  c2 <- cover(knw, c1)
  c3 <- cover(wbn, c2)
  c4 <- cover(sbn, c3)
  c5 <- cover(trn, c4)

  c5 <- unique(c5) # remove duplicate records

  return(c5)
}

mapping_lbg_to_flea <- function() {
  mapping <- googlesheets4::read_sheet(
    ss = "1O_6EFo9hK_VM4MEoBR_f9yufeVFlfw5hdhyPGY9N34Q",
    sheet = "LGP_Reclass"
  ) |>
    dplyr::select(
      GWSNAM_H, GWSCOD_H,
      FLEA = EXT_LEV2_CODE
    ) |>
    dplyr::mutate(
      FLEA = gsub("(^[2-9])0\\d{1}", "\\100", FLEA) |> as.numeric()
    ) |>
    dplyr::distinct()
  if (!any(mapping$GWSCOD_H == 1)) {
    mapping <- mapping |>
      dplyr::add_row(
        GWSCOD_H = 1,
        GWSNAM_H = "Stallen en gebouwen / bedrijfszetel",
        GWSGRPH_LB = "Landbouwinfrastructuur",
        FLEA = 101
      )
  }
  if (!any(mapping$GWSCOD_H == 2)) {
    mapping <- mapping |>
      dplyr::add_row(
        GWSCOD_H = 2,
        GWSNAM_H = "Bijgebouwen",
        GWSGRPH_LB = "Landbouwinfrastructuur",
        FLEA = 101
      )
  }
  return(mapping)
}


get_lbg_layernames <- function(path_to_lbg) {

  lyrs <- terra::vector_layers(path_to_lbg)
  lyrs <- lyrs[grepl(".+20\\d\\d.+", lyrs)]
  return(lyrs)
}

get_lbg <- function(
    path_to_lbg, layer, from_fields, where_field, where_values, flea_value) {

  assertthat::assert_that(assertthat::is.string(layer))
  assertthat::assert_that(is.character(from_fields))
  assertthat::assert_that(assertthat::is.string(where_field))
  assertthat::assert_that(is.numeric(where_values))

  query <- paste0(
    "SELECT ",
    paste0(from_fields, collapse = ","),
    " FROM ",
    layer,
    " WHERE ",
    where_field,
    " IN (",
    paste0("'", where_values, "'", collapse = ","),
    ")"
  )

  lbg <- vect(x = path_to_lbg,
              query = query,
              crs = "EPSG:31370")
  lbg$layer <- layer
  lbg$value <- flea_value
  # convert date(time) fields to ISO-8601 format
  time <- terra::datatype(lbg) == "time"
  datetimecols <- terra::names(lbg)[time]
  if (length(datetimecols) > 0) {
    for (i in datetimecols) {
      lbg[[i]] <-
        format(lbg[[i]], format = "%Y-%m-%dT%H:%M:%S.000Z", tz = "UTC")
    }
  }

  return(lbg)
}


#' Helper functions to deal with empty records in `tar_terra_vect`
.create_empty_geom <- function(x, type) {
  type <- match.arg(gsub("S$", "", toupper(type)), c("POINT", "LINE", "POLYGON"))
  if (nrow(x) == 0) {
    empty <- terra::vect(paste(type, "EMPTY"), crs = terra::crs(x))
    cols <- as.data.frame(x)[1, ]
    x <- cbind(empty, cols)
  }
  x
}

.filter_empty_geom <- function(x) {
  x[!is.na(x), ]
}




spatvector_crop <- function(x, y) {
  assertthat::assert_that(inherits(x, "SpatVector"))
  assertthat::assert_that(inherits(y, "SpatVector"))

  # spatial subset
  x <- x[y]
  # topology fix if needed
  x <- terra::makeValid(x)
  # if y contains overlapping polygons dissolve them
  # this was needed in case validation polygons overlapped
  y <- terra::aggregate(y)
  out <- terra::crop(x, y)
  return(out)
}

download_watersurfaces <- function(path_flea_data, meta) {
  n2khab::fileman_folders(path = path_flea_data)
  path <- file.path(
    path_flea_data,
    "n2khab_data", "10_raw",
    paste0("watersurfaces_", meta$version)
  )
  fs::dir_create(path)
  n2khab::download_zenodo(
    doi = meta$doi,
    path = path,
    quiet = TRUE,
    parallel = FALSE
  )
  return(path)
}


download_habitatmap_terr <- function(path_flea_data, meta) {
  n2khab::fileman_folders(path = path_flea_data)
  path <- file.path(
    path_flea_data,
    "n2khab_data", "20_processed", meta$version
  )
  fs::dir_create(path)
  n2khab::download_zenodo(
    doi = meta$doi,
    path = path,
    quiet = TRUE,
    parallel = FALSE
  )
  return(path)
}


get_watersurfaces <- function(path_version, polygons, meta) {
  #https://inbo.github.io/n2khab/reference/read_watersurfaces.html

  file_version <- switch(
    meta$version,
    "v1.0" = file.path(path_version, "watersurfaces.shp"),
    "v1.1" = file.path(path_version, "watersurfaces.gpkg"),
    "v1.2" = file.path(path_version, "watersurfaces.gpkg"),
    "v2024" = file.path(path_version, "watersurfaces.gpkg")
  )

  ws <- n2khab::read_watersurfaces(
    file = file_version,
    version = basename(path_version),
    fix_geom = TRUE
    )
  ws <- vect(ws)

  ws <- spatvector_crop(x = ws, y = polygons)
  ws$year_flea <- meta$year_flea
  ws$layer <- basename(path_version)
  ws$value <- NA_real_
  ws$area_name <- NULL
  ws$wfd_type_certain <- NULL

  return(ws)
}

get_types_heath <- function() {
  n2khab::read_types() |>
    filter(
      typeclass_name %in% c("Temperate heath and scrub", "Sclerophyllous scrub")
    ) |>
    pull(type)
}
get_types <- function(flea_id) {
  if (flea_id == 500) {
    n2khab::read_types() |>
      filter(
        typeclass_name %in% c("Temperate heath and scrub", "Sclerophyllous scrub")
      ) |>
      pull(type)
  }
}


get_habitatmap_terr <- function(
    path_version, polygons, meta, types, min_phab, flea_value
) {
  assertthat::assert_that(is.character(path_version), is.factor(types))
  assertthat::assert_that(inherits(polygons, "SpatVector"))
  assertthat::assert_that(is.data.frame(meta))
  assertthat::assert_that(
    is.numeric(min_phab), min_phab > 0, min_phab <= 100,
    is.numeric(flea_value)
  )

  file_version <- switch(
    meta$version,
    "habitatmap_terr_2018_v2" = file.path(path_version, "habitatmap_terr.gpkg"),
    "habitatmap_terr_2020_v2" = file.path(path_version, "habitatmap_terr.gpkg"),
    "habitatmap_terr_2023_v1" = file.path(path_version, "habitatmap_terr.gpkg")
  )

  hmt <- n2khab::read_habitatmap_terr(
    file = file_version,
    version = basename(path_version),
    keep_aq_types = FALSE,
    drop_7220 = TRUE
  )
  hmtp <- hmt$habitatmap_terr_polygons
  hmtt <- hmt$habitatmap_terr_types

  # keep only selected types
  hmtt <- hmtt |> filter(type %in% types, phab >= min_phab)
  hmtp <- hmtp |> inner_join(hmtt, by = "polygon_id")


  hmtp <- vect(hmtp)
  hmtp <- spatvector_crop(x = hmtp, y = polygons)
  if (nrow(hmtp) > 0) {
    hmtp$year_flea <- meta$year_flea
    hmtp$layer <- basename(path_version)
    hmtp$value <- flea_value
  }
  return(hmtp)
}



combine_grb_inbo_water <- function(grb_water, inbo_water, meta) {
  inbo_water <- vect(inbo_water)
  inbo_water <- inbo_water[
    inbo_water$layer == paste0("watersurfaces_",meta$version), ]

  # cover: values of x that overlap with y are replaced by y
  grb_water <- aggregate(x = grb_water, by = names(grb_water))
  grb_water <- makeValid(grb_water)
  inbo_water <- makeValid(inbo_water)
  water <- safe_cover(x = grb_water, y = inbo_water)
  water <- unique(water)
  water <- makeValid(water)
  return(water)
}




combine_sources <- function(
    water, settlements, polygons,
    lbg_101, lbg_104, lbg_200, lbg_300, lbg_400, lbg_500, lbg_900,
    terr_500
    ) {

  assertthat::assert_that(inherits(water, "SpatVector")) # a branch
  assertthat::assert_that(inherits(settlements, "SpatVector")) # a target
  assertthat::assert_that(inherits(lbg_101, "list")) # a pattern
  assertthat::assert_that(inherits(lbg_104, "list")) # a pattern
  assertthat::assert_that(inherits(lbg_200, "list")) # a pattern
  assertthat::assert_that(inherits(lbg_300, "list")) # a pattern
  assertthat::assert_that(inherits(lbg_400, "list")) # a pattern
  assertthat::assert_that(inherits(lbg_500, "list")) # a pattern
  assertthat::assert_that(inherits(lbg_900, "list")) # a pattern
  assertthat::assert_that(inherits(terr_500, "list")) # a pattern
  assertthat::assert_that(inherits(polygons, "list")) # a pattern

  vp <- vect(polygons)
  lbg_101 <- vect(lbg_101) |> unique()
  lbg_104 <- vect(lbg_104) |> unique()
  lbg_200 <- vect(lbg_200) |> unique()
  lbg_300 <- vect(lbg_300) |> unique()
  lbg_400 <- vect(lbg_400) |> unique()
  lbg_500 <- vect(lbg_500) |> unique()
  lbg_900 <- vect(lbg_900) |> unique()

  terr_500 <- vect(terr_500) |> unique()

  # get the validation year
  year_to_validate <- unique(water$year_flea)
  year_to_validate <- year_to_validate[!is.na(year_to_validate)]

  # filter the layers to only the validation year
  lbg_101 <- lbg_101[grepl(year_to_validate, x = lbg_101$layer), ]
  lbg_104 <- lbg_104[grepl(year_to_validate, x = lbg_104$layer), ]
  lbg_200 <- lbg_200[grepl(year_to_validate, x = lbg_200$layer), ]
  lbg_300 <- lbg_300[grepl(year_to_validate, x = lbg_300$layer), ]
  lbg_400 <- lbg_400[grepl(year_to_validate, x = lbg_400$layer), ]
  lbg_500 <- lbg_500[grepl(year_to_validate, x = lbg_500$layer), ]
  lbg_900 <- lbg_900[grepl(year_to_validate, x = lbg_900$layer), ]

  terr_500 <- terr_500[
    terr_500$year_flea == year_to_validate,
    c("polygon_id", "year_flea", "layer", "value")
  ]

  all_types <- data.frame(
    colname = c(
      "grts_rank", "cell", "stratum_name",
      "changecat", "gml_id", "layer", "jaar", "lbl", "value", "year_flea",
      "agg_n", "polygon_id", "wfd_code", "hyla_code", "name", "wfd_type",
      "depth_class", "connectivity", "usage", "wfd_type_alternative",
      "water_level_management", "GWSCOD_H", "GWSNAM_H"
    ),
    type = c(
      "numeric",
      "numeric", "character", "character", "character", "character",
      "numeric", "character", "numeric", "numeric", "numeric", "character",
      "character", "numeric", "character", "character", "character",
      "character", "character", "character", "character", "character",
      "character"
    )
  )

  vplist <- vector("list", nrow(vp))
  for (i in seq_along(vp)) {
    print(sprintf("%s out of %s done", i, nrow(vp)))
    vp_ <- vp[i]
    vp_ <- makeValid(vp_) # Ensure base polygon is valid

    # Process Water
    w_ <- water[vp_]
    w_ <- crop(w_, vp_)
    if (nrow(w_) > 0) {
      w_ <- subset(w_, expanse(w_) > 1)
      w_ <- subset(w_, w_$grts_rank == vp_$grts_rank | is.na(w_$grts_rank))
      w_$value <- as.numeric(w_$value)
      w_ <- makeValid(w_)
    }

    # Process Settlements
    s_ <- settlements[vp_]
    s_ <- crop(s_, vp_)
    s_ <- s_[s_$grts_rank == vp_$grts_rank, ]
    s_ <- subset(s_, !(s_$jaar > year_to_validate & grepl("^GRB", s_$layer)))
    if (nrow(s_) > 0) s_ <- makeValid(s_)

    # Process Agriculture
    lbg_101_ <- crop(lbg_101[vp_], vp_)
    lbg_104_ <- crop(lbg_104[vp_], vp_)
    lbg_200_ <- crop(lbg_200[vp_], vp_)
    lbg_300_ <- crop(lbg_300[vp_], vp_)
    lbg_400_ <- crop(lbg_400[vp_], vp_)
    lbg_500_ <- crop(lbg_500[vp_], vp_)
    lbg_900_ <- crop(lbg_900[vp_], vp_)

    lbg_ <- rbind(
      lbg_101_, lbg_104_, lbg_200_, lbg_300_, lbg_400_, lbg_500_, lbg_900_
    )
    if (nrow(lbg_) > 0) lbg_ <- makeValid(lbg_)

    # Process habitatmap_terr
    terr_500_ <- crop(terr_500[vp_], vp_)

    # Hierarchical Overlay (Z-Order)
    # 1. Initialize the base layer
    c2 <- vp_

    # 2. Sequentially cover only if the layer actually contains data
    if (nrow(lbg_) > 0) {
      c2 <- cover(c2, lbg_)
      c2 <- makeValid(c2)
    }

    if (nrow(terr_500_) > 0) {
      c2 <- cover(c2, terr_500_)
      c2 <- makeValid(c2)
    }

    if (nrow(s_) > 0) {
      c2 <- cover(c2, s_)
      c2 <- makeValid(c2)
    }

    if (nrow(w_) > 0) {
      c2 <- cover(c2, w_)
      c2 <- makeValid(c2)
    }

    # 3. Finally, disaggregate the combined result
    c2 <- disagg(c2)

    # Pre-dissolve geometry flattening
    if (nrow(c2) > 1) {
      c2 <- unique(c2)
      c2 <- makeValid(c2)

      # possibly still overlapping polygons (with slight differences in geom)
      if (nrow(c2) > 1) {
        result <- c2[1, ]
        for (k in 2:nrow(c2)) {
          part_k <- makeValid(c2[k, ])
          result <- cover(result, part_k)
          result <- makeValid(result)
        }
        c2 <- result
      }
    }

    out <- c2
    out <- disagg(out)

    out_areas <- expanse(out)
    out <- subset(out, !is.na(out_areas) & out_areas > 0.01)
    out <- makeValid(out)

    # Enforce base vp_ attributes globally on the resulting pieces
    out$grts_rank <- vp_$grts_rank
    out$cell <- vp_$cell
    out$stratum_name <- vp_$stratum_name
    out$changecat <- vp_$changecat
    out$year_flea <- year_to_validate

    stopifnot(nrow(out) >= 1)
    vplist[[i]] <- out
  }

  all_cols <- all_types$colname

  # Column Standardization
  vplist_complete <- lapply(vplist, function(x) {
    missing_cols <- setdiff(all_cols, names(x))
    if (length(missing_cols) > 0) {
      for (col in missing_cols) {
        target_type <- all_types$type[all_types$colname == col]
        x[[col]] <- as(NA, target_type)
      }
    }
    x <- x[, all_cols]  # Reorder to match

    x_vals <- values(x)
    for (col in all_cols) {
      target_type <- all_types$type[all_types$colname == col]
      x_vals[[col]] <- as(x_vals[[col]], target_type)
    }
    values(x) <- x_vals

    return(x)
  })

  vp_wa_se <- vect(vplist_complete)

  # Final unique pass
  vp_wa_se <- terra::unique(vp_wa_se)
  vp_wa_se <- makeValid(vp_wa_se)
  vp_wa_se <- terra::disagg(vp_wa_se)

  # remove any microscopic artifacts < 0.1 sq meters remaining
  areas <- expanse(vp_wa_se)
  vp_wa_se <- subset(vp_wa_se, !is.na(areas) & areas > 0.1)

  return(vp_wa_se)
}

postprocess_prelabelling <- function(
    water_grb_lbg,
    settlement_mask,
    mmu = 30
) {
  # add column indicating fraction of poly within settlement
  settl_prop <- exactextractr::exact_extract(
    x = settlement_mask,
    y = st_as_sf(water_grb_lbg),
    fun = 'mean',
    progress = FALSE
  )
  water_grb_lbg$prop_in_settlement <- settl_prop

  ws <- sf::st_as_sf(water_grb_lbg)

  # get the validation year
  year_to_validate <- unique(ws$year_flea)
  year_to_validate <- year_to_validate[!is.na(year_to_validate)]

  ws <- ws |>
    dplyr::select(
      grts_rank,
      layer,
      year_grb = jaar,
      stratum_name,
      changecat,
      value,
      year_flea,
      prop_in_settlement
    ) |>
    dplyr::group_by(grts_rank) |>
    dplyr::mutate(
      stratum_name = ifelse(
        is.na(stratum_name), dplyr::first(stratum_name), stratum_name),
      changecat = ifelse(
        is.na(changecat), dplyr::first(changecat), changecat),
      label = dplyr::case_when(
        !is.na(value) ~ as.character(value),
        layer %in% c("GRB:GBA", "GRB:GBG") ~ "101",
        layer %in% c("GRB:WBN", "GRB:SBN", "GRB:KNW", "GRB:TRN") ~ "102",
        grepl("watersurfaces", layer) & prop_in_settlement > 0.5 ~ "105",
        grepl("watersurfaces", layer) & prop_in_settlement <= 0.5 ~ "900",
        layer == "GRB:WTZ" ~ "water",
        TRUE ~ "other"
      ),
      year_flea = year_to_validate
    ) |>
    dplyr::ungroup()

  # --------------------------------------------------------------------
  # MMU Enforcement: Longest Shared Border Assignment
  # --------------------------------------------------------------------

  ws <- terra::vect(ws)
  group_cols <- c(
    "grts_rank", "stratum_name", "changecat", "year_flea", "label"
  )
  aggregate_cols <- setdiff(names(ws), group_cols)
  # need to isolate each validation polygons
  # to avoid slivers from one stratum get dissolved into another stratum
  stratum_cols <- c("grts_rank", "stratum_name", "changecat", "year_flea")

  # Get unique strata
  strata_df <- unique(values(ws)[, stratum_cols])

  # Apply dissolve function to each stratum
  processed_list <- vector(mode = "list", length = nrow(strata_df))
  for (strat in seq_len(nrow(strata_df))) {
    print(sprintf("Stratum number %s out of %s", strat, nrow(strata_df)))
    # Subset to just this stratum
    sub_poly <- filter_spatvector(ws, strata_df[strat, ])
    # Run dissolve function
    # terra_dissolve_boundaries <- purrr::possibly(
    #   terra_dissolve_boundaries, otherwise =  sub_poly
    # )
    processed_list[[strat]] <- terra_dissolve_boundaries(
      polygons = sub_poly,
      group_cols = group_cols,
      mmu = mmu
    )
  }
  # Combine the results back into a single SpatVector
  ws <- do.call(rbind, processed_list)
  ws <- terra::unique(ws)
  return(ws)
}


filter_spatvector <- function(v, df) {
  # 1. Validate inputs
  assertthat::assert_that(inherits(v, "SpatVector"),
                          msg = "Input 'v' must be a terra SpatVector.")
  assertthat::assert_that(is.data.frame(df) && nrow(df) == 1,
                          msg = "Input 'df' must be a 1-row data.frame.")

  filter_cols <- names(df)

  # 2. Check that all filter columns actually exist in the SpatVector
  missing_cols <- setdiff(filter_cols, names(v))
  assertthat::assert_that(
    length(missing_cols) == 0,
    msg = sprintf(
      "Columns not found in SpatVector: %s",
      paste(missing_cols, collapse = ", ")
    )
  )

  # 3. Extract attributes for fast base-R comparison
  v_attr <- terra::values(v)

  # 4. Build and combine the logical masks
  mask_list <- lapply(filter_cols, function(col) {
    v_attr[[col]] == df[[col]][1]
  })

  final_mask <- Reduce(`&`, mask_list)

  # 5. Handle NAs to prevent subsetting errors
  final_mask[is.na(final_mask)] <- FALSE

  # 6. Subset and return
  return(v[final_mask, ])
}

single_wsp <- function(wsp_pattern) {
  out <- wsp_pattern |>
    vect()
  return(out)
}

distinct_grts_strata <- function(wsp_target) {
  out <- wsp_target |>
    st_as_sf() |>
    st_drop_geometry() |>
    distinct(grts_rank, stratum_name, changecat)
  return(out)
}

# helper functions to dissolve boundaries
# https://github.com/r-spatial/sf/issues/2422#issuecomment-2398718177
poly_2_nb_id <- function(x,
                         snap = NULL,
                         queen = TRUE,
                         quiet = TRUE,
                         ...) {
  rlang::check_installed("spdep")
  fn <- invisible
  if (quiet) {
    fn <- suppressWarnings
  }

  rlang::try_fetch(
    fn({
      nb <- spdep::poly2nb(x, snap = snap, queen = queen, ...)
      comp_nb <- spdep::n.comp.nb(nb)
      comp_nb[["comp.id"]]
    }),
    error = \(cnd) {
      rep_len(0, length(x))
    }
  )
}
st_dissolve_by <- function(x,
                           ...,
                           .by = NULL,
                           do_union = TRUE,
                           .data_key = "data",
                           .dissolve_key = "group.comp.id") {
  stopifnot(
    !rlang::has_name(x, .dissolve_key),
    is.data.frame(x)
  )

  # Handle tidyselect style .by arguments
  by <- rlang::enquo(.by)
  if (!dplyr::is_grouped_df(x) && !rlang::quo_is_null(by)) {
    x <- dplyr::group_by(x, dplyr::across(!!by))
    .by <- NULL
  }

  x_group_vars <- NULL

  if (dplyr::is_grouped_df(x)) {
    x_group_vars <- dplyr::group_vars(x)
    .by <- x_group_vars
    x <- dplyr::ungroup(x)
  }

  sf_column_nm <- attr(x, "sf_column")

  # Create dissolve key with poly_2_nb_id
  x <- x |>
    dplyr::mutate(
      "{.dissolve_key}" := paste0(
        dplyr::cur_group_id(), ".",
        poly_2_nb_id(.data[[sf_column_nm]], ...)
      ),
      .by = .by
    )

  # Use st_combine or st_union (if `do_union = TRUE`)
  sf_summarise_fn <- sf::st_combine
  if (do_union) {
    sf_summarise_fn <- sf::st_union
  }

  x_dissolve <- x |>
    dplyr::summarise(
      # Keep unique values for grouping variables (if supplied)
      dplyr::across(
        tidyselect::any_of(x_group_vars),
        unique
      ),
      # Combine geometry with sf summary function
      dplyr::across(
        tidyselect::all_of(sf_column_nm),
        sf_summarise_fn
      ),
      .by = tidyselect::all_of(.dissolve_key)
    )

  return(x_dissolve)
}


terra_dissolve_boundaries <- function(
    polygons,
    group_cols, # aggregate by
    mmu # minimum mappable unit
) {
  assertthat::assert_that(inherits(polygons, "SpatVector"))
  assertthat::assert_that(is.character(group_cols))
  assertthat::assert_that(
    all(group_cols %in% names(polygons)),
    msg = "Not all group_cols are present in the polygon attributes."
  )
  assertthat::assert_that(is.numeric(mmu), mmu >= 0)
  assertthat::assert_that(all(!is.na(terra::values(polygons[, group_cols]))))

  # Reject overlapping polygons input
  area_sum <- sum(terra::expanse(polygons))
  area_union <- sum(terra::expanse(terra::aggregate(polygons)))
  assertthat::assert_that(
    abs(area_sum - area_union) < 0.1,
    msg = "Input contains overlapping polygons, which are not allowed."
  )

  # Dissolve everything by final classification first
  aggregate_cols <- setdiff(names(polygons), group_cols)
  polygons_df <- values(polygons)
  polygons <- terra::aggregate(
    polygons[, group_cols], #drop non-grouping cols; handled separately
    by = group_cols,
    dissolve = TRUE,
    count = FALSE
  )
  # Summarize attributes with correct type handling
  aggregated_attribs <- polygons_df |>
    group_by(across(all_of(group_cols))) |>
    summarise(
      across(
        all_of(aggregate_cols),
        \(x) if (is.numeric(x)) {
          mean(x, na.rm = TRUE)
        } else {
          paste(unique(na.omit(x)), collapse = "|")
        }
      ),
      .groups = "drop"
    )

  # Join the summarized attributes back to the geometries
  polygons <- merge(polygons, aggregated_attribs, by = group_cols)

  polygons <- terra::makeValid(polygons)
  polygons <- terra::disagg(polygons)

  if (mmu > 0) {
    # Identify Slivers vs Keepers
    area_sqm <- terra::expanse(polygons)
    sliver_mask <- !is.na(area_sqm) & area_sqm < mmu

    if (any(sliver_mask)) {
      slivers <- polygons[sliver_mask, ]
      keepers <- polygons[!sliver_mask, ]

      # Edge case safety: Only do this if there is at least one polygon >= MMU
      if (nrow(keepers) > 0) {

        slv_df <- terra::values(slivers)

        # Buffer ALL slivers at once
        slivers_buf <- terra::buffer(slivers, width = 0.1)

        # Prevent column renaming during intersect
        terra::values(slivers_buf) <- NULL

        # Build a spatial index using relate.
        # pairs = TRUE returns a 2-column matrix of just the intersecting
        # indices (id.x = sliver, id.y = keeper)
        overlap_pairs <- terra::relate(
          slivers_buf, keepers, relation = "intersects", pairs = TRUE
        )

        # Iterate through slivers
        for (i in seq_len(nrow(slivers))) {

          # Find ONLY the keepers that touch this specific sliver
          candidate_idx <- overlap_pairs[overlap_pairs[, 1] == i, 2]

          if (length(candidate_idx) > 0) {

            # Subset down to the local candidates
            local_keepers <- keepers[candidate_idx, ]
            slv_buf_single <- slivers_buf[i, ]

            overlaps <- terra::intersect(local_keepers, slv_buf_single)

            if (nrow(overlaps) > 0) {
              overlap_areas <- terra::expanse(overlaps)
              best_idx <- which.max(overlap_areas)
              # Dynamically adopt the attributes
              slv_df[i, group_cols] <-
                terra::values(overlaps)[best_idx, group_cols]
            }
          }
        }

        # Put the updated attributes back into the slivers SpatVector
        terra::values(slivers) <- slv_df

        # Recombine the updated slivers with the keepers
        polygons <- rbind(keepers, slivers)

        # Final Dissolve
        aggregate_cols <- setdiff(names(polygons), group_cols)
        polygons_df <- values(polygons)
        polygons <- terra::aggregate(
          polygons[, group_cols],
          by = group_cols,
          dissolve = TRUE,
          count = FALSE
        )
        # Summarize attributes with correct type handling
        aggregated_attribs <- polygons_df |>
          group_by(across(all_of(group_cols))) |>
          summarise(
            across(
              all_of(aggregate_cols),
              \(x) if (is.numeric(x)) {
                mean(x, na.rm = TRUE)
              } else {
                paste(unique(na.omit(x)), collapse = "|")
              }
            ),
            .groups = "drop"
          )

        # Join the summarized attributes back to the geometries
        polygons <- merge(polygons, aggregated_attribs, by = group_cols)
        polygons <- terra::makeValid(polygons)
        polygons <- terra::disagg(polygons)
      }
    }
  }
  return(polygons)
}

intersect_validation_polygons <- function(
    wsp_target, lu_changecat, input_years, mmu = 30) {
  out <- wsp_target |>
    st_as_sf() |>
    filter(stratum_name == lu_changecat) |>
    mutate(
      year_flea2 = year_flea
    ) |>
    st_make_valid() |>
    group_by(grts_rank, stratum_name, changecat, year_flea2) |>
    tidyr::nest()

  if (nrow(out) == 0) return(vect())

  out <- out |>
    mutate(
      data = purrr::map(data, function(x) {
        x %>% rename_with(
          .cols = c(layer, year_grb, value, label, prop_in_settlement),
          .fn = \(x) paste0(x, "_", .$year_flea[1], recycle0 = TRUE)
        )
      }
      )
    ) |>
    tidyr::pivot_wider(
      id_cols = c(grts_rank, stratum_name, changecat),
      names_from = year_flea2,
      names_prefix = "data_",
      values_from = data
    ) |>
    tidyr::pivot_longer(
      cols = starts_with("data_"),
      names_to = "year_flea",
      names_prefix = "data_") |>
    mutate(
      value = lapply(value, vect),
      value = lapply(value, makeValid)
    )

  out <- out |>
    summarize(
      intersected_data = list(
        purrr::reduce(value, terra::intersect)),
      .groups = "drop"
    ) |>
    mutate(
      intersected_data = lapply(intersected_data, st_as_sf)
    ) |>
    tidyr::unnest(intersected_data) |>
    st_as_sf(crs = 31370) |>
    select(!starts_with("year_flea")) |>
    vect() |>
    unique()

  # extra columns in case of doubt over label (second choice label)
  new_cols <- setNames(
    rep(list("none"), length(input_years)),
    paste0("label_", input_years, "_2")
  )

  out <- st_as_sf(out) |>
    rowwise() |>
    mutate(
      labels = paste(
        c_across(all_of(paste0("label_", input_years))),
        collapse = "-")
    ) |>
    ungroup() |>
    mutate(!!!new_cols) |>
    st_make_valid() |>
    st_cast("GEOMETRYCOLLECTION") |>
    st_collection_extract(type = "POLYGON") |>
    vect()
  out <- .filter_empty_geom(out)

  # --------------------------------------------------------------------
  # MMU Enforcement: Longest Shared Border Assignment
  # --------------------------------------------------------------------

  group_cols <- c(
    "grts_rank", "stratum_name", "changecat", "labels",
    paste0("label_", input_years),
    paste0("label_", input_years, "_2")
  )
  # need to isolate each validation polygons
  # to avoid slivers from one stratum get dissolved into another stratum
  stratum_cols <- c("grts_rank", "stratum_name", "changecat")

  # Get unique strata
  strata_df <- unique(values(out)[, stratum_cols])

  # Apply dissolve function to each stratum
  processed_list <- vector(mode = "list", length = nrow(strata_df))
  for (strat in seq_len(nrow(strata_df))) {
    print(sprintf("Stratum number %s out of %s", strat, nrow(strata_df)))
    # Subset to just this stratum
    sub_poly <- filter_spatvector(out, strata_df[strat, ])
    # Run dissolve function
    # terra_dissolve_boundaries <- purrr::possibly(
    #   terra_dissolve_boundaries, otherwise =  sub_poly
    # )
    processed_list[[strat]] <- terra_dissolve_boundaries(
      polygons = sub_poly,
      group_cols = group_cols,
      mmu = mmu
    )
  }
  # Combine the results back into a single SpatVector
  out_combined <- do.call(rbind, processed_list)
  out_combined <- terra::unique(out_combined)

  return(out_combined)
}

crop_labeled_polygons <- function(
  pvp,
  crop_with
) {
  spsub <- pvp[crop_with]
  spsub <- st_as_sf(spsub) |> st_make_valid()
  cropped <- st_intersection(spsub, st_as_sf(crop_with))
  cropped <- cropped |>
    filter(grts_rank == grts_rank.1, changecat == changecat.1,
           stratum_name == stratum_name.1) |>
    select(-ends_with(".1")) %>%
    filter(st_geometry_type(.) %in% c("POLYGON", "MULTIPOLYGON")) |>
    st_cast("MULTIPOLYGON") |>
    st_cast("POLYGON") |>
    vect()

  # removed code that cleaned up small intersections
  # the pvp object is already cleaned (geoms < mmu removed)
  # any geoms < mmu resulting from cropping the 90x90 to 50x50 are kept

  return(cropped)
}

# helper function to combine small polygons with larger polygons
# based on shared boundary length
combineGeoms_ <- purrr::possibly(terra::combineGeoms)
combine_small_intersections <- function(
    spatvec,
    mmu = 10) {

  assertthat::assert_that(is.numeric(mmu) && mmu > 0)

  # 1. Standardize and cleanup
  spatvec <- terra::snap(spatvec, tolerance = 0.1)
  spatvec <- terra::makeValid(spatvec)
  spatvec <- disagg(spatvec)

  # Calculate area using a temporary safe column name
  spatvec$tmp_area_calc <- expanse(spatvec)

  # Return early if clean
  if (all(spatvec$tmp_area_calc >= mmu)) {
    spatvec$tmp_area_calc <- NULL
    return(spatvec)
  }

  # 2. Split Data
  small_spatvec <- spatvec[spatvec$tmp_area_calc < mmu]
  large_spatvec <- spatvec[spatvec$tmp_area_calc >= mmu]

  # Clean up the area column
  large_spatvec$tmp_area_calc <- NULL
  small_spatvec$tmp_area_calc <- NULL

  # 3. Try combineGeoms
  small_for_combine <- small_spatvec
  values(small_for_combine) <- NULL

  new <- combineGeoms_(
    large_spatvec, small_for_combine,
    overlap = TRUE, boundary = TRUE, distance = FALSE,
    append = FALSE, dissolve = TRUE, erase = TRUE
  )

  # 4. Fallback Strategy: Assign & Aggregate
  if (is.null(new)) {
    message("combineGeoms failed. Switching to 'Assign and Aggregate' fallback.")

    # A. Find the nearest large polygon for each small polygon
    nearest_neighbors <- terra::nearby(
      small_spatvec, large_spatvec,
      centroids = FALSE
    )
    nearest_neighbors <- large_spatvec[nearest_neighbors[, "k1"]]

    # B. Transfer attributes
    # We overwrite the small polygon attributes with their large neighbor's
    # attributes
    # This effectively "merges" them logically, even if the border isn't
    # dissolved yet
    if (nrow(small_spatvec) == nrow(nearest_neighbors)) {
      values(small_spatvec) <- values(nearest_neighbors)
    } else {
      message(
        paste0("Nearest neighbor failed.",
        " Returning only polygons for which area larger than threshold.")
        )
      return(large_spatvec)
    }

    # C. Bind them back together
    # We now have a list of polygons where the small ones look identical
    # (data-wise) to their neighbors
    combined <- rbind(large_spatvec, small_spatvec)

    # D. Aggregate (Dissolve)
    # This melts the boundary lines between polygons that share the same
    # attributes.
    # It automatically handles the "overlap" issue because they become one
    # geometry.
    all_cols <- names(combined)
    label_cols <- all_cols[grep("^label", all_cols)]
    final <- terra::aggregate(combined, by = label_cols, dissolve = TRUE)
    final$mean_agg_n <- NULL
    names(final) <- all_cols

    # E. Final Cleanup
    final <- disagg(final) # cast multipolygon to polygon
    final <- terra::makeValid(final)

    return(final)

  } else {
    new <- terra::makeValid(new)
    return(new)
  }
}


unnest_spatvector <- function(df, col_name) {

  # 1. Extract the list of SpatVectors
  sv_list <- df[[col_name]]

  # 2. Calculate how many polygons are in each list element
  n_polys <- vapply(sv_list, function(x) {
    if (is.null(x)) 0 else nrow(x)
  }, FUN.VALUE = numeric(1))

  # 3. Combine all SpatVectors into one 'long' SpatVector
  combined_sv <- vect(sv_list)

  # 4. Expand the original dataframe metadata
  # We repeat each row of the dataframe 'n_polys' times
  df_expanded <- df[rep(1:nrow(df), n_polys), ]

  # 5. Remove the heavy list-columns from the metadata
  cols_to_remove <- names(df_expanded)[
    purrr::map_chr(df_expanded, class) == "list"
  ]
  df_expanded <- df_expanded[
    , !names(df_expanded) %in% cols_to_remove,
    drop = FALSE]

  # 6. Merge metadata with SpatVector attributes
  values(combined_sv) <- cbind(values(combined_sv), df_expanded)

  return(combined_sv)
}

raster_to_polygons <- function(
    tm, # temporal raster
    vp, # validation polygons
    input_years
) {
  assertthat::assert_that(inherits(tm, "SpatRaster"))
  assertthat::assert_that(inherits(vp, "SpatVector"))
  assertthat::assert_that(is.character(input_years))

  # mask raster using validation polygons
  masked <- mask(tm, vp)
  # first as.points, directly as.polygons runs into memory issues
  masked_p <- as.points(masked, values = TRUE, na.rm = TRUE)

  # add polygon attributes back ensuring every validation poly is joined with
  # the raster points
  masked_pp <- masked_p |>
    st_as_sf() |>
    st_join(
      st_as_sf(vp),
      join = st_intersects,
      left = TRUE
    ) |>
    vect()

  # convert points to squares 10x10
  masked_v <- point_to_gridcell(masked_pp, cell_width_m = 10)

  # extra columns in case of doubt over label (second choice label)
  new_cols <- setNames(
    rep(list("none"), length(input_years)),
    paste0("label_", input_years, "_2")
  )

  masked_f <- masked_v |>
    st_as_sf() |>
    mutate(
      label = gsub("_", "-", label)
    ) |>
    tidyr::separate_wider_delim(
      label,
      "-",
      names = input_years,
      names_sep = "_",
      cols_remove = FALSE
    ) |>
    mutate(!!!new_cols) |>
    st_as_sf() |>
    vect()
  return(masked_f)
}


safe_cover <- function(x, y, ...) {
  common <- intersect(names(x), names(y))
  x_types <- sapply(values(x), function(col) paste(class(col), collapse = "_"))
  y_types <- sapply(values(y), function(col) paste(class(col), collapse = "_"))
  x_check <- x_types[names(x_types) %in% common][common]
  y_check <- y_types[names(y_types) %in% common][common]
  if (!identical(x_check, y_check)) {
    stop(
      sprintf(
        "Type mismatch found for variables %s",
        paste(names(x_check)[x_check != y_check], collapse = ", ")
      )
    )
  }

  cover(x, y, ...)
}
