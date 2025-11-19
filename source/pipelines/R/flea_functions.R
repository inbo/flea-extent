#' Convert a single point location to a grid cell polygon
#'
#' @param xy a SpatVector with geometry type points
#' @param cell_width_m cell width in meter, default 500
#' @param point_position default centre of grid cell
#' @param crs default EPSG code 31370
#'
#' @return a SpatVector with geometry type polygon
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
#' @return SpatVector object containing the extracted points with their values,
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
#' @return A SpatRaster
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
#' @return A SpatVector
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



#' Get data from the WFS service for GRB
#'
#' @param layer A string. Should be one of the WFS layers available in the
#' service
#' @param bbox A SpatExtent or an object from which a SpatExtent can be
#' determined. The bbox values should be in CRS 31370.
#'
#' @return A SpatVector containing GRB layer objects that intersect bbox.
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
      value = NA
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
        layer %in% c("GRB:GBG", "GRB:GBA") ~ 101,
        layer %in% c(
          "GRB:WBN",
          "GRB:SBN",
          "GRB:KNW") |
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


#' Helper functions to deal with empty records in tar_terra_vect
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
  # catch case empty records
  # can be removed when terra or geotargets deals with this natively
  # https://github.com/ropensci/geotargets/issues/187
  if (nrow(x) == 0) {
    out <- .create_empty_geom(x, type = terra::geomtype(y))
    return(out)
  }
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
  ws$value <- NA
  ws$area_name <- NULL
  ws$wfd_type_certain <- NULL

  return(ws)
}


combine_grb_inbo_water <- function(grb_water, inbo_water, meta) {
  inbo_water <- vect(inbo_water)
  inbo_water <- inbo_water[
    inbo_water$layer == paste0("watersurfaces_",meta$version), ]

  # cover: values of x that overlap with y are replaced by y
  grb_water <- aggregate(x = grb_water, by = names(grb_water))
  water <- cover(x = grb_water, y = inbo_water)
  water <- unique(water)
  return(water)
}

combine_water_settlements <- function(
    water, settlements, polygons,
    lbg_101, lbg_104) {

  assertthat::assert_that(inherits(water, "SpatVector")) # a branch
  assertthat::assert_that(inherits(settlements, "SpatVector")) # a target
  assertthat::assert_that(inherits(lbg_101, "list")) # a pattern
  assertthat::assert_that(inherits(lbg_104, "list")) # a pattern
  assertthat::assert_that(inherits(polygons, "list")) # a pattern

  vp <- vect(polygons)
  lbg_101 <- vect(lbg_101) |> unique() # this combines multiple years
  lbg_104 <- vect(lbg_104) |> unique() # this combines multiple years

  # get the validation year
  year_to_validate <- unique(water$year_flea)
  year_to_validate <- year_to_validate[!is.na(year_to_validate)]

  # filter the lbg layers to only the validation year
  lbg_101 <- lbg_101[grepl(year_to_validate, x = lbg_101$layer), ]
  lbg_104 <- lbg_104[grepl(year_to_validate, x = lbg_104$layer), ]

  #vp <- vp[1:200,] # testing only
  vplist <- vector("list", nrow(vp))
  for (i in seq_along(vp)) {
    print(sprintf("%s out of %s done", i, nrow(vp)))
    vp_ <- vp[i] # selecteert 1 validatie-polygoon
    w_ <- water[vp_]
    w_area <- expanse(w_)
    w_ <- subset(w_, w_area > 1)
    w_ <- crop(w_, vp_)
    w_ <- subset(w_, w_$grts_rank == vp_$grts_rank | is.na(w_$grts_rank))
    s_ <- settlements[vp_]
    s_ <- crop(s_, vp_)
    s_ <- s_[s_$grts_rank == vp_$grts_rank, ]
    lbg_101_ <- lbg_101[vp_]
    lbg_101_ <- crop(lbg_101_, vp_)
    lbg_104_ <- lbg_104[vp_]
    lbg_104_ <- crop(lbg_104_, vp_)
    lbg_ <- rbind(lbg_101_, lbg_104_)
    c1 <- cover(s_, lbg_)
    c1_area <- expanse(c1)
    c1 <- subset(c1, c1_area > 1)
    c2 <- cover(w_, c1)
    c2_area <- expanse(c2)
    c2 <- subset(c2, c2_area > 1)
    # verwijder polygonen waar layer GRB en jaar (van GRB) > year_to_validate
    c2 <- subset(c2, !(c2$jaar > year_to_validate & grepl("^GRB", c2$layer)))
    if (nrow(c2) > 1) {
      c2 <- unique(c2)
      # possibly overlapping polygons (with slight differences in geom)
      c2 <- aggregate(
        x = c2,
        by = names(c2),
        dissolve = TRUE
      )
    }
    out <- cover(vp_, c2)
    # Repeat until areas match within tolerance
    # give up after 4 attempts
    j <- 1
    while (abs(expanse(vp_) - sum(expanse(out))) > 1) {
      out <- cover(vp_, out)
      j <- j + 1
      if (j == 4) {
        break
      }
    }
    out_area <- expanse(out)
    out <- subset(out, out_area > 1)
    out$grts_rank <- vp_$grts_rank
    out$cell <- vp_$cell
    out$stratum_name <- vp_$stratum_name
    out$changecat <- vp_$changecat
    out$year_flea <- year_to_validate
    out <- disagg(out) #casts multipolygon to polygon
    vplist[[i]] <- out
  }
  # # Ensure all SpatVectors have the same attributes
  # common_cols <- Reduce(intersect, lapply(vplist, names))
  #
  # vplist_std <- lapply(vplist, function(x) {
  #   x[, common_cols]
  # })

  all_types <- data.frame(
    colname = c(
      "grts_rank", "cell", "stratum_name",
      "changecat", "gml_id", "layer", "jaar", "lbl", "value", "year_flea",
      "agg_n", "polygon_id", "wfd_code", "hyla_code", "name", "wfd_type",
      "depth_class", "connectivity", "usage", "wfd_type_alternative",
      "water_level_management", "GWSCOD_H", "GWSNAM_H"),
    type = c(
      "numeric",
      "numeric", "character", "character", "character", "character",
      "numeric", "character", "numeric", "numeric", "numeric", "character",
      "character", "numeric", "character", "character", "character",
      "character", "character", "character", "character", "character",
      "character")
  )


  all_cols <- all_types$colname

  vplist_complete <- lapply(vplist, function(x) {
    missing_cols <- setdiff(all_cols, names(x))
    if (length(missing_cols) > 0) {
      # Add missing columns with NA of appropriate type
      for (col in missing_cols) {
        target_type <- all_types$type[all_types$colname == col]
        # Create NA of the correct type
        x[[col]] <- as(NA, target_type)
      }
    }
    x <- x[, all_cols]  # Reorder to match
    for (col in all_cols) {
      target_type <- all_types$type[all_types$colname == col]
      values(x)[[col]] <- as(values(x)[[col]], target_type)
    }
    return(x)
  })


  vp_wa_se <- vect(vplist_complete)

  # make sure records are unique
  vp_wa_se <- terra::unique(vp_wa_se) |> terra::disagg()

  # remove tiny areas
  areas <- expanse(vp_wa_se)
  vp_wa_se <- subset(vp_wa_se, areas > 1)

  return(vp_wa_se)
}

postprocess_water_settlements <- function(water_settlements) {
  ws <- st_as_sf(water_settlements)

  # get the validation year
  year_to_validate <- unique(ws$year_flea)
  year_to_validate <- year_to_validate[!is.na(year_to_validate)]

  ws <- ws |>
    select(
      grts_rank,
      layer,
      year_grb = jaar,
      stratum_name,
      changecat,
      value,
      year_flea
    ) |>
    group_by(grts_rank) |>
    mutate(
      stratum_name = ifelse(
        is.na(stratum_name), first(stratum_name), stratum_name),
      changecat = ifelse(
        is.na(changecat), first(changecat), changecat),
      label = case_when(
        value == 101 ~ "101",
        value == 102 ~ "102",
        value == 104 ~ "104",
        layer %in% c("GRB:GBA", "GRB:GBG") ~ "101",
        layer %in% c(
          "GRB:WBN",
          "GRB:SBN",
          "GRB:KNW",
          "GRB:TRN") ~ "102",
        grepl("Landbouwgebruikspercelen", layer) ~ "101",
        grepl("watersurfaces", layer) ~ "water",
        layer == "GRB:WTZ" ~ "water",
        TRUE ~ "other"
      ),
      year_flea = year_to_validate
    )
  ws <- vect(ws) |> unique()
  return(ws)
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




intersect_validation_polygons <- function(
    wsp_target, lu_changecat) {
  out <- wsp_target |>
    st_as_sf() |>
    filter(stratum_name == lu_changecat) |>
    mutate(year_flea2 = year_flea) |>
    st_make_valid() |>
    group_by(grts_rank, stratum_name, changecat, year_flea2) |>
    tidyr::nest()

  if (nrow(out) == 0) return(vect())

  out <- out |>
    mutate(
      data = purrr::map(data, function(x) {
        x %>% rename_with(
          .cols = c(layer, year_grb, value, label),
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

  out_agg <- out |>
    st_as_sf() |>
    mutate(
      labels = paste(label_2016, label_2019, label_2022, sep = "-")
    ) |>
    st_dissolve_by(
      .by = c(
        grts_rank, stratum_name, changecat, labels
        # , layer_2016, year_grb_2016, value_2016,
        # label_2016, layer_2019, year_grb_2019, value_2019,
        # label_2019, layer_2022, year_grb_2022,
        # value_2022, label_2022
      )
    ) |>
    vect() |>
    .filter_empty_geom()

  return(out_agg)
}

