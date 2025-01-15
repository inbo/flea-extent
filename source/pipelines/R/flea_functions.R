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
#' @return A SpatVector
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

  return(grb)
}
