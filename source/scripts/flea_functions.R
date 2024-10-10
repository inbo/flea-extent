#' Convert a single point location to a grid cell polygon
#'
#' @param xy an object of class POINT
#' @param cell_width_m cell width in meter, default 500
#' @param point_position default center of grid cell
#' @param crs default EPSG code 31370
#'
#' @return
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

  stopifnot(sf::st_is(xy, "POINT"))
  xy_df <- sf::st_drop_geometry(xy)
  xy <- sf::st_geometry(xy)

  # buffer with 1 point per quandrant
  halflength <- cell_width_m / 2
  xy_buffer <- sf::st_buffer(x = xy,
                             dist = sqrt(2 * halflength^2),
                             nQuadSegs = 1)

  # rotate 45 degrees around centroid
  rot <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
  pl <- (xy_buffer - xy) * rot(pi/4) + xy
  pl <- sf::st_sf(data.frame(xy_df, pl), crs = crs)
  return(pl)
}

#' Extract a sample of points from a raster
#'
#' This function extracts a specified number of points from a raster,
#' selecting the lowest values after sorting.
#'
#' @param rast A raster object (terra SpatRaster) containing GRTS rank numbers
#' @param n The number of points to extract
#'
#' @return An sf object containing the extracted points with their values,
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
  require(terra)
  require(sf)

  # Extract values, exclude NA
  extracted <- terra::extract(
    x = rast,
    y = as.points(rast, na.rm = TRUE),
    cells = TRUE,
    xy = TRUE,
    ID = FALSE)

  # Sort and select the lowest n
  sorted_indices <- order(extracted[[1]])[1:n]
  selected <- extracted[sorted_indices, ]
  selected$stratum_name <- names(selected)[1]
  names(selected)[1] <- "grts_rank"

  # Convert to sf
  selected <- selected |>
    st_as_sf(coords = c("x", "y"), crs = 31370)

  return(selected)
}


#' Extract stratified sample from a raster
#'
#' This function extracts a stratified sample from a raster based on temporal
#' stratification and GRTS (Generalized Random Tessellation Stratified)
#' rankings.
#'
#' @param stratum_raster A raster object representing the stratification
#' @param fleagrts A raster object containing GRTS rankings
#' @param stratum_name Character string specifying the name of the stratum
#' @param ntot Integer specifying the total number of samples to extract
#' @param nmin Integer specifying the minimum number of samples per stratum
#' @param min_stratum_size Integer specifying the minimum stratum size in terms
#' of number of raster cells
#'
#' @return A list of sf objects, each containing extracted samples for a change
#' category
#'
#' @importFrom terra `activeCat<-` classify segregate crop global
#' @importFrom dplyr %>% mutate
#' @importFrom tibble as_tibble
#'
#' @export
extract_sample <- function(
    stratum_raster,
    fleagrts,
    stratum_name,
    ntot,
    nmin,
    min_stratum_size) {
  require(terra)
  require(dplyr)
  require(tibble)

  # assertions
  assertthat::assert_that(!missing(stratum_raster))
  assertthat::assert_that(inherits(stratum_raster, "SpatRaster"))
  assertthat::assert_that(!missing(fleagrts))
  assertthat::assert_that(inherits(fleagrts, "SpatRaster"))
  assertthat::assert_that(terra::compareGeom(stratum_raster, fleagrts))
  assertthat::assert_that(assertthat::is.string(stratum_name))
  assertthat::assert_that(is.numeric(ntot), ntot > 0)
  assertthat::assert_that(is.numeric(nmin), nmin < ntot)

  cat(stratum_name)

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
    changecat = 1:4,
    label = levels(as.factor(levelvec))
  )

  # make separate layer for each change cat
  ts2 <- ts2 %>% segregate(classes = NULL, other = NA)
  names(ts2) <- levels(as.factor(levelvec))

  # make separate layers containing grts rankings for each change cat
  fleagrts_ts2 <- crop(fleagrts, ts2, mask = TRUE)
  names(fleagrts_ts2) <- levels(as.factor(levelvec))

  # determine stratum population sizes
  popsize <- global(fleagrts_ts2, fun = "notNA") |>
    as_tibble(rownames = "layername")

  # check for strata that are too small
  remove_me <- popsize$layername[popsize$notNA < min_stratum_size]
  if (length(remove_me) > 0) {
    cat(
      sprintf(
        "Change category %s of %s is too small and will be removed",
        remove_me,
        stratum_name
      )
    )
    popsize <- popsize[!(popsize$layername %in% remove_me), ]
    fleagrts_ts2 <- terra::subset(fleagrts_ts2, remove_me, negate = TRUE)
  }

  # determine sample size allocation
  # first distribute nmin to each stratum,
  # remaining allocate proportional to stratum size
  allocation <- popsize |>
    mutate(
      n_h = nmin + round((ntot - nmin * n()) * (notNA / sum(notNA)))
    )

  sample_ts2 <- vector(mode = "list", length = nlyr(fleagrts_ts2))
  sample_ts2 <- setNames(sample_ts2, names(fleagrts_ts2))
  for (i in names(sample_ts2)) {
    sample_ts2[[i]] <- extract_sample_helper(
      fleagrts_ts2[[i]],
      allocation$n_h[allocation$layername == i]
    )
  }

  return(sample_ts2)
}


