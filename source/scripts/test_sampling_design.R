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
    point_position = c("center", "lowerleft", "upperleft", "lowerright", "upperright"),
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





library(terra)
library(sf)
library(dplyr)
library(mapview)

flea_data <- gsub(
  pattern = "flea-extent", replacement = "flea-data", x = here::here())

lg2013 <- rast(file.path(flea_data, "data", "2013", "LG2013_finaal_update.tif"))
lg2016 <- rast(file.path(flea_data, "data", "2016", "LG2016_finaal_update.tif"))
lg2019 <- rast(file.path(flea_data, "data", "2019", "LG2019_finaal_update.tif"))

qgisprocess::qgis_show_help("slyr:lyrtoqml")

qgisprocess::qgis_run_algorithm(
  "slyr:lyrtoqml",
  INPUT = file.path(flea_data, "data", "2013", "LG2013_finaal_update.lyr"),
  OUTPUT = file.path(flea_data, "data", "2013", "LG2013_finaal_update.qml"))

slyr <- xml2::read_xml(
  file.path(flea_data, "data", "2013", "LG2013_finaal_update.qml"))

catstable <- xml2::xml_find_all(x = slyr, ".//pipe//rasterrenderer//colorPalette") |>
  xml2::xml_contents() |>
  purrr::map(xml2::xml_attrs) |>
  purrr::map_df(~as.list(.)) |>
  dplyr::relocate(value, label) |>
  dplyr::mutate(
    value = as.numeric(value),
    color = toupper(color))

apply_cats <- function(x, cats = catstable, name) {
  xc <- as.factor(x)
  names(cats)[names(cats) == "label"] <- name
  levels(xc) <- cats
  coltab(xc) <- cats |> dplyr::select(value, color) |> as.data.frame()
  names(xc) <- name
  return(xc)
}

lg2013 <- apply_cats(lg2013, name = "lg2013")
lg2016 <- apply_cats(lg2016, name = "lg2016")
lg2019 <- apply_cats(lg2019, name = "lg2019")


lg2013_selectie <- crop(lg2013, ext(200000, 205000, 200000, 205000))
lg2016_selectie <- crop(lg2016, ext(200000, 205000, 200000, 205000))
lg2019_selectie <- crop(lg2019, ext(200000, 205000, 200000, 205000))


# calculate for each pixel the dominant land cover inside a 9x9 block centered
# on the focal pixel
lg2013_stratification <- focal(lg2013_selectie, w = 9, fun = "modal") |>
  apply_cats(name = "lg2013")
lg2016_stratification <- focal(lg2016_selectie, w = 9, fun = "modal") |>
  apply_cats(name = "lg2016")
lg2019_stratification <- focal(lg2019_selectie, w = 9, fun = "modal") |>
  apply_cats(name = "lg2019")

rasterstoplot <- c(lg2013_selectie, lg2013_stratification)
names(rasterstoplot) <- c("2013 original", "2013 stratification")
plot(rasterstoplot)

# combine the spatial strata into temporal strata

temporal_stratification <- purrr::reduce(
  list(lg2013_stratification, lg2016_stratification, lg2019_stratification),
  concats)

binary_change <- function(data, lg) {
  binary <- vector("list", length = length(lg))
  binary <- setNames(binary, lg)
  for (i in lg) {
    binary[[i]] <- paste0(
      stringr::str_detect(data$lg2013, i) %>% as.numeric(),
      stringr::str_detect(data$lg2016, i) %>% as.numeric(),
      stringr::str_detect(data$lg2019, i) %>% as.numeric())
  }
  bind_cols(data, binary)
}

lg <- gsub(pattern = "^\\d\\s-\\s", replacement = "", x = catstable$label)

additional_levels <- freq(temporal_stratification) %>%
  as_tibble() %>%
  tidyr::separate(
    value, into = c("lg2013", "lg2016", "lg2019"),
    sep = "_", remove = FALSE) %>%
  binary_change(lg = lg) %>%
  rowwise() %>%
  mutate(stable = ifelse(
    all(lg2013 == lg2016, lg2016 == lg2019),
    "stable", "changed") %>%
      as.factor()) %>%
  ungroup()

join_levels <- cats(temporal_stratification)[[1]] %>%
  inner_join(
    additional_levels,
    by = join_by("lg2013_lg2016_lg2019" == "value"))
levels(temporal_stratification) <- join_levels
coltab(temporal_stratification) <- NULL


# map showing stable (TRUE) vs changed (FALSE)
activeCat(temporal_stratification) <- "stable"
plot(temporal_stratification)
additional_levels %>%
  group_by(stable) %>%
  summarise(n_pixels = sum(count))

binary_levels <- data.frame(
  value = c(0,1,10,11,100,101,110,111)) %>%
  mutate(
    label = stringr::str_pad(value, side = "left", pad = "0", width = 3),
    label2 = c(
      "not present", "gained", "dynamic",
      "gained", "lost", "dynamic", "lost", "stable"))


ts <- terra::catalyze(temporal_stratification) %>%
  terra::subset(subset = lg) %>%
  sapp(fun = function(x, cats = binary_levels) {
    nm <- names(x)
    xc <- as.factor(x)
    names(cats)[names(cats) == "label"] <- nm
    levels(xc) <- cats
    names(xc) <- nm
    return(xc)
  })

plot(ts)
ts2 <- ts %>%
  sapp(fun = function(x) {
    nm <- names(x)
    activeCat(x) <- "label2"
    names(x) <- nm
    return(x)
  })
plot(ts2)


# lazy conversion to points
lg2013_strat_points <- as.points(lg2013_stratification)

# conversion to sf brings into memory (but see lazysf)
lg2013_strat_points_sf <- st_as_sf(lg2013_strat_points)
lg2013_strat_points_df <- lg2013_strat_points_sf %>%
  bind_cols(st_coordinates(.)) %>%
  st_drop_geometry()

library(SamplingBigData)
sample_size <- 100
minimum_n_strat <- 5
lg2013_sample_strat <- vector(
  "list", length = length(unique(lg2013_strat_points_df$lg2013)))
lg2013_sample_strat <- setNames(
  lg2013_sample_strat, nm = unique(lg2013_strat_points_df$lg2013))
set.seed(214)
for (i in unique(lg2013_strat_points_df$lg2013)) {
  size_pop <- nrow(lg2013_strat_points_df)
  df <- lg2013_strat_points_df %>%
    filter(lg2013 == i)
  size_strat <- nrow(df)
  n_strat <- max(minimum_n_strat, round(size_strat / size_pop * sample_size))
  ips <- rep(n_strat / size_strat, size_strat)
  rowindex <- lpm2_kdtree(
    prob = ips,
    x = df[,c("X", "Y")],
    inOrder = TRUE)
  lg2013_sample_strat[[i]] <- df %>%
    slice(rowindex) %>%
    mutate(order = 1:n())
}

lg2013_sample_strat <- bind_rows(lg2013_sample_strat)

lg2013_sample_strat_sf <- st_as_sf(
  lg2013_sample_strat, coords = c("X", "Y"),
  crs = st_crs(lg2013_strat_points_sf))

lg2013_sample_strat_vect <- vect(lg2013_sample_strat_sf)

testje <- terra::crop(
  lg2013_stratification,
  lg2013_sample_strat_vect[1] %>% buffer(5))
plot(testje)
points(lg2013_sample_strat_vect[1])
polys(lg2013_sample_strat_vect[1] %>% buffer(5))

lg2013_sample_strat_sf_block <- lg2013_sample_strat_sf %>%
  point_to_gridcell(cell_width_m = 90, crs = 31370) %>%
  mutate(lg2013 = factor(lg2013, levels = catstable$label))

testje2 <- terra::crop(
  lg2013_stratification,
  lg2013_sample_strat_vect[1] %>% buffer(45 + 10))
plot(testje2)
points(lg2013_sample_strat_vect[1])
polys(lg2013_sample_strat_vect[1] %>% buffer(45))
polys(vect(lg2013_sample_strat_sf_block)[1])



plot(lg2013_stratification)
points(lg2013_sample_strat_vect)


write_sf(
  lg2013_sample_strat_sf_block,
  file.path(flea_data, "data", "2013", "lg2013_sample_strat_sf_block.gpkg"),
  delete_dsn = TRUE)

palette_inbo <- leaflet::colorFactor(
  palette = catstable$color, levels = catstable$label)


mapview(lg2013_selectie, alpha.regions = 0.3, maxpixels = 1e6) +
  mapview(
    lg2013_sample_strat_sf_block,
    zcol = "lg2013",
    alpha.regions = 0.5,
    color = palette_inbo(catstable$label),
    col.regions = palette_inbo(catstable$label))


