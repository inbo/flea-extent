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
library(ggplot2)
library(ggsankey)

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

apply_cats <- function(x, cats = catstable, name, coltab = TRUE) {
  xc <- as.factor(x)
  names(cats)[names(cats) == "label"] <- name
  levels(xc) <- cats
  if (coltab) {
    coltab(xc) <- cats |> dplyr::select(value, color) |> as.data.frame()
  }
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
##################################################

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

# all transitions

df <- join_levels %>%
  ggsankey::make_long(
    lg2013, lg2016, lg2019,
    value = count
  )

df2 <-  df %>%
  group_by(x, node) %>%
  summarise(n = sum(value))

df3 <- df %>%
  left_join(df2)

p <- df3 %>%
  ggplot(aes(x = x,
             next_x = next_x,
             node = node,
             next_node = next_node,
             fill = factor(node),
             label = paste0(node,": n = ", n),
             value = value)) +
  geom_sankey(alpha = 0.5) +
  geom_sankey_label(alpha = 0.5, colour = "black") +
  theme_sankey() +
  theme(legend.position = "none")

p

# only changes

df <- join_levels %>%
  filter(stable == "changed") %>%
  ggsankey::make_long(
    lg2013, lg2016, lg2019,
    value = count
  )

df2 <-  df %>%
  group_by(x, node) %>%
  summarise(n = sum(value))

df3 <- df %>%
  left_join(df2)

p <- df3 %>%
  ggplot(aes(x = x,
             next_x = next_x,
             node = node,
             next_node = next_node,
             fill = factor(node),
             label = paste0(node,": n = ", n),
             value = value)) +
  geom_sankey(alpha = 0.5) +
  geom_sankey_label(alpha = 0.5, colour = "black") +
  theme_sankey() +
  theme(legend.position = "none")

p



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
      "gained", "lost", "dynamic", "lost", "stable"),
    col = c("lightgrey", "lightblue", "orange1", "darkblue", "darkred",
            "orange2", "red1",
             "green"),
    col2 = c("lightgrey", "blue", "orange", "blue", "red", "orange", "red",
            "green"))

ts <- terra::catalyze(temporal_stratification) %>%
  terra::subset(subset = lg) %>%
  sapp(fun = function(x, cats = binary_levels) {
    nm <- names(x)
    xc <- as.factor(x)
    names(cats)[names(cats) == "label"] <- nm
    levels(xc) <- cats
    names(xc) <- nm
    coltab(xc) <- cats[, c("value", "col")]
    return(xc)
  })

plot(ts)
ts2 <- ts %>%
  sapp(fun = function(x, cats = binary_levels) {
    nm <- names(x)
    activeCat(x) <- "label2"
    names(x) <- nm
    coltab(x) <- cats %>% select(value, col = col2)
    return(x)
  })
plot(ts2)

# repeat above steps with following differences:
# create difference maps without applying focal filter on status maps
# apply focal filter afterwards
#####################################################


temporal_maps <- purrr::reduce(
  list(lg2013_selectie, lg2016_selectie, lg2019_selectie),
  concats)

tm_additional_levels <- freq(temporal_maps) %>%
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

tm_join_levels <- cats(temporal_maps)[[1]] %>%
  inner_join(
    tm_additional_levels,
    by = join_by("lg2013_lg2016_lg2019" == "value"))
levels(temporal_maps) <- tm_join_levels
coltab(temporal_maps) <- NULL

# all transitions

tm_df <- tm_join_levels %>%
  ggsankey::make_long(
    lg2013, lg2016, lg2019,
    value = count
  )

tm_df2 <-  tm_df %>%
  group_by(x, node) %>%
  summarise(n = sum(value))

tm_df3 <- tm_df %>%
  left_join(tm_df2)

p <- tm_df3 %>%
  ggplot(aes(x = x,
             next_x = next_x,
             node = node,
             next_node = next_node,
             fill = factor(node),
             label = paste0(node,": n = ", n),
             value = value)) +
  geom_sankey(alpha = 0.5) +
  geom_sankey_label(alpha = 0.5, colour = "black") +
  theme_sankey() +
  theme(legend.position = "none")

p

# only changes

tm_df <- tm_join_levels %>%
  filter(stable == "changed") %>%
  ggsankey::make_long(
    lg2013, lg2016, lg2019,
    value = count
  )

tm_df2 <-  tm_df %>%
  group_by(x, node) %>%
  summarise(n = sum(value))

tm_df3 <- tm_df %>%
  left_join(tm_df2)

p <- tm_df3 %>%
  ggplot(aes(x = x,
             next_x = next_x,
             node = node,
             next_node = next_node,
             fill = factor(node),
             label = paste0(node,": n = ", n),
             value = value)) +
  geom_sankey(alpha = 0.5) +
  geom_sankey_label(alpha = 0.5, colour = "black") +
  theme_sankey() +
  theme(legend.position = "none")

p



# map showing stable (TRUE) vs changed (FALSE)
activeCat(temporal_maps) <- "stable"
plot(temporal_maps)
tm_additional_levels %>%
  group_by(stable) %>%
  summarise(n_pixels = sum(count))

tm_binary_levels <- data.frame(
  value = c(0,1,10,11,100,101,110,111)) %>%
  mutate(
    label = stringr::str_pad(value, side = "left", pad = "0", width = 3),
    label2 = c(
      "not present", "gained", "dynamic",
      "gained", "lost", "dynamic", "lost", "stable"),
    col = c("lightgrey", "lightblue", "orange1", "darkblue", "darkred",
            "orange2", "red1",
            "green"),
    col2 = c("lightgrey", "blue", "orange", "blue", "red", "orange", "red",
             "green"))

tm_ts <- terra::catalyze(temporal_maps) %>%
  terra::subset(subset = lg) %>%
  sapp(fun = function(x, cats = tm_binary_levels) {
    nm <- names(x)
    xc <- as.factor(x)
    names(cats)[names(cats) == "label"] <- nm
    levels(xc) <- cats
    names(xc) <- nm
    coltab(xc) <- cats[, c("value", "col")]
    return(xc)
  })
plot(tm_ts)

tm_ts2 <- tm_ts %>%
  sapp(fun = function(x, cats = tm_binary_levels) {
    nm <- names(x)
    activeCat(x) <- "label2"
    names(x) <- nm
    coltab(x) <- cats %>% select(value, col = col2)
    return(x)
  })
plot(tm_ts2)

# apply majority filter, use 3 by 3 block
activeCat(temporal_maps) <- "lg2013_lg2016_lg2019"
temporal_maps_majority <- focal(
  temporal_maps, w = 3, fun = "modal") %>%
  apply_cats(cats = tm_join_levels,
             name = "lg2013_lg2016_lg2019",
             coltab = FALSE)
temporal_maps_majority_twice <- focal(
  temporal_maps, w = 3, fun = "modal") %>%
  focal(w = 3, fun = "modal") %>%
  apply_cats(cats = tm_join_levels,
             name = "lg2013_lg2016_lg2019",
             coltab = FALSE)

# map showing stable (TRUE) vs changed (FALSE)
activeCat(temporal_maps_majority) <- "stable"
plot(temporal_maps_majority)
activeCat(temporal_maps_majority_twice) <- "stable"
plot(temporal_maps_majority_twice)


tm_tsm <- terra::catalyze(temporal_maps_majority) %>%
  terra::subset(subset = lg) %>%
  sapp(fun = function(x, cats = tm_binary_levels) {
    nm <- names(x)
    xc <- as.factor(x)
    names(cats)[names(cats) == "label"] <- nm
    levels(xc) <- cats
    names(xc) <- nm
    coltab(xc) <- cats[, c("value", "col")]
    return(xc)
  })
plot(tm_tsm)
tm_tsm2 <- tm_tsm %>%
  sapp(fun = function(x, cats = tm_binary_levels) {
    nm <- names(x)
    activeCat(x) <- "label2"
    names(x) <- nm
    coltab(x) <- cats %>% select(value, col = col2)
    return(x)
  })
plot(tm_tsm2)

# Vergelijking
##############
activeCat(temporal_maps) <- "stable"
rasterstoplot <- c(
  temporal_stratification,
  temporal_maps,
  temporal_maps_majority,
  temporal_maps_majority_twice)
names(rasterstoplot) <- c(
  "majority filter 9x9 on input maps",
  "no filter",
  "majority filter 3x3 on change maps",
  "majority filter 3x3 on change maps\napplied twice"
  )
plot(rasterstoplot)



# Sample selection
##################

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


