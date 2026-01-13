library(spbal)
library(ggplot2)
library(dplyr)
library(sf)
library(terra)




git_root <- rprojroot::find_root(rprojroot::is_git_root)

flea_data <- gsub(
  pattern = "flea-extent", replacement = "flea-data", x = git_root
)

lg2013 <- rast(file.path(flea_data, "data", "2013", "LG2013_finaal_update.tif"))
plot(lg2013)

bbox_lg2013 <- terra::ext(lg2013)
temp <- lg2013
non_na_cells <- which(!is.na(values(temp)))
values(temp)[non_na_cells] <- 1
plot(temp)
sf_lg2013 <- as.polygons(temp) %>% st_as_sf()

# the following creates a very large object (> 7 Gb) and takes > 10 minutes
# sf_lg2013 <- as.polygons(lg2013) %>% st_as_sf() # nolint

# Vertically aligned master sample bounding box.
bb <- spbal::BoundingBox(shapefile = sf_lg2013)
plot(bb)

set.seed(511)
n_samples <- 20
result <- spbal::BAS(
  shapefile = sf_lg2013,
  n = n_samples,
  boundingbox = bb
)
bas20 <- result$sample

plot(lg2013)
points(vect(bas20))

# pass the result$seed to increase sample size
n_samples <- 50
result2 <- spbal::BAS(
  shapefile = sf_lg2013,
  n = n_samples,
  boundingbox = bb,
  seeds = result$seed
)
bas50 <- result2$sample
bas50



### maybe better to convert raster to points sf object and then use
# ?spbal::HIP


# # try sgsR::sample_balanced - takes too much time?
# test <- sgsR::sample_balanced(
#   temp,
#   200
# )



testraster <- terra::rast(
  nrows = 500, ncols = 1000, crs = "EPSG:31370"
)
values(testraster) <- sample(
  c(NA, 1), size = 500000, prob = c(0.4, 0.6), replace = TRUE
)
plot(testraster)


# Get coordinates of valid cells
microbenchmark::microbenchmark(
  {
    valid_cells <- which(!is.na(values(testraster)))
    xyFromCell(testraster, valid_cells)
  },
  {
    valid_cells <- which(!is.na(values(testraster)))
    rowColFromCell(testraster, valid_cells)
  },
  {
    crds(testraster, na.rm = TRUE)
  }
)

spatially_balanced_sample <- function(rast, n, method = "HIP", seed = NULL) {
  require("terra")
  require("sf")
  require("spbal")

  if (!is.null(seed)) set.seed(seed)

  # Perform sampling based on the chosen method
  if (method == "HIP") {
    # Get non-NA cell indices
    valid_cells <- which(!is.na(values(rast)))
    # Get coordinates of valid cells
    coords <- xyFromCell(rast, valid_cells)
    # Create an sf object with valid cell centroids
    points_sf <- st_as_sf(
      as.data.frame(coords), coords = c("x", "y"),
      crs = crs(rast)
    )
    sample_result <- HIP(points_sf, n)
    sample_result$Population <- NULL
  } else if (method == "BAS") {
    polys_sf <- as.polygons(rast) %>% st_as_sf()
    bb <- spbal::BoundingBox(shapefile = polys_sf)
    sample_result <- BAS(shapefile = polys_sf, n = n, boundingbox = bb)
  } else {
    stop("Unsupported sampling method. Please use 'HIP' or 'BAS'.")
  }

  return(sample_result)
}
microbenchmark::microbenchmark(
  hip <- spatially_balanced_sample(testraster, n = 100, method = "HIP"),
  bas <- spatially_balanced_sample(testraster, n = 100, method = "BAS"),
  times = 5L
)

bas_balanced_sample <- function(rast, n, stratum, seed = NULL) {
  require("terra")
  require("sf")
  require("spbal")

  if (!is.null(seed)) set.seed(seed)

  polys_sf <- as.polygons(rast, round = FALSE, aggregate = FALSE) %>%
    st_as_sf()
  bb <- spbal::BoundingBox(shapefile = polys_sf)
  sample_result <- BAS(
    shapefile = polys_sf,
    n = n,
    stratum = stratum,
    boundingbox = bb
  )

  return(sample_result)
}

lg2013_selectie <- crop(lg2013, ext(200000, 205000, 200000, 205000))
as.polygons(lg2013_selectie) %>% st_as_sf()
nvec <- rep(10, 9)
names(nvec) <- as.character(1:9)
test <- bas_balanced_sample(
  rast = lg2013_selectie,
  n = nvec,
  stratum = "LG2013_finaal_update"
)
test$sample |>
  ggplot() +
  geom_sf(aes(colour = LG2013_finaal_update))

# lazy conversion to points
lg2013_strat_points <- as.points(lg2013_selectie)

# conversion to sf brings into memory (but see lazysf)
lg2013_strat_points_sf <- st_as_sf(lg2013_strat_points)
lg2013_strat_points_df <- lg2013_strat_points_sf %>%
  bind_cols(st_coordinates(.)) %>%
  st_drop_geometry()

library(SamplingBigData)
sample_size <- 100
minimum_n_strat <- 5
lg2013_sample_strat <- vector(
  "list",
  length = length(unique(lg2013_strat_points_df$LG2013_finaal_update))
)
lg2013_sample_strat <- setNames(
  lg2013_sample_strat,
  nm = unique(lg2013_strat_points_df$LG2013_finaal_update)
)
set.seed(214)
for (i in unique(lg2013_strat_points_df$LG2013_finaal_update)) {
  size_pop <- nrow(lg2013_strat_points_df)
  df <- lg2013_strat_points_df %>%
    filter(LG2013_finaal_update == i)
  size_strat <- nrow(df)
  n_strat <- max(minimum_n_strat, round(size_strat / size_pop * sample_size))
  ips <- rep(n_strat / size_strat, size_strat)
  rowindex <- lpm2_kdtree(
    prob = ips,
    x = df[, c("X", "Y")],
    inOrder = TRUE
  )
  lg2013_sample_strat[[i]] <- df %>%
    slice(rowindex) %>%
    mutate(order = seq_len(n()))
}

lg2013_sample_strat <- bind_rows(lg2013_sample_strat)

lg2013_sample_strat_sf <- st_as_sf(
  lg2013_sample_strat,
  coords = c("X", "Y"),
  crs = st_crs(lg2013_strat_points_sf)
)



lg2013_sample_strat_vect <- vect(lg2013_sample_strat_sf)

# dim = 500 * 500 < 2^9 * 3^6 # nolint
haltonframe <- spbal::HaltonFrame(
  N = 1,
  J = base::c(9, 6),
  bases = base::c(2, 3),
  boundingbox = NULL,
  shapefile = as.polygons(lg2013_selectie) %>% st_as_sf(),
  panels = NULL,
  panel_overlap = NULL,
  seeds = NULL,
  stratum = NULL,
  verbose = TRUE
)

haltonframe$J
haltonframe$bb
haltonframe$seeds
haltonframe$hf.pts.shp %>%
  filter(spbalSeqID < 100) %>%
  ggplot() +
  geom_sf()
