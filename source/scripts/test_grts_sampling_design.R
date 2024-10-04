# setup

library(grtsdb)
library(ggplot2)
library(dplyr)
library(sf)
library(terra)
git_root <- rprojroot::find_root(rprojroot::is_git_root)
flea_data <- gsub(
  pattern = "flea-extent", replacement = "flea-data", x = git_root)

# read data

fleagrts <- rast(file.path(flea_data, "data/c-mon/flea_cmon_level15.tiff"))
lg2013 <- rast(file.path(flea_data, "data", "2013", "LG2013_finaal_update.tif"))
lg2016 <- rast(file.path(flea_data, "data", "2016", "LG2016_finaal_update.tif"))
lg2019 <- rast(file.path(flea_data, "data", "2019", "LG2019_finaal_update.tif"))

# preprocess raster data

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

lg2013 <- resample(lg2013, fleagrts) # needed because of slightly different origin
lg2016 <- resample(lg2016, fleagrts)
lg2019 <- resample(lg2019, fleagrts)

terra::compareGeom(lg2013, lg2016, lg2019)
terra::compareGeom(lg2013, fleagrts)

# combine the spatial strata into temporal strata

datatype(lg2013) # integers 0 to 255 INT1U
# unique(c(), as.raster = TRUE) of numeric layers is safer than concats()

if (file.exists(
  file.path(flea_data, "data/2013_2016_2019", "temporal_stratification.tif"))
  ) {
  temporal_stratification <-
    rast(
      file.path(
        flea_data, "data/2013_2016_2019", "temporal_stratification.tif"
        )
      )
} else {
  temporal_stratification <- unique(
    c(lg2013, lg2016, lg2019) * 1,
    as.raster = TRUE
  )
  # deal with NAs
  temporal_stratification <- mask(temporal_stratification, lg2013)
  temporal_stratification <- mask(temporal_stratification, lg2016)
  temporal_stratification <- mask(temporal_stratification, lg2019)
  temporal_stratification <- droplevels(temporal_stratification)

  binary_change <- function(data, lg) {
    binary <- vector("list", length = length(lg))
    binary <- setNames(binary, lg)
    for (i in lg) {
      binary[[i]] <- paste0(
        stringr::str_detect(data$lg2013_label, i) %>% as.numeric(),
        stringr::str_detect(data$lg2016_label, i) %>% as.numeric(),
        stringr::str_detect(data$lg2019_label, i) %>% as.numeric())
    }
    bind_cols(data, binary)
  }

  categorize_land_use_change <- function(B, simple = TRUE) {
    if (simple) {
      case_when(
        # Stable conditions
        grepl("^0+$", B) ~ "Stable absence",
        grepl("^1+$", B) ~ "Stable presence",

        # Simple changes
        grepl("^0+1+$", B) ~ "Gain",
        grepl("^1+0+$", B) ~ "Loss",

        # Default case
        TRUE ~ "Other complex pattern")
    } else {
      case_when(
        # Stable conditions
        grepl("^0+$", B) ~ "Stable absence",
        grepl("^1+$", B) ~ "Stable presence",

        # Simple changes
        grepl("^0+1+$", B) ~ "Gain",
        grepl("^1+0+$", B) ~ "Loss",

        # Complex changes
        grepl("^0+1+0+$", B) ~ "Temporary gain",
        grepl("^1+0+1+$", B) ~ "Temporary loss",
        grepl("^0+1+0+1+$", B) ~ "Intermittent presence (starting absent)",
        grepl("^1+0+1+0+$", B) ~ "Intermittent presence (starting present)",

        # Oscillating changes
        grepl("^(01)+0?$", B) ~ "Oscillating (starting absent)",
        grepl("^(10)+1?$", B) ~ "Oscillating (starting present)",

        # Other complex patterns
        grepl("01.*1$", B) & !grepl("^0+1+$", B) ~ "Complex gain",
        grepl("10.*0$", B) & !grepl("^1+0+$", B) ~ "Complex loss",

        # Default case
        TRUE ~ "Other complex pattern"
      )
    }
  }

  lg <- gsub(pattern = "^\\d\\s-\\s", replacement = "", x = catstable$label)

  additional_levels <- freq(temporal_stratification) %>%
    as_tibble() %>%
    tidyr::separate(
      value, into = c("lg2013", "lg2016", "lg2019"),
      sep = "_", remove = FALSE) %>%
    left_join(
      catstable %>%
        mutate(
          value = as.character(value),
          lg2013_label = label,
          .keep = "none"
        ) ,
      by = join_by(lg2013 == value)) %>%
    left_join(
      catstable %>%
        mutate(
          value = as.character(value),
          lg2016_label = label,
          .keep = "none"
        ) ,
      by = join_by(lg2016 == value)) %>%
    left_join(
      catstable %>%
        mutate(
          value = as.character(value),
          lg2019_label = label,
          .keep = "none"
        ) ,
      by = join_by(lg2019 == value)) %>%
    binary_change(lg = lg) %>%
    rowwise() %>%
    mutate(stable = ifelse(
      all(lg2013 == lg2016, lg2016 == lg2019),
      "stable", "changed") %>%
        as.factor()) %>%
    ungroup() %>%
    mutate(
      across(
        all_of(lg),
        \(x) categorize_land_use_change(x),
        .names = "{.col}_changecat"
      )
    )

  join_levels <- cats(temporal_stratification)[[1]] %>%
    mutate(across(starts_with("lg"), as.character)) %>%
    inner_join(
      additional_levels,
      by = join_by(lg2013, lg2016, lg2019, label == value))
  levels(temporal_stratification) <- join_levels
  coltab(temporal_stratification) <- NULL

  writeRaster(
    x = temporal_stratification,
    filename =
      file.path(flea_data, "data/2013_2016_2019", "temporal_stratification.tif"),
    overwrite = FALSE
  )
}

plot(temporal_stratification)
plot(`activeCat<-`(temporal_stratification, "stable"))
plot(`activeCat<-`(temporal_stratification, "Urbaan"))
plot(`activeCat<-`(temporal_stratification, "Urbaan_changecat"))


# calculate for each pixel the dominant temporal stratum inside a 9x9 block centered
# on the focal pixel
temporal_stratification_modal9 <- focal(
  temporal_stratification, w = 9, fun = "modal")


# Source C++ function to count unique values in window
library(Rcpp)
sourceCpp(here::here("source/scripts/unique-landuse-count.cpp"))
# Then use it with focalCpp
temporal_stratification_countunique9 <- terra::focalCpp(
  x = temporal_stratification, w = 9, count_unique_landuse)
hist(temporal_stratification_countunique9, maxcell = 1e7)
plot(temporal_stratification_countunique9)

jointable <- cats(temporal_stratification)[[1]] |>
  as_tibble() |>
  select(-count, -layer)

names(temporal_stratification_modal9) <- "value"
temporal_stratification_modal9 <- as.factor(temporal_stratification_modal9)

m9catstable <- freq(temporal_stratification_modal9) |>
  as_tibble() |>
  inner_join(
    jointable %>% mutate(ID = as.character(ID)),
    by = join_by(value == ID)
  ) |>
  select(-layer) |>
  mutate(value = as.numeric(value))

levels(temporal_stratification_modal9) <- m9catstable
temporal_stratification_modal9
plot(`activeCat<-`(temporal_stratification_modal9, "label"))
plot(`activeCat<-`(temporal_stratification_modal9, "Urbaan_changecat"))

# sample extraction
# for each land use category

# 1. Set the active category to the desired factor
lu <- "Open natuur_changecat"
ts1 <- `activeCat<-`(temporal_stratification, lu)
# 3. Adjust the number of levels
levelvec <- levels(ts1)[[1]][[2]]
levelvec[levelvec == "Stable absence"] <- NA
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
plot(ts2)

ts2 <- ts2 %>% segregate(classes = NULL, other = NA)
names(ts2) <- levels(as.factor(levelvec))
plot(ts2)
freq(ts2)

fleagrts_ts2 <- crop(fleagrts, ts2, mask = TRUE)
names(fleagrts_ts2) <- levels(as.factor(levelvec))
plot(fleagrts_ts2)


extract_sample <- function(rast, n) {
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

  return(selected)
}
sample20 <- lapply(fleagrts_ts2, extract_sample, n = 20)
gain <- st_as_sf(sample20[[1]], coords = c("x", "y"), crs = 31370)
plot(fleagrts_ts2[[1]])
points(vect(gain))

one_point <- gain %>%
  slice(4) %>%
  terra::vect()
one_cell <- one_point %>%
  terra::buffer(5) %>%
  terra::crop(x = lg2019, y = .)
plot(one_cell)
points(one_point)
polys(one_point %>% buffer(5))
