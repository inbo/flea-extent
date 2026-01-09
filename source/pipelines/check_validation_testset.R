library(googlesheets4)
library(dplyr)
library(ggplot2)
library(sf)
library(here)
flea_data <- gsub("extent", "data", here())


notes <- read_sheet(ss = "1KkhhJgw9Gkru3XgA17d7YNz9fCG50f0c6yLysLVag6s")
glimpse(notes)
summary(notes)

validation_polys_po <- file.path(
  flea_data,
  "validation",
  "test_validation_sample_30_originalworkflow_patrik_20260106.gpkg"
) |>
  read_sf(
    layer = "test_validation_sample_30_originalworkflow_patrik_20260106"
  )
glimpse(validation_polys_po)
validation_polys_po <- validation_polys_po |>
  rename(
    prelabels = labels
  ) |>
  mutate(
    newlabels = paste(label_2016, label_2019, label_2022, sep = "-"),
    validator = "PO",
    method = "polygons"
  )

# prop polygonen met juist prelabel
sum(validation_polys_po$prelabels == validation_polys_po$newlabels) / nrow(validation_polys_po)

# check polygonen onvolledig gevalideerd
onvolledig <- grepl("other", validation_polys_po$newlabels)
validation_polys_po[
  onvolledig, c("grts_rank", "label_2016", "label_2019", "label_2022")
] |>
  st_drop_geometry()

# check oppervlakte
validation_polys_po %>%
  mutate(
    area_m2 = st_area(.)
  ) |>
  st_drop_geometry() |>
  group_by(grts_rank) |>
  summarise(area_m2 = sum(area_m2)) |>
  filter(
    !between(as.numeric(area_m2), 2490, 2501)
  )

# validation cells PO
validation_cells_po <- file.path(
  flea_data,
  "validation",
  "test_validation_30_rasterizedsamples_patrik_20260108.gpkg"
) |>
  read_sf(
    layer = "test_validation_30_rasterizedsamples_patrik_20260108"
  )

validation_cells_po <- validation_cells_po |>
  rename(
    prelabels = label_label
  ) |>
  mutate(
    newlabels = paste(label_2016, label_2019, label_2022, sep = "-"),
    validator = "PO",
    method = "cells"
  )

# prop polygonen met juist prelabel
sum(is.na(validation_cells_po$prelabels))

sum(
  validation_cells_po$prelabels == validation_cells_po$newlabels, na.rm = TRUE
) /
  nrow(validation_cells_po)

# check polygonen onvolledig gevalideerd
onvolledig <- grepl("other", validation_cells_po$newlabels)
validation_cells_po[
  onvolledig, c("grts_rank", "label_2016", "label_2019", "label_2022")
] |>
  st_drop_geometry()

# check oppervlakte
validation_cells_po %>%
  mutate(
    area_m2 = st_area(.)
  ) |>
  st_drop_geometry() |>
  group_by(grts_rank) |>
  summarise(area_m2 = sum(area_m2)) |>
  dplyr::filter(
    !between(as.numeric(area_m2), 2490, 2501)
  )



# validation polys Steven
validation_polys_sds <- file.path(
  flea_data,
  "validation",
  "steven_ds",
  "FLEA_testset_30_validatie_202512_sds.shp"
) |>
  read_sf()


glimpse(validation_polys_sds)
validation_polys_sds <- validation_polys_sds |>
  rename(
    prelabels = labels
  ) |>
  mutate(
    newlabels = paste(label_2016, label_2019, label_2022, sep = "-"),
    validator = "SDS",
    method = "polygons"
  )

# prop polygonen met juist prelabel
sum(validation_polys_sds$prelabels == validation_polys_sds$newlabels) /
  nrow(validation_polys_sds)

# check polygonen onvolledig gevalideerd
onvolledig <- grepl("other|water", validation_polys_sds$newlabels)
validation_polys_sds[
  onvolledig, c("grts_rank", "label_2016", "label_2019", "label_2022")
] |>
  st_drop_geometry()

# check oppervlakte
validation_polys_sds %>%
  mutate(
    area_m2 = st_area(.)
  ) |>
  st_drop_geometry() |>
  group_by(grts_rank) |>
  summarise(area_m2 = sum(area_m2)) |>
  filter(
    !between(as.numeric(area_m2), 2490, 2501)
  )


# validation cells Steven
validation_cells_sds <- file.path(
  flea_data,
  "validation",
  "steven_ds",
  "FLEA_test_cells_30_validatie_202512_sds.shp"
) |>
  read_sf()

validation_cells_sds <- validation_cells_sds |>
  rename(
    prelabels = label_labe
  ) |>
  mutate(
    newlabels = paste(label_2016, label_2019, label_2022, sep = "-"),
    validator = "SDS",
    method = "cells"
  )
# prop polygonen met juist prelabel
sum(validation_cells_sds$prelabels == validation_cells_sds$newlabels) /
  nrow(validation_cells_sds)

# check polygonen onvolledig gevalideerd
onvolledig <- grepl("other", validation_cells_sds$newlabels)
validation_cells_sds[
  onvolledig, c("grts_rank", "label_2016", "label_2019", "label_2022")
] |>
  st_drop_geometry()

# check oppervlakte
validation_cells_sds %>%
  mutate(
    area_m2 = st_area(.)
  ) |>
  st_drop_geometry() |>
  group_by(grts_rank) |>
  summarise(area_m2 = sum(area_m2)) |>
  filter(
    !between(as.numeric(area_m2), 2490, 2501)
  )

# check joins
notes |>
  anti_join(
    bind_rows(validation_polys_po, validation_polys_sds),
    by = join_by(validator, method, grts_rank)
  )

validation_cells_po <- validation_cells_po |>
  dplyr::filter(!st_is_empty(validation_cells_po))

nrow(validation_cells_po) == 30 * 25
nrow(validation_cells_sds) == 30 * 25

validation_cells <- bind_rows(
  validation_cells_po,
  validation_cells_sds |>
    select(-fid_) |>
    rename(
      geom = geometry,
      stratum_name = stratum_na,
      area_labeled = area_label,
      label_2016_2 = label_2017,
      label_2019_2 = label_2020,
      label_2022_2 = label_2023
    )
)


validation_polys_sds <- validation_polys_sds |>
  select(-fid_)
validation_polys_sds_old <- validation_polys_sds
names(validation_polys_sds) <- names(validation_polys_po)
validation_polys_sds <- st_set_geometry(validation_polys_sds, "geom")

validation_polys <- bind_rows(
  validation_polys_po,
  validation_polys_sds
)

notes <- notes |>
  left_join(validation_polys_po |>
              distinct(grts_rank, stratum_name, changecat, area_labeled))

############################################################################

notes |>
  ggplot() +
  geom_point(
    aes(x = method, y = time_minutes, group = grts_rank)
  ) +
  geom_line(
    aes(x = method, y = time_minutes, group = grts_rank)
  ) +
  stat_summary(
    fun.data = mean_cl_boot,
    aes(x = method, y = time_minutes),
    colour = INBOtheme::inbo_felrood
  ) +
  facet_wrap(~validator)


notes |>
  ggplot() +
  geom_point(
    aes(x = method, y = time_minutes, group = grts_rank, colour = stratum_name),
    position = position_dodge(width = 0.2), alpha = 0.3
  ) +
  geom_line(
    aes(x = method, y = time_minutes, group = grts_rank, colour = stratum_name),
    position = position_dodge(width = 0.2), alpha = 0.3
  ) +
  stat_summary(
    fun.data = mean_cl_boot,
    aes(x = method, y = time_minutes, colour = stratum_name),
    position = position_dodge(width = 0.2)
  ) +
  facet_wrap(~validator)

#######################################
library(tidyr)
compare_cells <- validation_cells |>
  st_drop_geometry() |>
  select(
    validator, method, grts_rank, label_2016, label_2019, label_2022
  ) |>
  bind_cols(
    as_tibble(st_coordinates(st_centroid(st_geometry(validation_cells))))
  ) |>
  pivot_longer(
    cols = starts_with('label_'),
    names_to = "year",
    names_prefix = "label_",
    values_to = "label"
  ) |>
  pivot_wider(
    id_cols = c(grts_rank, X, Y, method, year),
    names_from = validator,
    values_from = label
  ) |>
  count(PO, SDS) |>
  mutate(
    PO = readr::parse_number(PO),
    SDS = readr::parse_number(SDS),
    beoordeling = ifelse(PO == SDS, "zelfde", "verschillend"),
    percentage =  n / sum(n) * 100
  ) |>
  arrange(beoordeling, -n)
compare_cells |>
  knitr::kable(digits = 1)


compare_cells |>
  group_by(beoordeling) |>
  summarise(
    n = sum(n),
    percentage = sum(percentage)
  ) |>
  dplyr::filter(!is.na(beoordeling)) |>
  knitr::kable(digits = 0)

# inschatting tijd nodig voor volledige validatie
library(targets)
Sys.setenv(TAR_PROJECT = "validation_sample")

tar_load(validation_polygons_50)
prelabeled <- tar_read(prelabeled_validation_polygons_50)
rasterlabeled <- tar_read(raster_labeled_validation_polygons_50)

library(terra)
validation_polygons_50 <- vect(validation_polygons_50)

time_estimates <- values(validation_polygons_50) |>
  as_tibble() |>
  count(stratum_name, changecat) |>
  mutate(
    stratum_name = stringr::str_extract(stratum_name, "\\d+"),
    level = ifelse(stratum_name %in% c(101:106), 2, 1),
    days = ifelse(level == 1, 5 * n, 10 * n) / 60 / 7.36
  )

time_estimates |>
  dplyr::filter(
    !stratum_name %in% c("600", "1000")
  ) |>
  summarise(
    n_halfway = sum(n) / 2,
    n_all = sum(n),
    days_halfway = sum(days) / 2,
    days_all = sum(days),
    pm_halfway = days_halfway / (160/12),
    pm_all = days_all / (160/12)
  )


######################################
validation_polys_po
?terra::rasterize()

one_po <- validation_polys_po |>
  dplyr::filter(grts_rank == 8820)
one_sds <- validation_polys_sds |>
  dplyr::filter(grts_rank == 8820)

plot(st_geometry(one_sds))
plot(st_geometry(one_po))

one_sds |>
  st_intersection(one_po) %>%
  mutate(
    area = st_area(.),
    id = seq_len(n())
  ) |>
  st_drop_geometry() |>
  select(
    id,
    area,
    label_2016_sds = label_2016,
    label_2019_sds = label_2019,
    label_2022_sds = label_2022,
    label_2016_po = label_2016.1,
    label_2019_po = label_2019.1,
    label_2022_po = label_2022.1
  ) |>
  pivot_longer(
    cols = starts_with('label_'),
    names_to = c("year", "validator"),
    names_prefix = "label_",
    names_sep = "_",
    values_to = "label"
  ) |>
  pivot_wider(
    id_cols = c(id, area, year),
    names_from = validator,
    values_from = label
  ) |>
  mutate(
    beoordeling = ifelse(sds == po, "zelfde", "verschillend")
  ) |>
  group_by(
    year, beoordeling
  ) |>
  summarise(
    area = sum(area)
  )

compare_polys <- st_make_valid(validation_polys_sds) |>
  st_intersection(
    st_make_valid(validation_polys_po)
  ) %>%
  mutate(
    area = st_area(.),
    id = seq_len(n())
  ) |>
  st_drop_geometry() |>
  select(
    id,
    grts_rank,
    area,
    label_2016_sds = label_2016,
    label_2019_sds = label_2019,
    label_2022_sds = label_2022,
    label_2016_po = label_2016.1,
    label_2019_po = label_2019.1,
    label_2022_po = label_2022.1
  ) |>
  pivot_longer(
    cols = starts_with('label_'),
    names_to = c("year", "validator"),
    names_prefix = "label_",
    names_sep = "_",
    values_to = "label"
  ) |>
  pivot_wider(
    id_cols = c(id, grts_rank, area, year),
    names_from = validator,
    values_from = label
  ) |>
  mutate(
    beoordeling = ifelse(sds == po, "zelfde", "verschillend")
  ) |>
  group_by(
    grts_rank, year, beoordeling
  ) |>
  summarise(
    area = sum(area)
  )


library(units)
compare_polys |>
  ggplot() +
  geom_col(
    aes(x = year, y = area, fill = beoordeling)
  ) +
  facet_wrap(~ grts_rank)


