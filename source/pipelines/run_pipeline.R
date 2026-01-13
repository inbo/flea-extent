################
# run pipeline #
################
library(targets)

Sys.setenv(TAR_PROJECT = "validation_sample")

# check status
tar_visnetwork(targets_only = TRUE, physics = TRUE)
#tar_manifest() |> View()

# check reason if any is outdated
tar_sitrep() |> dplyr::filter(dplyr::if_any(.cols = !c(name)))

# run the pipeline
px <- tar_make(as_job = TRUE, use_crew = TRUE)

# # run with profiling
# results <- profvis::profvis(
#   targets::tar_make(
#     callr_function = NULL, # Do not run the pipeline behind a callr::r() process.
#     use_crew = FALSE, # Disable parallel computing with crew (optional)
#     as_job = FALSE # Do not run the pipeline in a Posit Workbench / RStudio background job.
#   )
# )
# print(results, aggregate = TRUE) # aggregate = TRUE is crucial for interpretable flame graphs.


tar_progress_summary()
tar_poll()

# stop process
# ps::ps_kill(px$as_ps_handle())

####################
# inspect pipeline #
####################

#targets::tar_prune()

mt <- targets::tar_meta(
  fields = error,
  complete_only = TRUE
)
mt
#View(mt)
wn <- targets::tar_meta(fields = warnings, complete_only = TRUE)
wn
#View(wn)
targets::tar_visnetwork(
  label = c("description", "time", "size"),
  targets_only = TRUE, physics = TRUE
)


# logging
library(autometric)
log_file <- "log.txt"
log_data <- log_read(log_file)
log_plot(log_data, metric = "resident")

#tar_load_globals()
#tar_load(names = c(mapnames, catstable, grts_ext, grts_origin))

tar_read(mapnames)
ct <- tar_read(catstable)
ml <- tar_read(maps)


ml[[1]]
my_cats <- terra::cats(ml[[1]])
my_cats
coltab <- my_cats[[1]] |> dplyr::select(value, red, green, blue)
terra::coltab(ml[[1]]) <- coltab
terra::plot(ml[[1]], colNA = "white")
terra::values(ml[[1]], row = 5000, nrows = 1)
terra::datatype(ml[[1]])
terra::NAflag(ml[[1]])
ft <- terra::freq(ml[[3]])
ft |>
  dplyr::mutate(prop = round(count / sum(count), 4))





grts <- tar_read(fleagrts)
grts
terra::datatype(grts)
terra::compareGeom(grts, ml[[1]])
waldo::compare(
  terra::crs(grts),
  terra::crs(ml[[1]])
  )

tm <- targets::tar_read(temporal_map)
tm
terra::plot(tm)

terra::cats(tm)[[1]] |> head()
terra::cats(tm)[[1]] |> tail()

tms <- targets::tar_read(temporal_map_strata)
tms
terra::plot(tms)
terra::values(tms, row = 1, nrows = 1)

tms_cats <- terra::cats(tms)
terra::activeCat(tms) <- "stable"
terra::plot(tms, colNA = "white")

targets::tar_read(lu_changecats)

sg <- targets::tar_read(separate_grts)
all(purrr::map(sg, ~inherits(.x, "SpatRaster")) |> unlist())
sg[[1]]


vs <- targets::tar_read(validation_sample)
terra::vect(vs) |>
  sf::st_as_sf(crs = 31370) |>
  sf::st_drop_geometry() |>
  dplyr::count(grts_rank, name = "times_selected") |>
  dplyr::count(times_selected) |>
  dplyr::mutate(prop = n/ sum(n))

terra::vect(vs) |>
  sf::st_as_sf(crs = 31370) |>
  sf::st_drop_geometry() |>
  dplyr::count(stratum_name, changecat) |>
  tidyr::pivot_wider(names_from = changecat, values_from = n) |>
  View()

vp <- targets::tar_read(validation_polygons)
lapply(vp, nrow) |> unlist() |> sum()
# =3248 + 824*2 + 48*3 of 3248 + 824 + 48 = 4120 unieke rastercellen

library(ggplot2)
library(sf)
terra::vect(vs) |>
  st_as_sf() |>
  dplyr::mutate(
    value = stringr::str_extract(
      stratum_name, "^bc_(.*)_changecat", group = 1
    ) |> as.numeric()
  ) |>
  dplyr::inner_join(ct, by = "value") |>
  ggplot() +
  geom_sf(aes(colour  = changecat), alpha = 0.2) +
  facet_wrap(~ stringr::str_wrap(label, 15))

terra::vect(vs) |>
  st_as_sf() |>
  st_drop_geometry() |>
  dplyr::count(stratum_name, changecat)

grb_parc <- tar_read(grb_parcels)
grb_parc <- terra::vect(grb_parc) |>
  st_as_sf() |>
  dplyr::mutate(source = "parcels")
grb_set <- tar_read(grb_settlements)
grb_set <- terra::vect(grb_set) |>
  st_as_sf() |>
  dplyr::mutate(source = "settlements")
grb_water <- tar_read(grb_waterways)
grb_water <- terra::vect(grb_water) |>
  st_as_sf() |>
  dplyr::mutate(source = "water")

grb_water_set <- dplyr::bind_rows(
  grb_water,
  grb_set)

grb_water_set |>
  mapview::mapview(zcol = "source", alpha.regions = 0.2) +
  mapview::mapview(grb_parc, alpha.region = 0)

lbg_101 <- tar_read(lbg_101_cropped)
lbg_104 <- tar_read(lbg_104_cropped)

mapview::mapview(terra::vect(lbg_101), alpha.regions = 0.2
                 , col.regions = "orange") +
  mapview::mapview(terra::vect(lbg_104), alpha.regions = 0.2,
                   col.regions = "yellow") +
  mapview::mapview(terra::vect(vp), alpha.regions = 0)

prelabeled <- tar_read(prelabeled_validation_polygons)

# 102
mapview::mapview(prelabeled$prelabeled_validation_polygons_b0bae0239048f770,
                 zcol = "labels")

# 106
mapview::mapview(prelabeled$prelabeled_validation_polygons_c66dc9ac1a4d5ea7,
                 zcol = "labels")

# 200
mapview::mapview(prelabeled$prelabeled_validation_polygons_0692faf9b384816b,
                 zcol = "labels")



##################
# debug pipeline #
##################

targets::tar_load_globals()
debugonce(process_settlement)
targets::tar_load(grb_settlements)
process_settlement(grb = grb_settlements)


debugonce(get_grb_by_row)
test <- get_grb_by_row(
  layer = "GRB:ADP",
  polygons = tar_read(validation_polygons_6e7d3123e950eb4d)[1:2,]
)
targets::tar_load_globals()
targets::tar_workspace("vp_water_settlements_5319be99c3d05901")
debugonce(combine_water_settlements)
test <- combine_water_settlements(
  water = vp_water,
  settlements = grb_settlements_processed,
  polygons = validation_polygons,
  lbg_101 = lbg_101_cropped,
  lbg_104 = lbg_104_cropped
)

targets::tar_load_globals()
targets::tar_load(
  names = c(
  grb_settlements_processed,
  vp_water,
  validation_polygons,
  lbg_101_cropped,
  lbg_104_cropped))

debugonce(combine_water_settlements)
test <-  combine_water_settlements(
  water = vp_water$vp_water_f39e59863b95dc86,
  settlements = grb_settlements_processed,
  lbg_101 = lbg_101_cropped,
  lbg_104 = lbg_104_cropped,
  polygons = validation_polygons
)


targets::tar_load_globals()
targets::tar_workspace("watersurfaces_processed_d90ebef693247dc3")
debugonce(get_watersurfaces)
test <- get_watersurfaces(
  path_version = zenodo_watersurface,
  polygons = validation_polygons,
  meta = watersurfaces_meta
)


targets::tar_workspace("prelabeled_validation_polygons_50_04ec7c30fc6c9565")
debugonce(crop_labeled_polygons)
debug(combine_small_intersections)
test <- crop_labeled_polygons(
  pvp = prelabeled_validation_polygons,
  crop_with = validation_polygons_50
)

targets::tar_workspace("vp_water_grb_lbg_5319be99c3d05901")

debugonce(combine_water_grb_lbg)

result <- combine_water_grb_lbg(
  water = vp_water,
  settlements = grb_settlements_processed,
  lbg_101 = lbg_101_cropped,
  lbg_104 = lbg_104_cropped,
  lbg_200 = lbg_200_cropped,
  lbg_300 = lbg_300_cropped,
  polygons = validation_polygons
)

targets::tar_workspace(prelabeled_validation_polygons_6fdb9209e87c593c)
debugonce(intersect_validation_polygons)
debug(combine_small_intersections)
result <- intersect_validation_polygons(
  wsp_target = vp_water_grb_lbg_singletarget,
  lu_changecat = lu_changecats,
  input_years = input_years,
  area_too_small = 10
)

hist(expanse(result))
sum(expanse(result) < 10)

targets::tar_workspace(vp_water_grb_lbg_5319be99c3d05901)
debugonce(combine_water_grb_lbg)
result <- combine_water_grb_lbg(
  water = vp_water,
  settlements = grb_settlements_processed,
  lbg_list = list(!!!lbg_symbols),
  polygons = validation_polygons
)

tar_manifest(fields = "command", names = vp_water_grb_lbg) |> View()

targets::tar_workspace(raster_labeled_validation_polygons_50_de2d47ebdd431bdb)
debugonce(raster_to_polygons)
result <- raster_to_polygons(
  tm = temporal_map_strata, # temporal raster
  vp = validation_polygons_50, # validation polygons
  input_years = as.character(input_years)
)
