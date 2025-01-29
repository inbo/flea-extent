################
# run pipeline #
################

library(targets)

Sys.setenv(TAR_PROJECT = "validation_sample")

# check status
tar_visnetwork(level_separation = 5000)

# check reason if any is outdated
tar_sitrep() |> dplyr::filter(dplyr::if_any(.cols = !c(name)))

# run the pipeline
tar_make()

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
targets::tar_meta(fields = warnings, complete_only = TRUE)
targets::tar_visnetwork(label = c("description", "time", "size"))


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
terra::cats(ml[[1]])
settlement_mask <- app(
  ml[[1]],
  fun = function(x) {
    x[!x %in% c(101, 102, 105, 106)] <- NA
    return(x)
  }
)
plot(settlement_mask, colNA = "snow4")

terra::plot(ml[[1]], colNA = "orange")
terra::values(ml[[1]], row = 5000, nrows = 1)
terra::coltab(ml[[1]])
terra::cats(ml[[1]])
terra::datatype(ml[[1]])
terra::NAflag(ml[[1]]) # not preserved!
ft <- terra::freq(ml[[3]])
ft |>
  mutate(prop = round(count / sum(count), 4))

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


terra::activeCat(tms) <- "stable"
terra::plot(tms, colNA = "orange")

targets::tar_read(lu_changecats)

sg <- targets::tar_read(separate_grts)
all(purrr::map(sg, ~inherits(.x, "SpatRaster")) |> unlist())

vs <- targets::tar_read(validation_sample)
terra::vect(vs) |> sf::st_as_sf(crs = 31370) |>
  sf::st_drop_geometry() |>
  dplyr::count(grts_rank) |> dplyr::count(n)
terra::vect(vs) |> sf::st_as_sf(crs = 31370) |>
  sf::st_drop_geometry() |>
  dplyr::count(stratum_name, changecat) |>
  tidyr::pivot_wider(names_from = changecat, values_from = n)

vp <- targets::tar_read(validation_polygons)
lapply(vp, nrow) |> unlist() |> sum()


library(ggplot2)
library(sf)
terra::vect(vs) |>
  st_as_sf() |>
  ggplot() +
  geom_sf(aes(colour  = changecat), alpha = 0.2) +
  facet_wrap(~ stratum_name)

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

dplyr::bind_rows(
  grb_water,
  grb_set) |>
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

mapview::mapview(prelabeled$prelabeled_validation_polygons_b0bae0239048f770)

prelabeled$prelabeled_validation_polygons_b0bae0239048f770 |>
  sfst_as_sf()


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
  polygons = validation_polygons
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
targets::tar_workspace("vp_water_settlements_cleaned_57abadb85e4ad4bc")
debugonce(postprocess_water_settlements)
test <- postprocess_water_settlements(
  vp_water_settlements
)


targets::tar_load_globals()
targets::tar_workspace("prelabeled_validation_polygons_e9e42187a231237d")
debugonce(intersect_validation_polygons)
test <- intersect_validation_polygons(
 vp_water_settlements_singletarget, lu_changecats
)
