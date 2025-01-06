library(targets)

Sys.setenv(TAR_PROJECT = "validation_sample")

tar_make()

# debugging and inspection

targets::tar_meta(
  fields = error,
  complete_only = TRUE
)

targets::tar_visnetwork()

log <- autometric::log_read("log.txt")
library(ggplot2)
log |>
  ggplot() +
  geom_line(aes(x = time, y = cpu, colour = factor(pid)))
log |>
  ggplot() +
  geom_line(aes(x = time, y = resident, colour = factor(pid)))


#tar_load_globals()
#tar_load(names = c(mapnames, catstable, grts_ext, grts_origin))

tar_read(mapnames)
tar_read(catstable)
ml <- tar_read(maps)

ml[[1]]
terra::plot(ml[[1]])
terra::values(ml[[1]], row = 5000, nrows = 1)
terra::coltab(ml[[1]])
terra::cats(ml[[1]])
terra::datatype(ml[[1]])

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

tms <- targets::tar_read(temporal_map_strata)
tms
terra::plot(tms)
terra::activeCat(tms) <- "stable"
terra::plot(tms)


# develop
targets::tar_load_globals()
tar_load(names = c(mapnames, catstable, temporal_map))
debug(add_changecats_tempstrat)
test <- add_changecats_tempstrat(
  tempstrat = temporal_map, cats = catstable, mapnames = mapnames
)

targets::tar_load_globals()
targets::tar_workspace("temporal_map_strata")
debugonce(binary_change)
test <- add_changecats_tempstrat(
  tempstrat = temporal_map, cats = catstable, mapnames = mapnames
)

