library(targets)

Sys.setenv(TAR_PROJECT = "validation_sample")

tar_make()

# debugging and inspection

targets::tar_meta(
  fields = error,
  complete_only = TRUE
)

targets::tar_visnetwork(label = c("description", "time", "size"))


# logging
library(autometric)
log_file <- "log.txt"
log_data <- log_read(log_file)
log_plot(log_data, metric = "resident")

#tar_load_globals()
#tar_load(names = c(mapnames, catstable, grts_ext, grts_origin))

tar_read(mapnames)
tar_read(catstable) |> tail()
ml <- tar_read(maps)

ml[[1]]
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

cc <- targets::tar_read(changecats)

# develop
targets::tar_load_globals()
tar_load(names = c(separate_grts))
debugonce(get_changecats)
get_changecats(separate_grts)

targets::tar_load_globals()
targets::tar_workspace("separate_grts_03c5fe21ea8598b0")
debugonce(separate_grts_strata)
test <- separate_grts_strata(
    stratum_raster = temporal_map_strata,
    fleagrts = fleagrts,
    stratum_name = lu_changecats)




