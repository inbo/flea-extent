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
terra::ext(grts)
terra::ext(ml[[1]])

