library(targets)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
Sys.setenv(TAR_PROJECT = "validation_sample")


prelabeled <- tar_read(prelabeled_validation_polygons_50)
rasterlabeled <- tar_read(raster_labeled_validation_polygons_50)

flea_data <- "C:/R/GitRepositories/flea-data"
fs::dir_create(
  file.path(flea_data, "validation")
)

prelabeled_sf <- prelabeled |>
  lapply(st_as_sf)

rasterlabeled_sf <- rasterlabeled |>
  lapply(st_as_sf)

prelabeled_sf |>
  purrr::walk(
    function(x) {
      if (nrow(x) == 0) return(invisible(NULL))

      layername <- unique(x$stratum_name)
      changecats <- unique(x$changecat)

      for (i in changecats) {

        xi <- x |>
          dplyr::filter(changecat == i)
        write_sf(
          obj = xi,
          dsn = file.path(
            flea_data, "validation",
            paste0(gsub("-", "", Sys.Date()), "_test_validation_sample.gpkg")),
          layer = paste0(layername, "_", i),
          delete_layer = TRUE)
      }
    }
  )

rasterlabeled_sf |>
  purrr::walk(
    function(x) {
      if (nrow(x) == 0) return(invisible(NULL))

      layername <- unique(x$stratum_name)
      layername <- gsub("bc_", "", layername)
      layername <- gsub("_changecat", "", layername)
      changecats <- unique(x$changecat)

      for (i in changecats) {

        xi <- x |>
          dplyr::filter(changecat == i) |>
          rename(label = label_label) |>
          relocate(
            label,
            matches("label_\\d{4}$"),
            ends_with("_2")
          )

        write_sf(
          obj = xi,
          dsn = file.path(
            flea_data, "validation",
            paste0(gsub("-", "", Sys.Date()), "_test_validation_sample.gpkg")),
          layer = paste0("cells_", layername, "_", i),
          delete_layer = TRUE)
      }
    }
  )

# enkel validatiepolygoon contour
# toevoegen aan geopackage
validation_poly_50 <- tar_read(name = validation_polygons_50)
validation_poly_50_sf <- validation_poly_50 |>
  lapply(st_as_sf)

validation_poly_50_sf |>
  purrr::walk(
    function(x) {
      if (nrow(x) == 0) return(invisible(NULL))

      layername <- unique(x$stratum_name)
      layername <- paste0("square50m_", stringr::str_extract(layername, "\\d+"))
      changecats <- unique(x$changecat)

      for (i in changecats) {

        xi <- x |>
          dplyr::filter(changecat == i)
        write_sf(
          obj = xi,
          dsn = file.path(
            flea_data, "validation",
            paste0(gsub("-", "", Sys.Date()), "_test_validation_sample.gpkg")),
          layer = paste0(layername, "_", i),
          delete_layer = TRUE)
      }
    }
  )

# sample als puntenlaag toevoegen
vs <- tar_read(name = validation_sample)

vs |>
  vect() |>
  st_as_sf() |>
  write_sf(
    dsn = file.path(
      flea_data, "validation",
      paste0(gsub("-", "", Sys.Date()), "_test_validation_sample.gpkg")),
    layer = "validation_sample_all",
    delete_layer = TRUE
  )

# settlement masks wegschrijven
sett <- tar_read(settlement_masks)
sett <- lapply(
  sett,
  \(x) subst(x, 1, 1, NA)
)
lapply(
  sett,
  \(x) {
    year <- stringr::str_extract(names(x), "\\d{4}")
    writeRaster(
      x = x,
      filename = file.path(
        flea_data, "validation",
        paste0("settlements_", year, ".tiff")
      ),
      overwrite = TRUE
    )
  }
)


# 30 polygonen voor test
# Settlement: 104 (low green) Loss: 6
# Bos: 400 (Forest and woodland) Loss: 6
# Landbouw: 200 (cropland), Gain: 6
# Natuur: 500 (heathland and shrub), Gain: 6
# water: 900 (lakes and reservoirs), Gain:  6

pl <- vect(prelabeled) |>
  st_as_sf() |>
  filter(
    (stratum_name == "bc_200_changecat" & changecat == "Gain") |
      (stratum_name == "bc_104_changecat" & changecat == "Loss") |
      (stratum_name == "bc_400_changecat" & changecat == "Loss") |
      (stratum_name == "bc_500_changecat" & changecat == "Gain") |
      (stratum_name == "bc_900_changecat" & changecat == "Gain")
  )

set.seed(42)
selectie <- pl %>%
  mutate(
    area = st_area(geometry),
    is_labeled = labels != "other-other-other"
  ) |>
  st_drop_geometry() |>
  group_by(stratum_name, changecat, grts_rank) |>
  summarise(
    area_labeled = sum(area[is_labeled]),
    n_polys = n(),
    n_labels = n_distinct(labels)
  ) |>
  arrange(
    stratum_name, changecat, area_labeled
  ) |>
  mutate(
    jump = floor(n() / 6),
    start = sample(seq_len(jump[1]), size = 1)
  ) |>
  slice(
    seq(start[1], start[1] + 5 * jump[1], jump[1])
  )

pl |>
  inner_join(selectie) |>
  select(-jump, -start) |>
  write_sf(
    dsn = file.path(
      flea_data, "validation",
      "test_validation_sample.gpkg"),
    layer = "testset_30_validaties",
    delete_layer = TRUE
  )


vect(rasterlabeled) |>
  st_as_sf() |>
  inner_join(selectie) |>
  write_sf(
    dsn = file.path(
      flea_data, "validation",
      paste0(gsub("-", "", Sys.Date()), "_test_validation_sample.gpkg")),
    layer = "test_cells_30_validaties",
    delete_layer = TRUE
  )


plotjes <- pl |>
  inner_join(selectie) |>
  group_by(grts_rank, stratum_name) |>
  tidyr::nest() |>
  mutate(
    plot = purrr::map2(data, stratum_name, \(x, y) {
      x |> ggplot() +
        geom_sf(aes(fill = labels)) +
        labs(title = y)
    })
  )

plotjes$plot


pl |>
  inner_join(selectie) |>
  select(-jump, -start) |>
  st_drop_geometry() |>
  distinct(grts_rank, stratum_name, changecat) |>
  group_by(stratum_name, changecat) |>
  summarise(
    grts_ranks = xfun::join_words(grts_rank, and = " en ", oxford_comma = FALSE)
  ) |>
  knitr::kable()

