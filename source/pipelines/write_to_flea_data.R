library(targets)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
Sys.setenv(TAR_PROJECT = "validation_sample")


prelabeled <- tar_read(prelabeled_validation_polygons_50)
rasterlabeled <- tar_read(raster_labeled_validation_polygons_50)

flea_data <- "C:/R/GitRepositories/flea-data" # nolint: absolute_path_linter.
flea_validation <- fs::dir_create(
  file.path(flea_data, "validation_test2")
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
            flea_validation,
            paste0(gsub("-", "", Sys.Date()), "_test_validation_sample.gpkg")
          ),
          layer = paste0(layername, "_", i),
          delete_layer = TRUE
        )
      }
    }
  )

# uit eerste test bleek dat dit niet goed werkte
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
            flea_validation,
            paste0(gsub("-", "", Sys.Date()), "_test_validation_sample.gpkg")
          ),
          layer = paste0("cells_", layername, "_", i),
          delete_layer = TRUE
        )
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
            flea_validation,
            paste0(gsub("-", "", Sys.Date()), "_test_validation_sample.gpkg")
          ),
          layer = paste0(layername, "_", i),
          delete_layer = TRUE
        )
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
      flea_validation,
      paste0(gsub("-", "", Sys.Date()), "_test_validation_sample.gpkg")
    ),
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
        flea_validation,
        paste0("settlements_", year, ".tiff")
      ),
      overwrite = TRUE
    )
  }
)

# bronnenkaarten wegschrijven
write_categorical_raster <- function(
  x, filename, datatype = "INT1U", overwrite = TRUE
) {
  writeRaster(
    x, filename = filename, datatype = datatype, overwrite = overwrite
  )

  rat <- cats(x)[[1]]
  if (is.null(rat)) return(invisible(filename))

  label_col <- setdiff(
    names(rat), c("value", "color", "red", "green", "blue")
  )[1]

  # Prepend a transparent/nodata entry for position 0
  color_entries <- c(
    "      <Entry c1='0' c2='0' c3='0' c4='0'/>",
    mapply(
      function(r, g, b) {
        sprintf("      <Entry c1='%d' c2='%d' c3='%d' c4='255'/>",
                r, g, b)
      },
      rat$red, rat$green, rat$blue,
      SIMPLIFY = TRUE
    )
  )

  rat_rows <- mapply(
    function(i, val, lbl) {
      sprintf('        <Row index="%d"><F>%d</F><F>%s</F></Row>',
              i,        # 0-based row position
              val,      # pixel value in the Value field
              xml_escape(lbl))
    },
    seq_len(nrow(rat)) - 1L,   # i: 0, 1, 2, ...
    rat$value,                  # val: 1, 2, 3, ...
    rat[[label_col]],
    SIMPLIFY = TRUE
  )
  xml <- paste0(
    "<PAMDataset>\n",
    "  <PAMRasterBand band='1'>\n",
    "    <GDALRasterAttributeTable tableType='thematic'>\n",
    paste0(
      "      <FieldDefn index='0'><Name>Value</Name>",
      "<Type>0</Type><Usage>0</Usage></FieldDefn>\n"
    ),
    "      <FieldDefn index='1'><Name>",
    label_col,
    "</Name><Type>2</Type><Usage>2</Usage></FieldDefn>\n",
    paste(rat_rows, collapse = "\n"), "\n",
    "    </GDALRasterAttributeTable>\n",
    "    <ColorTable palette='RGB'>\n",
    paste(color_entries, collapse = "\n"), "\n",
    "    </ColorTable>\n",
    "  </PAMRasterBand>\n",
    "</PAMDataset>\n"
  )

  writeLines(xml, paste0(filename, ".aux.xml"))
  write_qml(rat, filename, label_col)   # new: QGIS reads this first
  invisible(filename)
}

write_qml <- function(rat, filename, label_col) {
  items <- mapply(function(val, lbl, col) {
    sprintf(
      "        <paletteEntry alpha='255' color='%s' label='%s' value='%d'/>",
      col, html_escape(lbl), val
    )
  }, rat$value, rat[[label_col]], rat$color, SIMPLIFY = TRUE)

  qml <- paste0(
    "<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>\n",
    "<qgis version=\"3.0\" styleCategories=\"Symbology\">\n",
    "  <pipe>\n",
    "    <rasterrenderer band=\"1\" type=\"paletted\" opacity=\"1\">\n",
    "      <colorPalette>\n",
    paste(items, collapse = "\n"), "\n",
    "      </colorPalette>\n",
    "    </rasterrenderer>\n",
    "  </pipe>\n",
    "</qgis>\n"
  )

  writeLines(qml, paste0(tools::file_path_sans_ext(filename), ".qml"))
}

html_escape <- function(x) {
  x <- gsub("&",  "&amp;",  x)
  x <- gsub("<",  "&lt;",   x)
  x <- gsub(">",  "&gt;",   x)
  x <- gsub("\"", "&quot;", x)
  x
}

xml_escape <- function(x) {
  x <- gsub("&",  "&amp;",  x)
  x <- gsub("<",  "&lt;",   x)
  x <- gsub(">",  "&gt;",   x)
  x <- gsub("\"", "&quot;", x)
  x
}

maps_sources <- tar_read(maps_sources)
lapply(
  maps_sources,
  \(x) {
    year <- stringr::str_extract(names(x), "\\d{4}")
    write_categorical_raster(
      x = x,
      filename = file.path(
        flea_validation,
        paste0("sources_", year, ".tiff")
      ),
      datatype = "INT1U",
      overwrite = TRUE
    )
  }
)
maps_sources

plot(maps_sources$maps_sources_d16cafd084c7d512)


maps_v13 <- tar_read(maps)
plot(maps_v13$maps_0f833c2d4198addf)
names(maps_v13$maps_0f833c2d4198addf)
terra::freq(maps_v13$maps_0f833c2d4198addf)
terra::minmax(maps_v13$maps_0f833c2d4198addf)

lapply(
  maps_v13,
  \(x) {

    ct <- coltab(x)
    readr::write_csv(
      ct[[1]],
      file.path(
        flea_validation,
        paste0(names(x), ".csv")
      )
    )

    write_categorical_raster(
      x = x,
      filename = file.path(
        flea_validation,
        paste0(names(x), ".tif")
      ),
      datatype = "INT2U",
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


set.seed(20260707)
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
      flea_validation,
      paste0(gsub("-", "", Sys.Date()), "_test_validation_sample.gpkg")
    ),
    layer = "testset_30_validaties",
    delete_layer = TRUE
  )

# nolint start
# vect(rasterlabeled) |>
#   st_as_sf() |>
#   inner_join(selectie) |>
#   write_sf(
#     dsn = file.path(
#       flea_data, "validation",
#       paste0(gsub("-", "", Sys.Date()), "_test_validation_sample.gpkg")),
#     layer = "test_cells_30_validaties",
#     delete_layer = TRUE
#   )
# nolint end


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


selected_table_data <- pl |>
  inner_join(selectie) |>
  select(-jump, -start) |>
  st_drop_geometry() |>
  distinct(grts_rank, stratum_name, changecat)

selected_table_data |>
  group_by(stratum_name, changecat) |>
  summarise(
    grts_ranks = xfun::join_words(grts_rank, and = " en ", oxford_comma = FALSE)
  ) |>
  knitr::kable()

# nolint start
# # check for one
# plotjes |>
#   ungroup() |>
#   slice(1) |> pull(plot)
#
# sf_sample_3363 <- plotjes |>
#   ungroup() |>
#   slice(1) |>
#   inner_join(selected_table_data) |>
#   pull(data)
#
# sf_sample_3363[[1]] |> View()
#
# cellnr <- sf_sample_3363[[1]]$cell |> unique()
# tar_load(temporal_map_strata)
# temporal_map_strata[cellnr]
#
# library(ggplot2)
# sapply(rasterlabeled_sf, \(x) any(unique(x$cell) == cellnr))
# rasterlabeled_sf$raster_labeled_validation_polygons_50_544dab2123b30aa6 |>
#   filter(grts_rank == 3363) |>
#   ggplot() +
#   geom_sf() +
#   geom_sf_text(aes(label = label_label)) +
#   geom_sf(
#     data = sf_sample_3363[[1]],
#     aes(fill = labels), alpha = 0.2
#   )
#
# mapview::mapview(sf_sample_3363[[1]])
#
# tar_load(vp_water_grb_lbg_terr_singletarget)
#
# vp_water_grb_lbg_terr_singletarget |>
#   st_as_sf() |>
#   filter(stratum_name == "bc_104_changecat", changecat == "Loss", grts_rank == 3363)
# nolint end
