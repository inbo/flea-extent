library(terra)
git_root <- rprojroot::find_root(rprojroot::is_git_root)
flea_data <- gsub(
  pattern = "flea-extent", replacement = "flea-data", x = git_root
)

gdb <- "Z:/Projects/PRJ_FLEA/flea_data.gdb" # nolint: absolute_path_linter.
lyr <- "Z:/Projects/PRJ_FLEA/reclass_bwk2016.lyr"# nolint: absolute_path_linter.

file.exists(gdb)
file.exists(lyr)

bwk2016 <- rast(gdb, subds = "reclass_bwk2016")
bwk2020 <- rast(gdb, subds = "reclass_bwk2020")
bwk2023 <- rast(gdb, subds = "reclass_bwk2023")

setMinMax(bwk2016)
setMinMax(bwk2020)
setMinMax(bwk2023)

fs::dir_create(
  file.path(flea_data, "data", "2016_2020_2023")
)
output_qml <- file.path(
  flea_data, "data", "2016_2020_2023", "reclass_bwk2016.qml"
)
qgisprocess::qgis_run_algorithm(
  "slyr:lyrtoqml",
  INPUT = lyr,
  OUTPUT = output_qml
)

slyr <- xml2::read_xml(output_qml)

catstable <- xml2::xml_find_all(
  x = slyr, ".//pipe//rasterrenderer//colorPalette"
) |>
  xml2::xml_contents() |>
  purrr::map(xml2::xml_attrs) |>
  purrr::map_df(~ as.list(.)) |>
  dplyr::relocate(value, label) |>
  dplyr::mutate(
    value = as.numeric(value),
    color = toupper(color)
  )

apply_cats <- function(x, cats = catstable, name, coltab = TRUE) {
  xc <- as.factor(x)
  names(cats)[names(cats) == "label"] <- name
  levels(xc) <- cats
  if (coltab) {
    coltab(xc) <- cats |>
      dplyr::select(value, color) |> # nolint
      as.data.frame()
  }
  names(xc) <- name
  return(xc)
}

bwk2016 <- apply_cats(bwk2016, name = "bwk2016")
bwk2020 <- apply_cats(bwk2020, name = "bwk2020")
bwk2023 <- apply_cats(bwk2023, name = "bwk2023")

fleagrts <- rast(file.path(flea_data, "data/c-mon/flea_cmon_level15.tiff"))

# check same origin and resolution
origin(bwk2016)
origin(fleagrts)
res(bwk2016)
res(fleagrts)

# the difference in origin is very small (0.4 mm and 3.6 mm)
# we can just reassign the origin to c(0, 0)
# alternatively (and maybe safer), use resampling(),
# however currently that requires removing coltab() <- NULL

origin(bwk2016) <- origin(fleagrts)
origin(bwk2020) <- origin(fleagrts)
origin(bwk2023) <- origin(fleagrts)

# the value 0 should be assigned to NA
NAflag(bwk2016) <- 0
NAflag(bwk2020) <- 0
NAflag(bwk2023) <- 0

cats(bwk2016)
minmax(bwk2016)
is.factor(bwk2016)
plot(bwk2016)

# https://github.com/rspatial/terra/issues/1668
# workaround: remove colours before saving
coltab(bwk2016) <- NULL
writeRaster(
  bwk2016,
  filename = file.path(
    flea_data, "data", "2016_2020_2023",
    "bwk2016.tif"
  ),
  overwrite = TRUE
)

coltab(bwk2020) <- NULL
writeRaster(
  bwk2020,
  filename = file.path(
    flea_data, "data", "2016_2020_2023",
    "bwk2020.tif"
  ),
  overwrite = TRUE
)

coltab(bwk2023) <- NULL
writeRaster(
  bwk2023,
  filename = file.path(
    flea_data, "data", "2016_2020_2023",
    "bwk2023.tif"
  ),
  overwrite = TRUE
)

check <- rast(
  file.path(
    flea_data, "data", "2016_2020_2023",
    "bwk2016.tif"
  )
)

# re-apply colours (remove when / if issue is solved)
coltab(check) <- as.data.frame(catstable[, c("value", "color")])
plot(check)
minmax(check)
