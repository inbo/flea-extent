library(terra)

flea_data <- gsub(
  pattern = "flea-extent", replacement = "flea-data", x = here::here())

lg2013 <- rast(file.path(flea_data, "data", "2013", "LG2013_finaal_update.tif"))


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

lg2013 <- as.factor(lg2013)
levels(lg2013) <- catstable
coltab(lg2013) <- catstable |> dplyr::select(value, color) |> as.data.frame()
lg2013
plot(lg2013)

lg2013_selectie <- crop(lg2013, ext(200000, 205000, 200000, 205000))
plot(lg2013_selectie)

stratification <- focal(lg2013_selectie, w = 9, fun = "modal")
stratification <- as.factor(stratification)
levels(stratification) <- catstable
coltab(stratification) <- catstable |> dplyr::select(value, color) |> as.data.frame()

plot(stratification)




