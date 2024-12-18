input_maps <- function(names) {
  assert_that(is.character(names))
  return(names)
}

get_lyrinfo <- function(lyr) {
  assert_that(file.exists(lyr))

  output_qml <- tempfile(fileext = ".qml")
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
    ) |>
    dplyr::select(-alpha)
  return(catstable)
}

apply_cats <- function(x, cats, name) {
  assert_that(is.data.frame(cats))
  assert_that(all(c("label", "value", "color") %in% names(cats)))
  xc <- as.factor(x)
  names(cats)[names(cats) == "label"] <- name
  levels(xc) <- cbind(cats, t(col2rgb(cats$color)))
  names(xc) <- name
  return(xc)
}


get_map <- function(gdb, name, cats, origin, grts) {
  assert_that(file.exists(gdb))
  assert_that(is.string(name))
  assert_that(is.numeric(origin) && length(origin) == 2)
  map <- rast(gdb, subds = name)
  setMinMax(map)
  map <- apply_cats(x = map, cats = cats, name = name)
  NAflag(map) <- 0
  map <- writeRaster(
    map,
    tempfile(fileext = ".tif"),
    datatype = "INT2U", overwrite = TRUE)
  origin(map) <- origin
  map <- extend(map, grts)
  return(map)
}



get_grts <- function(path) {
  assert_that(file.exists(path))
  fleagrts <- rast(path)
  return(fleagrts)
}


binary_change <- function(data, lg) {
  binary <- vector("list", length = length(lg))
  binary <- setNames(binary, lg)
  for (i in lg) {
    binary[[i]] <- paste0(
      stringr::str_detect(data$lg2013_label, i) %>% as.numeric(),
      stringr::str_detect(data$lg2016_label, i) %>% as.numeric(),
      stringr::str_detect(data$lg2019_label, i) %>% as.numeric()
    )
  }
  bind_cols(data, binary)
}

categorize_land_use_change <- function(b, simple = TRUE) {
  if (simple) {
    case_when(
      # Stable conditions
      grepl("^0+$", b) ~ "Stable absence",
      grepl("^1+$", b) ~ "Stable presence",

      # Simple changes
      grepl("^0+1+$", b) ~ "Gain",
      grepl("^1+0+$", b) ~ "Loss",

      # Default case
      TRUE ~ "Other complex pattern"
    )
  } else {
    case_when(
      # Stable conditions
      grepl("^0+$", b) ~ "Stable absence",
      grepl("^1+$", b) ~ "Stable presence",

      # Simple changes
      grepl("^0+1+$", b) ~ "Gain",
      grepl("^1+0+$", b) ~ "Loss",

      # Complex changes
      grepl("^0+1+0+$", b) ~ "Temporary gain",
      grepl("^1+0+1+$", b) ~ "Temporary loss",
      grepl("^0+1+0+1+$", b) ~ "Intermittent presence (starting absent)",
      grepl("^1+0+1+0+$", b) ~ "Intermittent presence (starting present)",

      # Oscillating changes
      grepl("^(01)+0?$", b) ~ "Oscillating (starting absent)",
      grepl("^(10)+1?$", b) ~ "Oscillating (starting present)",

      # Other complex patterns
      grepl("01.*1$", b) & !grepl("^0+1+$", b) ~ "Complex gain",
      grepl("10.*0$", b) & !grepl("^1+0+$", b) ~ "Complex loss",

      # Default case
      TRUE ~ "Other complex pattern"
    )
  }
}

create_temporal_maps <- function(input_maps, cats) {
  map_stack <- rast(input_maps)
  temporal_stratification <- unique(
    map_stack * 1,
    as.raster = TRUE
  )
  # deal with NAs
  temporal_stratification <- mask(temporal_stratification, map_stack)
  temporal_stratification <- droplevels(temporal_stratification)

  lg <- gsub(pattern = "^\\d\\s-\\s", replacement = "", x = cats$label)

  additional_levels <- freq(temporal_stratification) %>%
    as_tibble() %>%
    tidyr::separate(
      value,
      into = c("lg2013", "lg2016", "lg2019"),
      sep = "_",
      remove = FALSE
    ) %>%
    left_join(
      cats %>%
        mutate(
          value = as.character(value),
          lg2013_label = label,
          .keep = "none"
        ),
      by = join_by(lg2013 == value)
    ) %>%
    left_join(
      cats %>%
        mutate(
          value = as.character(value),
          lg2016_label = label,
          .keep = "none"
        ),
      by = join_by(lg2016 == value)
    ) %>%
    left_join(
      catstable %>%
        mutate(
          value = as.character(value),
          lg2019_label = label,
          .keep = "none"
        ),
      by = join_by(lg2019 == value)
    ) %>%
    binary_change(lg = lg) %>%
    rowwise() %>%
    mutate(stable = ifelse(
      all(lg2013 == lg2016, lg2016 == lg2019),
      "stable", "changed"
    ) %>%
      as.factor()) %>%
    ungroup() %>%
    mutate(
      across(
        all_of(lg),
        \(x) categorize_land_use_change(x),
        .names = "{.col}_changecat"
      )
    )

  join_levels <- cats(temporal_stratification)[[1]] %>%
    mutate(across(starts_with("lg"), as.character)) %>%
    inner_join(
      additional_levels,
      by = join_by(lg2013, lg2016, lg2019, label == value)
    )
  levels(temporal_stratification) <- join_levels
  coltab(temporal_stratification) <- NULL


  return(temporal_stratification)
}
