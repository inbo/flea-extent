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
    dplyr::select(-alpha) |>
    dplyr::filter(value != 0)
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
  # replace 0 with NA
  map <- mask(map, map, maskvalues = 0)
  setMinMax(map)
  map <- apply_cats(x = map, cats = cats, name = name)
  origin(map) <- origin
  map <- extend(map, grts)
  crs(map) <- crs(grts)
  map <- writeRaster(
    map,
    tempfile(fileext = ".tif"),
    datatype = "INT2U", overwrite = TRUE)
  return(map)
}



get_grts <- function(path) {
  assert_that(file.exists(path))
  fleagrts <- rast(path)
  return(fleagrts)
}


binary_change <- function(data, lg_values, mapnames) {
  binary <- vector("list", length = length(lg_values))
  bc_colnames <- paste0("bc_", lg_values)
  binary <- setNames(binary, bc_colnames)
  colselect <- paste0("value_", mapnames)
  for (i in seq_along(lg_values)) {
    binary[[bc_colnames[i]]] <-
      purrr::map(colselect, ~{
        stringr::str_detect(
          data[[.x]],
          paste0("^", lg_values[i], "$")
        ) %>% as.numeric()
      }) %>%
      purrr::list_transpose() %>%
      purrr::map_chr(
        .f = \(x) paste(x, collapse = "")
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

create_temporal_maps <- function(input_maps) {
  namesvec <- purrr::map_vec(input_maps, names)
  map_stack <- rast(input_maps)
  names(map_stack) <- namesvec
  maskmap <- any(is.na(map_stack))
  map_stack <- mask(map_stack, maskmap, maskvalues = 1)
  temporal_stratification <- unique(
    map_stack * 1,
    as.raster = TRUE
  )
  return(temporal_stratification)
}

add_changecats_tempstrat <- function(tempstrat, cats, mapnames) {

  lg_values <- as.character(cats$value)

  additional_levels <- freq(tempstrat) %>%
    as_tibble() %>%
    tidyr::separate(
      value,
      into = mapnames,
      sep = "_",
      remove = FALSE
    ) %>%
    tidyr::pivot_longer(
      cols = all_of(mapnames),
      names_to = "mapname",
      values_to = "year_value"
    )

  additional_levels <- additional_levels %>%
    left_join(
      cats %>%
        mutate(
          value = as.character(value),
          year_label = label,
          .keep = "none"
        ),
      by = join_by(
         year_value == value
      )
    )

  additional_levels <- additional_levels %>%
    tidyr::pivot_wider(
      id_cols = c(layer, value, count),
      names_from = mapname,
      values_from = c(year_value, year_label),
      names_sort = TRUE,
      names_glue = "{gsub('year_','',.value)}_{mapname}"
    )

  bc_colnames <- paste0("bc_", lg_values)

  additional_levels <- additional_levels %>%
    binary_change(lg_values = lg_values, mapnames = mapnames) %>%
    rowwise() %>%
    mutate(stable = all(
      c_across(starts_with("value_")) == first(c_across(starts_with("value_")))
        ) %>%
        if_else("stable", "changed") %>%
        as.factor()
      ) %>%
    ungroup() %>%
    mutate(
      across(
        all_of(bc_colnames),
        \(x) categorize_land_use_change(x),
        .names = "{.col}_changecat"
      )
    )

  join_levels <- cats(tempstrat)[[1]] %>%
    mutate(across(all_of(mapnames), as.character)) %>%
    as_tibble() %>%
    inner_join(
      additional_levels,
      by = join_by(label == value)
    ) %>%
    select(-starts_with("value_"))
  levels(tempstrat) <- join_levels
  coltab(tempstrat) <- NULL

  return(tempstrat)
}


get_changecat_columns <- function(tempstrat) {
  changecat_columns <- names(cats(tempstrat)[[1]])
  changecat_columns <- changecat_columns[
    stringr::str_detect(changecat_columns, "_changecat$")
  ]
  return(changecat_columns)
}


calc_mask <- function(maps, values) {
  my_mask <- maps %in% values
  return(my_mask)
}


