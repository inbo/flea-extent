# Load packages required to define the pipeline:
library(targets)
library(geotargets)
# Load other packages as needed.

# Set target options:
controller <- crew::crew_controller_local(
  workers = 4,
  seconds_idle = 60,
  options_metrics = crew::crew_options_metrics(
    path = "worker_log_directory/",
    seconds_interval = 1
  )
)

if (tar_active()) {
  controller$start()
  autometric::log_start(
    path = "log.txt",
    seconds = 1,
    pids = controller$pids()
  )
}


tar_option_set(
  packages = c("tibble", "geotargets", "assertthat", "terra", "dplyr", "sf"),
  format = "qs",
  memory = "transient",
  garbage_collection = TRUE,
  controller = controller,
  storage = "worker",
  retrieval = "worker"
  #
)

targets_project_dir <- rprojroot::find_root(rprojroot::is_git_root) |>
  file.path("source/pipelines/")

tar_config_set(
  script = file.path(targets_project_dir, "script_validation_sample.R"),
  store = file.path(targets_project_dir, "store_validation_sample"),
  config = "_targets.yaml",
  project = "validation_sample",
  use_crew = TRUE
)

# Run the R scripts in the r/ folder with your custom functions:
tar_source(
  files = file.path(targets_project_dir, "r")
)
# tar_source("other_functions.R") # Source other scripts as needed # nolint

# inputs
git_root <- rprojroot::find_root(rprojroot::is_git_root)
flea_data <- gsub(
  pattern = "flea-extent", replacement = "flea-data", x = git_root
)
input_names <- c(
  "ecosysteemkaart_niv1_2016_v13",
  "ecosysteemkaart_niv1_2019_v13",
  "ecosysteemkaart_niv1_2022_v13"
)
bronnenkaarten <- c(
  "bronnenkaart_2016_v13",
  "bronnenkaart_2019_v13",
  "bronnenkaart_2022_v13"
)


input_years <- c(2016, 2019, 2022)
# nolint start
path_to_gdb <- "Z:/Projects/PRJ_FLEA/flea_output.gdb"
path_to_lbg <- "Z:/Projects/PRJ_FLEA/landbouwdata.gdb"
path_to_lyr <- "Z:/Projects/PRJ_FLEA/ecosysteemkaart_niv1_v12.lyr"
path_to_grts <- file.path(flea_data, "data/c-mon/flea_cmon_level15.tiff")
path_to_sources_lyr <- "Z:/Projects/PRJ_FLEA/bronnenkaart_v13.lyr"
# nolint end

layers_ruimtebeslag <- c(
  "GRB:WBN",
  "GRB:SBN",
  "GRB:GBG",
  "GRB:GBA",
  "GRB:KNW",
  "GRB:TRN"
)

# opm: GRB:WLAS zijn lijnobjecten
layers_water <- c(
  "GRB:WTZ"
)

layers_perceelgrens <- c(
  "GRB:ADP"
)

settlement_codes <- 101:106

# lbg mapping codes (for static branching)
lbg_mapping_flea_codes <- data.frame(
  flea_id = c(101, 104, 200, 300, 400, 500, 900)
)

habitatmap_terr_mapping_flea_codes <- data.frame(# nolint: object_length_linter.
  flea_id = c(500)
)


# to be changed later: download the raster files from zenodo

# target list:
sample_selection <- list(
  # names of maps
  tar_target(
    name = mapnames,
    command = input_maps(
      names = input_names
    )
  ),
  # get colour codes and class labels
  tar_target(
    name = catstable,
    command = get_lyrinfo(lyr = path_to_lyr)
  ),
  # read raster files
  tar_terra_rast(
    name = maps,
    command = get_map(
      gdb = path_to_gdb,
      name = mapnames,
      cats = catstable,
      origin = grts_origin,
      grts = fleagrts
    ),
    pattern = map(mapnames),
    preserve_metadata = "zip"
  ),
  tar_target(
    name = sourcenames,
    command = input_maps(
      names = bronnenkaarten
    ),
    description =
      "Name of the rasterlayers encoding the source used to
    assign an ecosystem type"
  ),
  tar_target(
    name = sourcestable,
    command = get_lyrinfo(lyr = path_to_sources_lyr),
    description = "categorical mapping for sources"
  ),
  tar_terra_rast(
    name = maps_sources,
    command = get_map(
      gdb = path_to_gdb,
      name = sourcenames,
      cats = sourcestable,
      #same as grts_origin
      origin = c(0, 0),
      grts = fleagrts
    ),
    pattern = map(sourcenames),
    preserve_metadata = "zip"
  ),
  # calculate masks
  tar_terra_rast(
    name = settlement_masks,
    command = calc_mask(maps = maps, values = settlement_codes),
    pattern = map(maps)
  ),
  tar_terra_rast(
    name = fleagrts,
    command = get_grts(
      path = path_to_grts
    )
  ),
  tar_target(
    name = grts_origin,
    command = terra::origin(fleagrts)
  ),
  # create (temporal) difference maps
  tar_terra_rast(
    name = temporal_map,
    command = create_temporal_maps(
      input_maps = maps
    ),
    preserve_metadata = "zip"
  ),
  # add change categories to temporal map
  tar_terra_rast(
    name = temporal_map_strata,
    command = add_changecats_tempstrat(
      tempstrat = temporal_map,
      cats = catstable,
      mapnames = mapnames
    ),
    preserve_metadata = "zip"
  ),
  # get names of columns with change categories for each land use type
  tar_target(
    name = lu_changecats,
    command = get_changecat_columns(tempstrat = temporal_map_strata)
  ),
  geotargets::tar_terra_rast(
    name = separate_grts,
    command = separate_grts_strata(
      stratum_raster = temporal_map_strata,
      fleagrts = fleagrts,
      stratum_name = lu_changecats
    ),
    pattern = map(lu_changecats),
    preserve_metadata = "zip",
    deployment = "main",
    garbage_collection = TRUE
  ),
  # extract grts sample for each land use change category
  geotargets::tar_terra_vect(
    name = validation_sample,
    command = extract_sample(
      separate_grts = separate_grts,
      stratum_name = lu_changecats,
      ntot = 40 * 2 * 4,
      nmin = 40,
      min_stratum_size = 1000 # 10 ha
    ),
    pattern = map(separate_grts, lu_changecats),
    deployment = "main",
    memory = "transient",
    garbage_collection = TRUE
  ),
  # convert sampling locations to square polygons of size 9 x 9
  geotargets::tar_terra_vect(
    name = validation_polygons,
    command = point_to_gridcell(
      xy = validation_sample,
      cell_width_m = 90,
      point_position = "center",
      crs = 31370
    ),
    pattern = map(validation_sample)
  ),
  geotargets::tar_terra_vect(
    name = validation_polygons_50,
    command = point_to_gridcell(
      xy = validation_sample,
      cell_width_m = 50,
      point_position = "center",
      crs = 31370
    ),
    pattern = map(validation_sample)
  )
)

prelabelling_sources <- list(
  # add grb waterways
  tar_target(
    name = lyrs_waterways,
    command = read_layernames(x = layers_water),
    description = "Name of GRB waterway layers"
  ),
  geotargets::tar_terra_vect(
    name = grb_waterways,
    command = get_grb_by_row(
      layer = lyrs_waterways,
      polygons = validation_polygons
    ),
    pattern = cross(lyrs_waterways, validation_polygons),
    filetype = "GPKG",
    description = "GRB waterways vector data for each validation polygon"
  ),
  # add grb settlements
  tar_target(
    name = lyrs_settlements,
    command = read_layernames(x = layers_ruimtebeslag),
    description = "Name of GRB settlement layers"
  ),
  geotargets::tar_terra_vect(
    name = grb_settlements,
    command = get_grb_by_row(
      layer = lyrs_settlements,
      polygons = validation_polygons
    ),
    pattern = cross(lyrs_settlements, validation_polygons),
    filetype = "GPKG",
    description = "GRB settlements vector data for each validation polygon"
  ),
  # Landbouwgebruikspercelen data
  targets::tar_target(
    name = lbg_mapping_df,
    command = mapping_lbg_to_flea(),
    description = "Mapping between GWSCOD_H, GWSNAME_H and flea landuse codes"
  ),
  targets::tar_target(
    name = lbg_layers,
    command = get_lbg_layernames(path_to_lbg),
    description = "Name of LBG layernames"
  ),
  tarchetypes::tar_map(
    values = lbg_mapping_flea_codes,
    names = "flea_id",
    # get the unique GWSCOD_H that map to flea_ids
    targets::tar_target(
      #The base name (becomes mapping_101, mapping_104, etc.)
      name = lbg_mapping,
      command = lbg_mapping_df$GWSCOD_H[
        lbg_mapping_df$FLEA == flea_id & !is.na(lbg_mapping_df$FLEA)
      ],
      description = "GWSCOD_H mapping to flea code"
    ),
    # extract them from the LBG layers
    geotargets::tar_terra_vect(
      name = lbg, # Becomes lbg_101, etc.
      command = get_lbg(
        path_to_lbg = path_to_lbg,
        layer = lbg_layers,
        from_fields = c("GWSCOD_H", "GWSNAM_H", "STAT_BGV"),
        where_field = "GWSCOD_H",
        # referencing 'lbg_mapping' here automatically resolves to
        # 'lbg_mapping_101'
        # because they are in the same tar_map scope!
        where_values = lbg_mapping,
        flea_value = flea_id
      ),
      # This applies dynamic branching to every static branch
      # in this case the lbg_layer for each year
      pattern = map(lbg_layers),
      description = "LBG vector data"
    ),
    # --- Step 3: Crop (Dynamic Branching: cross) ---
    # Result: lbg_cropped_101, lbg_cropped_104...
    geotargets::tar_terra_vect(
      name = lbg_cropped,
      command = {
        sc <- spatvector_crop(x = lbg, y = validation_polygons)
        sc$year_flea <- stringr::str_extract(sc$layer, "\\d{4}") |> as.numeric()
        sc
      },
      # 'lbg' here refers to lbg_101 (which is already branched).
      # 'cross' will multiply lbg_101 branches by validation_polygons branches.
      pattern = cross(lbg, validation_polygons),
      description = "LBG vector data cropped to validation polygons"
    )
  ),
  # settlements processing
  geotargets::tar_terra_vect(
    name = grb_settlements_processed,
    command = process_settlement(
      grb = grb_settlements
    ),
    description = "Cover and assign codes 101 or 102 to GRB settlements.
    101 = (GBG > GBA > KNW) > 102 = (WBN > SBN > TRN verhard)"
  ),
  # GRB WTZ processing
  geotargets::tar_terra_vect(
    name = grb_water_wtz_processed,
    command = process_water_wtz(
      grb_wtz = grb_waterways
    )
  ),
  targets::tar_target(
    name = watersurfaces_meta,
    command = data.frame(
      doi = c(
        "10.5281/zenodo.3386859",
        "10.5281/zenodo.4117543",
        "10.5281/zenodo.14203168"
      ),
      source = c(
        "ortho_2015_2016",
        "ortho_2018_2019",
        "ortho_2021_2023"
      ),
      year_flea = c(2016, 2019, 2022),
      version = c(
        "v1.0", "v1.1", "v2024"
      )
    )
  ),
  targets::tar_target(
    name = habitatmap_terr_meta,
    command = data.frame(
      doi = c(
        "10.5281/zenodo.3540740",
        "10.5281/zenodo.13861955",
        "10.5281/zenodo.13886579"
      ),
      version = c(
        "habitatmap_terr_2018_v2",
        "habitatmap_terr_2020_v2",
        "habitatmap_terr_2023_v1"
      ),
      year_flea = c(2016, 2019, 2022)
    )
  ),
  targets::tar_target(
    name = zenodo_habitatmap_terr,
    command = download_habitatmap_terr(
      path_flea_data = flea_data,
      meta = habitatmap_terr_meta
    ),
    pattern = map(habitatmap_terr_meta)
  ),
  targets::tar_target(
    name = zenodo_watersurface,
    command = download_watersurfaces(
      path_flea_data = flea_data,
      meta = watersurfaces_meta
    ),
    pattern = map(watersurfaces_meta)
  ),
  # read INBO watersurfaces maps and crop with validation polygons
  geotargets::tar_terra_vect(
    name = watersurfaces_processed,
    command = get_watersurfaces(
      path_version = zenodo_watersurface,
      polygons = validation_polygons,
      meta = watersurfaces_meta
    ),
    pattern = cross(
      map(zenodo_watersurface, watersurfaces_meta),
      validation_polygons
    )
  ),
  tarchetypes::tar_map(
    values = habitatmap_terr_mapping_flea_codes,
    names = "flea_id",
    targets::tar_target(
      name = types,
      command = get_types(flea_id = flea_id)
    ),
    # read INBO habitatmap_terr and crop with validation polygons
    geotargets::tar_terra_vect(
      name = terr_cropped,
      command = get_habitatmap_terr(
        path_version = zenodo_habitatmap_terr,
        polygons = validation_polygons,
        meta = habitatmap_terr_meta,
        types = types,
        min_phab = 50, #selecteert alle eerste eenheden (tweede eenheid max 30%)
        flea_value = flea_id
      ),
      pattern = cross(
        map(zenodo_habitatmap_terr, habitatmap_terr_meta),
        validation_polygons
      )
    )
  )
)

overlay_prelabel_sources <- list(
  # combine the GRB water layer with the INBO watersurfaces
  geotargets::tar_terra_vect(
    name = vp_water,
    command = combine_grb_inbo_water(
      grb_water = grb_water_wtz_processed,
      inbo_water = watersurfaces_processed,
      meta = watersurfaces_meta
    ),
    pattern = map(watersurfaces_meta)
  ),
  geotargets::tar_terra_vect(
    name = vp_water_grb_lbg_terr,
    command = combine_sources(
      water = vp_water,
      settlements = grb_settlements_processed,
      lbg_101 = lbg_cropped_101,
      lbg_104 = lbg_cropped_104,
      lbg_200 = lbg_cropped_200,
      lbg_300 = lbg_cropped_300,
      lbg_400 = lbg_cropped_400,
      lbg_500 = lbg_cropped_500,
      lbg_900 = lbg_cropped_900,
      terr_500 = terr_cropped_500,
      polygons = validation_polygons
    ),
    pattern = map(vp_water)
  ),
  geotargets::tar_terra_vect(
    name = vp_water_grb_lbg_terr_cleaned,
    command = postprocess_prelabelling(
      water_grb_lbg = vp_water_grb_lbg_terr,
      settlement_mask = settlement_masks
    ),
    pattern = map(vp_water_grb_lbg_terr, settlement_masks)
  )
  ,
  geotargets::tar_terra_vect(
    name = vp_water_grb_lbg_terr_singletarget,
    command = single_wsp(vp_water_grb_lbg_terr_cleaned)
  ),
  geotargets::tar_terra_vect(
    name = prelabeled_validation_polygons,
    command = intersect_validation_polygons(
      wsp_target = vp_water_grb_lbg_terr_singletarget,
      lu_changecat = lu_changecats,
      input_years = input_years,
      mmu = 30
    ),
    pattern = map(lu_changecats)
  ),
  geotargets::tar_terra_vect(
    name = prelabeled_validation_polygons_50,
    command = crop_labeled_polygons(
      pvp = prelabeled_validation_polygons,
      crop_with = validation_polygons_50
    ),
    pattern = map(prelabeled_validation_polygons, validation_polygons_50)
  ),
  geotargets::tar_terra_vect(
    name = raster_labeled_validation_polygons_50,
    command = raster_to_polygons(
      tm = temporal_map_strata, # temporal raster
      vp = validation_polygons_50, # validation polygons
      input_years = as.character(input_years)
    ),
    pattern = map(validation_polygons_50),
    deployment = "main"
  )
)

c(
  sample_selection,
  prelabelling_sources,
  overlay_prelabel_sources
)
