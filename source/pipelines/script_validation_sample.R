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
  error = "null",
  memory = "transient",
  garbage_collection = TRUE,
  controller = controller
  #
)

targets_project_dir <- rprojroot::find_root(rprojroot::is_git_root) |>
  file.path("source/pipelines/")

tar_config_set(
  script = file.path(targets_project_dir, "script_validation_sample.R"),
  store = file.path(targets_project_dir, "store_validation_sample"),
  config = "_targets.yaml",
  project = "validation_sample",
  use_crew = TRUE)

# Run the R scripts in the R/ folder with your custom functions:
tar_source(
  files = file.path(targets_project_dir, "R")
)
# tar_source("other_functions.R") # Source other scripts as needed # nolint

# inputs
git_root <- rprojroot::find_root(rprojroot::is_git_root)
flea_data <- gsub(
  pattern = "flea-extent", replacement = "flea-data", x = git_root
)
input_names <- c("reclass_bwk2016", "reclass_bwk2020", "reclass_bwk2023")
input_years <- c(2016, 2019, 2022)
path_to_gdb <- "Z:/Projects/PRJ_FLEA/flea_data.gdb"
path_to_lbg <- "Z:/Projects/PRJ_FLEA/landbouwdata.gdb"
path_to_lyr <- "Z:/Projects/PRJ_FLEA/reclass_bwk2016.lyr"
path_to_grts <- file.path(flea_data, "data/c-mon/flea_cmon_level15.tiff")

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


# to be changed later: download the raster files from zenodo

# target list:
list(
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
  # calculate masks
  tar_terra_rast(
    name = settlement_masks,
    command = calc_mask(maps = maps, values = c(101, 102, 105, 106)),
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
  # add grb waterways
  # add grb settlements
  # add grb parcel outlines
  tar_target(
    name = lyrs_waterways,
    command = read_layernames(x = layers_water)
  ),
  tar_target(
    name = lyrs_settlements,
    command = read_layernames(x = layers_ruimtebeslag)
  ),
  tar_target(
    name = lyrs_parcels,
    command = read_layernames(x = layers_perceelgrens)
  ),
  geotargets::tar_terra_vect(
    name = grb_waterways,
    command = get_grb_by_row(
      layer = lyrs_waterways,
      polygons = validation_polygons
    ),
    pattern = cross(lyrs_waterways, validation_polygons),
    filetype = "GPKG"
  ),
  geotargets::tar_terra_vect(
    name = grb_settlements,
    command = get_grb_by_row(
      layer = lyrs_settlements,
      polygons = validation_polygons
    ),
    pattern = cross(lyrs_settlements, validation_polygons),
    filetype = "GPKG"
  ),
  geotargets::tar_terra_vect(
    name = grb_parcels,
    command = get_grb_by_row(
      layer = lyrs_parcels,
      polygons = validation_polygons
    ),
    pattern = cross(lyrs_parcels, validation_polygons),
    filetype = "GPKG"
  ),
  targets::tar_target(
    name = lbg_layers,
    command = get_lbg_layernames(path_to_lbg)
  ),
  geotargets::tar_terra_vect(
    name = lbg_101,
    command = get_lbg(
      path_to_lbg = path_to_lbg,
      layer = lbg_layers,
      from_fields = c("GWSCOD_H", "GWSNAM_H"),
      where_field = "GWSCOD_H",
      where_values = c(1, 2, 11, 12 ,13, 14, 15, 16, 9536),
      flea_value = 101
    ),
    pattern = map(lbg_layers)
  ),
  geotargets::tar_terra_vect(
    name = lbg_104,
    command = get_lbg(
      path_to_lbg = path_to_lbg,
      layer = lbg_layers,
      from_fields = c("GWSCOD_H", "GWSNAM_H"),
      where_field = "GWSCOD_H",
      where_values = c(9),
      flea_value = 104
    ),
    pattern = map(lbg_layers)
  ),
  geotargets::tar_terra_vect(
    name = lbg_101_cropped,
    command = spatvector_crop(x = lbg_101, y = validation_polygons),
    pattern = cross(lbg_101, validation_polygons)
  ),
  geotargets::tar_terra_vect(
    name = lbg_104_cropped,
    command = spatvector_crop(x = lbg_104, y = validation_polygons),
    pattern = cross(lbg_104, validation_polygons)
  ),
  geotargets::tar_terra_vect(
    name = grb_settlements_processed,
    command = process_settlement(
      grb = grb_settlements)
  ),
  geotargets::tar_terra_vect(
    name = grb_water_wtz_processed,
    command = process_water_wtz(
      grb = grb_waterways)
  ),
  geotargets::tar_terra_vect(
    name = grb_parcels_processed,
    command = process_parcels(
      grb = grb_parcels)
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
    name = vp_water_settlements,
    command = combine_water_settlements(
      water = vp_water,
      settlements = grb_settlements_processed,
      polygons = validation_polygons
    ),
    pattern = map(vp_water)
  ),
  geotargets::tar_terra_vect(
    name = vp_water_settlements_cleaned,
    command = postprocess_water_settlements(
      vp_water_settlements
    ),
    pattern = map(vp_water_settlements)
  ),
  geotargets::tar_terra_vect(
    name = vp_water_settlements_singletarget,
    command = single_wsp(vp_water_settlements_cleaned)
  ),
  geotargets::tar_terra_vect(
    name = prelabeled_validation_polygons,
    command = intersect_validation_polygons(
      wsp_target = vp_water_settlements_singletarget,
      lu_changecat = lu_changecats
    ),
    pattern = map(lu_changecats)
  )
  #,
  #geotargets::tar_terra_vect(
  #  name = grb_waterways_processed,
  #  command = process_waterways(grb_waterways)
  #)



  # apply majority filter, use 3 by 3 block


  # merge the samples,check if locations were sampled > 1
  # deduplicate them, keeping all metadata


  # write out sampling polygons

  # any other stuff that could be done automatically
  # (maybe some validation steps)

)
