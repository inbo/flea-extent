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
  packages = c("tibble", "geotargets", "assertthat", "terra", "dplyr"),
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
path_to_gdb <- "Z:/Projects/PRJ_FLEA/flea_data.gdb"
path_to_lyr <- "Z:/Projects/PRJ_FLEA/reclass_bwk2016.lyr"
path_to_grts <- file.path(flea_data, "data/c-mon/flea_cmon_level15.tiff")


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
  )
  # convert sampling locations to square polygons of size 9 x 9


  # apply majority filter, use 3 by 3 block


  # merge the samples,check if locations were sampled > 1
  # deduplicate them, keeping all metadata


  # write out sampling polygons

  # any other stuff that could be done automatically
  # (maybe some validation steps)

)
