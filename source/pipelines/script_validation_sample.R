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
  packages = c("tibble", "geotargets", "assertthat", "terra"),
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
  )

  # apply majority filter, use 3 by 3 block

  # create table containing all occuring transitions for each land-use

  # reclassify transitions into stable, gained, lost, dynamic for each land-use

  # for each land use, draw a spatially balanced ordered sample
  # within each of the temporal classes
  # possibly exclude some strata that are rare such as "urban - lost"
  # sample sizes? Equal? Or less in stable, dynamic more in gained, lost?

  # merge the samples,check if locations were sampled > 1
  # deduplicate them, keeping all metadata

  # convert sampling locations to square polygons of size 9 x 9

  # write out sampling polygons

  # any other stuff that could be done automatically
  # (maybe some validation steps)

)
