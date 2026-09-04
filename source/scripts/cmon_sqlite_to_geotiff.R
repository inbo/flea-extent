library(grtsdb)
library(ggplot2)
library(dplyr)
library(sf)
library(terra)
git_root <- rprojroot::find_root(rprojroot::is_git_root)

flea_data <- gsub(
  pattern = "flea-extent", replacement = "flea-data", x = git_root
)

flea_bbox <- rbind(x = c(20000, 259000), y = c(153000, 250000))
flea_ext <- ext(as.numeric(flea_bbox), xy = TRUE)

con <- connect_db(file.path(flea_data, "data/c-mon/grts.sqlite"))
# use internal code from grtsdb::extract_sample
# to get SQL syntax
samplesize <- 100
bbox <- flea_bbox
cellsize <- 10
offset <- NULL
level <- n_level(bbox = bbox, cellsize = cellsize)
fields <- dbListFields(con, sprintf("level%02i", level))
fields <- fields[grep("^x[[:digit:]]*$", fields)]
center <- rowMeans(bbox)
midpoint <- 2^(level - 1) - 0.5
where <- sprintf(
  "%s %s %f", rep(fields, 2),
  rep(c(">=", "<="), each = length(center)),
  (as.vector(bbox) - center) / cellsize + midpoint
)
where <- paste(where, collapse = " AND ")
fields <- sprintf(
  "(%1$s - %2$f) * %3$f + %4$f AS %1$sc",
  fields, midpoint, cellsize, center
)

# amend the query to extract all cells
# order from topleft to bottomright
sql <- sprintf(
  "SELECT %s, ranking FROM level%02i WHERE %s ORDER BY -x2c, x1c",
  paste(fields, collapse = ", "), level, where
)

allcells <- RSQLite::dbGetQuery(con, sql) # 9Gb

dbDisconnect(con)

class(allcells)
head(allcells)
# S4 method for class 'data.frame'
# If the value is "xyz", the matrix or data.frame x must have at least two
# columns, the first with x (or longitude) and the second with y (or latitude)
# coordinates that represent the centers of raster cells.
# The additional columns are the values associated with the raster cells
rast(
  x = allcells,
  type = "xyz",
  crs = "EPSG:31370",
  digits = 6,
  extent = flea_ext
) |>
  terra::writeRaster(
    filename = file.path(flea_data, "data/c-mon/flea_cmon_level15.tiff"),
    datatype = "INT4U"
  )

fleagrts <- rast(file.path(flea_data, "data/c-mon/flea_cmon_level15.tiff"))
fleagrts
plot(fleagrts)
origin(fleagrts)
ext(fleagrts)

# convert to base 4 fractional representation
# adapted from n2khab:::convert_dec_to_base4frac
convert_dec_to_base4frac <- function(x, levels = 15) {
  sapply(x, function(x) {
    ifelse(
      is.na(x),
      NA,
      as.double(
        ifelse(# nolint start
          x > 0,
          {
            d <- floor(log(x, 4) + 1)
            paste(
              c("0", "1", "2", "3")[
                as.integer(
                  abs(
                    diff(
                      x %% 4^seq(d, 0)
                    )
                  ) %/% 4^seq(d - 1, 0) + 1)
              ],
              collapse = ""
            )
          },# nolint end
          "0")
      ) / 10^levels
    )
  })
}

sprintf("%.15f", convert_dec_to_base4frac(0:10))
# digits are ordered from deepest nested split level (first digit after decimal
# point = level 15) to first split level (level 1 on the right)

testset <- terra::spatSample(fleagrts, 100, xy = TRUE) |>
  mutate(
    base4frac = sprintf("%.15f", convert_dec_to_base4frac(ranking)),
    level1 = stringr::str_extract(base4frac, "\\d$")
  )

testset |>
  st_as_sf(coords = c("x", "y"), crs = 31370) |>
  ggplot() +
  geom_sf(aes(colour = level1))
