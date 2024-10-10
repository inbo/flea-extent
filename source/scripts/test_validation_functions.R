# code rmd NCA_validatingextent adapted, test only for nara map
git_root <- rprojroot::find_root(rprojroot::is_git_root)
source(file.path(git_root, "source/scripts/nca_functions.R"))
flea_data <- gsub(
  pattern = "flea-extent", replacement = "flea-data", x = git_root
)
library(tidyverse)
library(terra)
library(git2rdata)
library(plotly)
library(sf)
library(INBOtheme)

conflicted::conflicts_prefer(dplyr::filter)

# load validation.Rdata
# this file was produced by
# validation_data(data_root = flea_data)
# but it is currently not working, probably due to different package versions
# anyway, not really needed for this script
load(
  file.path(
    flea_data,
    "data/validation.Rdata"
  )
)

# preprocessing steps
#####################
ls <- list(
  nara2016 = terra::rast(
    file.path(flea_data, "data/2016/LG2016_finaal_ori.tif")
  ),
  nara2013 = terra::rast(
    file.path(flea_data, "data/2013/LG2013_finaal_ori.tif")
  )
)

data <- ls %>% map(function(x) {
  cleanmapdata(
    data = x,
    points_id = points_id,
    tbltrans = tbltrans,
    type =
      ifelse(str_detect(names(x), "ori"), "nara",
        ifelse(str_detect(names(x), "bos"), "bosw",
          ifelse(str_detect(names(x), "BAK"), "bak",
            ifelse(str_detect(names(x), "gras"),
              "gras", "all"
            )
          )
        )
      ),
    year = as.numeric(gsub("[^0-9.-]", "", names(x)))
  )
})

ls <- list(
  nara = file.path(flea_data, "data/2013_2016/LG2013_LG2016_finaal_ori.csv")
)
area <- ls %>% map(function(x) {
  Cleanchangeareadata(
    file = x,
    tbltrans = tbltrans,
    type =
      ifelse(str_detect(x, "ori"), "nara",
        ifelse(str_detect(x, "bos"), "bosw",
          ifelse(str_detect(x, "BAK"), "bak",
            ifelse(str_detect(x, "gras"),
              "gras", "all"
            )
          )
        )
      )
  )
})
rm(ls)

# note area$area is not yet an area, but the count of pixels
# and each pixel is 0.01 ha

ttl <- list(nara = "the original land use map")

this_map <- "nara"

title <- ttl[[this_map]]
mapdata1 <- data[[str_c(this_map, "2013")]]$valid_eng
refdata1 <- as.factor(points_id$lu13oord_eng)
mapdata2 <- data[[str_c(this_map, "2016")]]$valid_eng
refdata2 <- as.factor(points_id$lu16oord_eng)
maparea <- area[[this_map]]
# note maparea$area is not yet an area, but the count of pixels
# and each pixel is 0.01 ha

# calculations
##############

# map of 2013
res1 <- calculate_accuracy(mapdata1, refdata1)
# OA
round(100 * unname(res1$overall["Accuracy"]), digits = 0)
# UA PA
# UA = precision (this does NOT equal specificity) = positive predictive value
# PA = recall = sensitivity = true positive fraction
t(res1$byClass) |> round(digits = 2)

# map of 2016
res2 <- calculate_accuracy(mapdata2, refdata2)
# OA
round(100 * unname(res2$overall["Accuracy"]), digits = 2)
# UA PA
# UA = precision (this does NOT equal specificity) = positive predictive value
# PA = recall = sensitivity = true positive fraction
t(res2$byClass) |> round(digits = 2)

data.frame(
  users = c(res1$byClass[, 5], res2$byClass[, 5]),
  producers = c(res1$byClass[, 6], res2$byClass[, 6]),
  lu = rep(str_remove(rownames(res1$byClass), "Class: "), 2),
  year = c(
    rep(2013, nrow(res1$byClass)),
    rep(2016, nrow(res2$byClass))
  )
) %>%
  mutate(
    lu = as.factor(lu),
    year = as.factor(year)
  ) %>%
  ggplot() +
  geom_point(aes(x = users, y = producers, color = lu, shape = year)) +
  scale_color_discrete(name = "land use") +
  xlab("user's accuracy (precision)") +
  ylab("producer's accuracy (recall)") +
  theme_bw()


observed_changes <- data.frame(
  map = as.factor(str_c(mapdata1, mapdata2, sep = "-")),
  ref = as.factor(str_c(refdata1, refdata2, sep = "-")),
  mapchange = as.factor(ifelse(mapdata1 == mapdata2,
    "No change", "Change"
  )),
  refchange = as.factor(ifelse(refdata1 == refdata2,
    "No change", "Change"
  ))
) %>%
  mutate(
    map = factor(map, levels = sort(unique(c(levels(map), levels(ref))))),
    map = factor(map, levels = levels(map))
  )
observed_changes %>%
  pivot_longer(cols = 1:2, names_to = "type", values_to = "change") %>%
  group_by(change) %>%
  summarize(
    map = sum(type == "map"),
    ref = sum(type == "ref")
  ) %>%
  ungroup()

reschange1 <- calculate_accuracy(observed_changes$map, observed_changes$ref)
reschange2 <- calculate_accuracy(
  observed_changes$mapchange, observed_changes$refchange
)

A <- reschange1$byClass[, c(1, 2, 5, 6)]
rownames(A) <- str_remove(rownames(A), "Class: ")
A |>
  round(digits = 2) |>
  knitr::kable(
    caption = str_c(
      "Accuracy of land use changes  for ",
      title,
      " (only some accuracy measures are shown but others can be requested)."
    ),
    row.names = TRUE,
    booktabs = TRUE
  )

reschange2$table |>
  round(digits = 2) |>
  knitr::kable(
    caption = str_c(
      "Confusion matrix for land use changes,
      only comparing change / no change for ",
      title, ". "
    ),
    row.names = TRUE,
    booktabs = TRUE
  )

b <- reschange2$byClass[c(1, 2, 5, 6)]
b |> round(digits = 2)

ov <- validation_uncertainty(
  ma = as.data.frame.matrix(reschange1$table),
  maparea = maparea$area,
  pixelsize = 0.01
) # each cell is 100 square meters = 0.01ha

ov

plot_validation_data(ov)


##############################################################################
##############################################################################

# this function calculates the map-relevant version of error matrix
# map-relevant means that each p_ij is weighted by area proportions for each
# stratum (map class)
# so this is different from default caret::confusionMatrix
cm <- confusion_matrix(
  maparea = maparea$area,
  ma = as.data.frame.matrix(reschange1$table)
)
dim(cm)
cm |> round(digits = 4)

oa_df <- calc_oa(
  maparea = maparea$area,
  ma = as.data.frame.matrix(reschange1$table)
)
oa_df

ua_pa_df <- calc_ua_pa(
  maparea = maparea$area,
  ma = as.data.frame.matrix(reschange1$table)
)
ua_pa_df

ua_pa_df %>%
  separate(
    class,
    c("class_p1", "class_p2"),
    sep = "-",
    remove = FALSE
  ) %>%
  mutate(change = class_p1 != class_p2) %>%
  ggplot(aes(x = pa_est, y = ua_est, colour = change)) +
  geom_abline(alpha = 0.5) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = class), size = 2) +
  geom_errorbar(aes(ymin = ua_low, ymax = ua_high), alpha = 0.3) +
  geom_errorbarh(aes(xmin = pa_low, xmax = pa_high), alpha = 0.3) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))

areas_df <- calc_areas(
  maparea = maparea$area,
  ma = as.data.frame.matrix(reschange1$table),
  pixelsize = 0.01 # each cell is 100 square meters = 0.01ha
)
areas_df

areas_df %>%
  separate(
    class,
    c("class_p1", "class_p2"),
    sep = "-",
    remove = FALSE
  ) %>%
  mutate(change = class_p1 != class_p2) %>%
  mutate(
    class = reorder(class, area_est_ha)
  ) %>%
  ggplot() +
  geom_pointrange(
    aes(
      x = class,
      y = area_est_ha,
      ymin = area_low_ha,
      ymax = area_high_ha,
      colour = area_low_ha < 0
    )
  ) +
  scale_y_log10() +
  coord_flip() +
  facet_grid(paste0("Change: ", change) ~ ., scales = "free", space = "free")

# relative margins of error larger than 1 result in
# negative lower bound of design-based confidence interval
areas_df %>%
  separate(
    class,
    c("class_p1", "class_p2"),
    sep = "-",
    remove = FALSE
  ) %>%
  mutate(change = class_p1 != class_p2) %>%
  mutate(
    class = reorder(class, area_rme),
    fprop_map_rbias = cut(
      prop_map_rbias,
      breaks = c(
        min(prop_map_rbias) - 0.01,
        -0.1,
        0.1,
        1,
        max(prop_map_rbias) + 0.01
      ),
      labels = c(
        "Underestimation\nmore than 10%",
        "Relative bias\nbetween -10% and 10%",
        "Overestimation\nbetween 10% and 100%",
        "Overestimation\nmore than 100%"
      )
    )
  ) %>%
  ggplot() +
  geom_point(
    aes(
      x = class,
      colour = fprop_map_rbias,
      y = area_rme,
      size = abs(prop_map_rbias)
    )
  ) +
  scale_y_continuous(
    name = "Relative margin of error",
    labels = scales::percent
  ) +
  coord_flip() +
  facet_grid(paste0("Change: ", change) ~ ., scales = "free", space = "free")

# bias - variance tradeoff is always in favor of area estimation via sample
# except for field-field, which is better estimated from pixel counting
# as judged by mean squared error = variance + bias^2
# and assuming variance is zero for the map and bias is zero for the sample
areas_df %>%
  separate(
    class,
    c("class_p1", "class_p2"),
    sep = "-",
    remove = FALSE
  ) %>%
  mutate(change = class_p1 != class_p2) %>%
  mutate(
    class = reorder(class, prop_map_bias)
  ) %>%
  ggplot() +
  geom_point(
    aes(
      x = class,
      y = prop_map_bias,
      size = area_rme,
      colour = prop_mse_map < prop_mse_sample
    ),
    alpha = 0.3
  ) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(
    "Under (-) or over (+) estimation\nPercentage of area of interest",
    labels = scales::percent
  ) +
  #  scale_colour_gradient2(midpoint = 0, mid = "white") +
  coord_flip() +
  facet_grid(paste0("Change: ", change) ~ ., scales = "free", space = "free")
