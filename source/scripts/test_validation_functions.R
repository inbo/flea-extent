# BEGIN code rmd NCA_validatingextent/src/_evaluation.Rmd
# adapted, test only for nara map
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
# validation_data function with data_root = flea_data
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

a <- reschange1$byClass[, c(1, 2, 5, 6)]
rownames(a) <- str_remove(rownames(a), "Class: ")
a |>
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

# END code rmd NCA_validatingextent/src/_evaluation.Rmd

##############################################################################
##############################################################################

# function confusion_matrix() calculates the map-relevant version of error
# matrix
# map-relevant means that each p_ij is weighted by area proportions for each
# stratum (map class)
# so this is different from default caret::confusionMatrix which was used in
# function calculate_accuracy()
# see https://pages.cms.hu-berlin.de/EOL/gcg_eo/06_accuracy_assessment.html

# TODO check results against mapac R package
cm <- confusion_matrix(
  maparea = maparea$area,
  ma = as.data.frame.matrix(reschange1$table)
)
dim(cm)
cm |> round(digits = 4)

# note maparea$area is not yet an area, but the count of pixels
# and each pixel is 0.01 ha
n_h <- observed_changes |>
  count(map)
maparea <- maparea |>
  left_join(n_h, by = join_by(changecat == map)) |>
  mutate(
    ips = n / area
  )
observed_changes <- observed_changes |>
  left_join(maparea, by = join_by(map == changecat))
aa <- mapac::aa_card(
  data = observed_changes[, c("ref", "map")],
  w = maparea$area/sum(maparea$area),
  strata = levels(observed_changes$ref),
  area = sum(maparea$area),
  confusion_matrix = FALSE,
  olofsson = TRUE
)

waldo::compare(
  x = cm,
  y = aa$cmp,
  tolerance = 1e-10
) #OK

oa_df <- calc_oa(
  maparea = maparea$area,
  ma = as.data.frame.matrix(reschange1$table)
)
oa_df
waldo::compare(oa_df$oa_est, aa$accuracy[1])
waldo::compare(sqrt(oa_df$oa_var), aa$accuracy[2])


ua_pa_df <- calc_ua_pa(
  maparea = maparea$area,
  ma = as.data.frame.matrix(reschange1$table)
)

ua_pa_df
waldo::compare(unname(ua_pa_df$ua_est), aa$stats$ua, tolerance = 1e-10)
waldo::compare(unname(sqrt(ua_pa_df$ua_var)), aa$stats$ua_se, tolerance = 1e-3)
waldo::compare(unname(ua_pa_df$pa_est), aa$stats$pa, tolerance = 1e-10)
waldo::compare(unname(sqrt(ua_pa_df$pa_var)), aa$stats$pa_se, tolerance = 1e-3)


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

mapac::aa_class_accuracy_plot(aa)
mapac::aa_confusion_matrix_flextable(
  aa,
  proportion = TRUE,
  diagonal = TRUE,
  format.body = "%.2f",
  format.accuracy = "%.3f",
  rotate.header = TRUE
  )

areas_df <- calc_areas(
  maparea = maparea$area,
  ma = as.data.frame.matrix(reschange1$table),
  pixelsize = 0.01 # each cell is 100 square meters = 0.01ha
)
areas_df
waldo::compare(unname(areas_df$prop_est), aa$stats$p_i, tolerance = 1e-10)
waldo::compare(unname(areas_df$prop_est), aa$area$proportion, tolerance = 1e-10)
waldo::compare(
  unname(areas_df$area_est_ha),
  aa$area$area * 0.01,
  tolerance = 1e-10
)
waldo::compare(
  unname(areas_df$area_rme * areas_df$area_est_ha),
  aa$area$area_ci * 0.01,
  tolerance = 1e-5
)

# estimate areas using ReGenesees
#################################
prop_area_h <- observed_changes |>
  distinct(map, area) |>
  mutate(
    prop_area_h = area / sum(area)
  )
# area Flanders
sum(prop_area_h$area)
# check if total area matches
sum(aa$area$area) - sum(prop_area_h$area)

# create survey data.frame (long format)
# the variable oa is needed for estimation
svydata <- observed_changes |>
  mutate(
    weights = 1/ ips,
    ids = paste0("id_", 1:n()),
    ones = 1,
    oa = as.numeric(map == ref)
  ) |>
  mutate(
    n_ref = n(),
    .by = ref
  ) |>
  inner_join(
    prop_area_h |> select(-area), by = join_by(map)
  )

# create design object (could deal with more complex models)
design <- ReGenesees::e.svydesign(
  data = svydata,
  ids = ~ ids,
  strata = ~ map,
  weights = ~ weights,
  fpc = ~ ips)

# create dataframe containing marginal population totals
df.pop <- ReGenesees::pop.template(
  data = svydata,
  calmodel = ~ map - 1
) # see what the template should look like
ReGenesees::pop.desc(df.pop)
df.pop <- maparea |>
  select(changecat, area) |>
  pivot_wider(
    names_from = changecat,
    names_prefix = "map",
    values_from = area) |>
  as.data.frame()

# calibrate on map marginal totals
# this adds the calibration weights to the design
cal <- ReGenesees::e.calibrate(
  design = design,
  df.population = df.pop,
  calmodel = ~ map - 1
)
summary(cal)
# but because in this case, apparently,
# the ratio between calibrated weights and initial weights equals 1;
# so using cal or design in estimation will have the same result
# this is probably logical because inverse probability weights (N_h/n_h) already
# incorporate the map marginal totals per stratum
# this would not be the case if we calibrate on different strata or
# additional strata and other auxiliary information
# the calibration weights are stored in cal$variables$weights.cal
ReGenesees::g.range(cal)

# check calibration
ReGenesees::check.cal(cal)
summary(weights(cal))
ReGenesees::svystatTM(cal, ~ ones) # area of Flanders in 0.01 ha units

# estimate areas for ref
cal_ref_areas <- ReGenesees::svystatTM(
  design = cal,
  y =  ~ ref,
  estimator = "Total",
  conf.int = TRUE,
  deff = TRUE)
# note that above with y = ~ map would just recover known marginal totals (SE=0)
waldo::compare(aa$area$area, cal_ref_areas$Total, tolerance = 1e-10)
waldo::compare(aa$area$area_ci, 1.96*cal_ref_areas$SE, tolerance = 1e-4)


# can we estimate the areas also via Bayes formula?
# $P(\text{referentie}=A) = \frac{UA}{PA}P(\text{kaart} = A)$
# first a naive attempt for just the estimator without SE
waldo::compare(
  aa$area$area, aa$stats$ua / aa$stats$pa * maparea$area,
  tolerance = 1e-4
) # only difference is one NA resulting from division pa 0
plot(aa$area$area, aa$stats$ua / aa$stats$pa * maparea$area)
abline(0, 1)
# how to estimate $\frac{UA}{PA}P(\text{kaart} = A)$ with ReGenesees?
# oa = n_jj = variable equalling one if ref == map, 0 otherwise
# n = n_i. = n_map the sample size per map stratum
# n_ref = n_.k = the number of ref cases
# UA / PA = (n_jj / n_i. ) / (n_jj / n_.k ) simplifies to n_.k / n_i.
# P(\text{kaart} = i) = area_i
#  n_.k / n_i. * area_i = sum over i of n_ik / n_i. * area_i
# compare with eq 9 Olofsson:
# W_i = area proportion of map class i = area_i / area_tot
# n_ik = cell count map i ref k
# n_i. = row sum map i
# area estimator for ref class k =
# total map area times sum over i (from 1 to 25) of W_i * n_ik/n_i.
# which again equals sum over i (from 1 to 25) of area_i * n_ik/n_i.

# note that area_i/n_i. are inverse probability weights
# because our survey data are in long format, we can
# represent n_ik as a vector of ones in combi with by = ~ref
# to obtain the same as before using either an expression in svystatL or
cal_ref_areas_checkL <- ReGenesees::svystatL(
  design = cal,
  expr =  expression(ones),
  by = ~ ref,
  conf.int = TRUE,
  deff = TRUE)
cal_ref_areas_checkTM <- ReGenesees::svystatTM(
  design = cal,
  y =  ~ ones,
  by = ~ ref,
  conf.int = TRUE,
  deff = TRUE)
waldo::compare(cal_ref_areas_checkL$ones, cal_ref_areas$Total)
waldo::compare(cal_ref_areas_checkTM$Total.ones, cal_ref_areas$Total)
waldo::compare(cal_ref_areas_checkL$SE.ones, cal_ref_areas$SE)
waldo::compare(cal_ref_areas_checkTM$SE.Total.ones, cal_ref_areas$SE)

# can it also be used to calculate OA, UA and PA?
# OA: YES
cal_oa <- ReGenesees::svystatTM(
  design = cal,
  y =  ~ oa,
  estimator = "Mean",
  conf.int = TRUE,
  deff = TRUE)
waldo::compare(aa$accuracy[1], cal_oa$Mean, tolerance = 1e-10)
waldo::compare(aa$accuracy[2], cal_oa$SE, tolerance = 1e-5)

# UA: YES
# (which need be estimated by row of contingency table which means by "map")
cal_ua <- ReGenesees::svystatTM(
  design = cal,
  y =  ~ oa,
  by = ~ map,
  estimator = "Mean",
  conf.int = TRUE,
  deff = TRUE)
waldo::compare(aa$stats$ua, cal_ua$Mean.oa, tolerance = 1e-10)
waldo::compare(aa$stats$ua_se, cal_ua$SE.Mean.oa, tolerance = 1e-2) # more conservative
hist(aa$stats$ua_se -  cal_ua$SE.Mean.oa)
plot(aa$stats$ua_se, cal_ua$SE.Mean.oa)
abline(0, 1)

# PA
cal_pa <- ReGenesees::svystatTM(
  design = cal,
  y =  ~ oa,
  by = ~ ref,
  estimator = "Mean",
  conf.int = TRUE,
  deff = TRUE)
waldo::compare(aa$stats$pa, cal_pa$Mean.oa, tolerance = 1e-10)
waldo::compare(aa$stats$pa_se, cal_pa$SE.Mean.oa, tolerance = 1e-2)
hist(aa$stats$pa_se -  cal_pa$SE.Mean.oa)
plot(aa$stats$pa_se, cal_pa$SE.Mean.oa)
abline(0, 1)

areas_df %>%
  inner_join(
    aa$stats, by = join_by(class)
  ) %>%
  separate(
    class,
    c("class_p1", "class_p2"),
    sep = "-",
    remove = FALSE
  ) %>%
  mutate(change = class_p1 != class_p2) %>%
  mutate(
    class = reorder(
      sprintf(
        "%s\n(n = %s; ua = %s; pa = %s)",
        class, n_points, round(ua, 2), round(pa, 2)
        ),
      area_est_ha
    )
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

# estimate total sample size for stratified random sampling

n_est <- data.frame(oa_ses = seq(0.001, 0.01, 0.001)) |>
  mutate(
    n_tot = mapac::sample_size(
      oa_se = oa_ses,
      w = maparea$area / sum(maparea$area),
      ua = aa$stats$ua
    )
  )
n_est |>
  ggplot() +
  geom_point(
    aes(x = oa_ses, y = n_tot)
  ) +
  scale_y_log10() +
  labs(
    x = "Standard error overall accuracy",
    y = "Total sample size (stratified sampling)"
  )


# Allocation strategy: first allocate between 50 and 100 to the change classes
# allocate remainder proportional to the area of the stable classes

###############################################################################
###############################################################################
# recode into gain, loss, stable presence

binary_change <- function(
    data,
    lg,
    year1 = "lg2013_label",
    year2 = "lg2016_label") {
  binary <- vector("list", length = length(lg))
  binary <- setNames(binary, lg)
  for (i in lg) {
    binary[[i]] <- paste0(
      stringr::str_detect(data[[year1]], i) %>% as.numeric(),
      stringr::str_detect(data[[year2]], i) %>% as.numeric()
    )
  }
  bind_cols(data, binary)
}

categorize_land_use_change <- function(b) {
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
}

lg <- c("Field", "Urban", "High green", "Open nature", "Other")

simplechange <- c("Stable presence", "Stable absence", "Loss", "Gain")

mapdata <- maparea |>
  select(changecat, area, n) |>
  separate(
    changecat,
    c("lg2013_label", "lg2016_label"),
    sep = "-",
    remove = FALSE
  ) |>
  binary_change(lg = lg) %>%
  rowwise() %>%
  mutate(stable = ifelse(
    lg2013_label == lg2016_label,
    "stable", "changed"
  ) %>%
    as.factor()) %>%
  ungroup() %>%
  mutate(
    across(
      all_of(lg),
      \(x) categorize_land_use_change(x),
      .names = "{.col}_changecat"
    ),
    across(
      ends_with("_changecat"),
      \(x) factor(x, levels = simplechange)
    )
  )

# prepare validation data for simplified change cats
validationdata <- observed_changes |>
  as_tibble() |>
  select(map, ref) |>
  separate(
    map,
    c("lg2013_map", "lg2016_map"),
    sep = "-",
    remove = FALSE
  ) |>
  separate(
    ref,
    c("lg2013_ref", "lg2016_ref"),
    sep = "-",
    remove = FALSE
  ) |>
  binary_change(lg = lg, year1 = "lg2013_map", year2 = "lg2016_map") |>
  rename_with(.fn = ~paste0(.x, "_map"), .cols = all_of(lg)) |>
  binary_change(lg = lg, year1 = "lg2013_ref", year2 = "lg2016_ref") |>
  rename_with(.fn = ~paste0(.x, "_ref"), .cols = all_of(lg)) |>
  mutate(
    across(
      all_of(paste0(lg, "_map")),
      \(x) categorize_land_use_change(x),
      .names = "{.col}_changecat"
    ),
    across(
      all_of(paste0(lg, "_ref")),
      \(x) categorize_land_use_change(x),
      .names = "{.col}_changecat"
    ),
    across(
      ends_with("_changecat"),
      \(x) factor(x, levels = simplechange)
    )
  )
# Field strata info for simplified change cats
mapfield <- mapdata |>
  summarise(
    across(
      c(area, n),
      sum
      ),
    .by = Field_changecat
  ) |>
  mutate(
    ips = n / area
)

# calculate accuracies for field
aa_field <- mapac::aa_card(
  data = validationdata[, c("Field_ref_changecat", "Field_map_changecat")] |>
    as.data.frame(),
  w = (mapfield$area/sum(mapfield$area))[order(mapfield$Field_changecat)],
  strata = simplechange,
  area = sum(mapfield$area),
  confusion_matrix = FALSE,
  olofsson = TRUE
)

mapac::aa_confusion_matrix_flextable(aa_field)
mapac::aa_class_accuracy_plot(aa_field)

aa_field2 <- mapac::aa_stratified(
  stratum = validationdata$map,
  reference = validationdata$Field_ref_changecat,
  map = validationdata$Field_map_changecat,
  h = levels(validationdata$map),
  N_h = (maparea$n)[order(maparea$changecat)]
)

mapac::aa_confusion_matrix_flextable(aa_field2)
mapac::aa_class_accuracy_plot(aa_field2)


