# copy of NCA_functions.R in project NCA_validatingextent
# nolint start
library(assertthat)
library(caret) # Confusion matrix maken


cleanmapdata <- function(data = data, points_id, tbltrans, type, year) {
  maps <- terra::extract(
    x = data,
    y = terra::vect(
      st_as_sf(
        x = points_id[, c("POINT_X", "POINT_Y")],
        coords = c("POINT_X", "POINT_Y"),
        crs = "EPSG:31370"
      )
    )
  ) # extract de waardes van de referentiedata
  names(maps)[2] <- "landgebruik"
  maps <- maps %>%
    mutate(
      x = points_id[, "POINT_X"],
      y = points_id[, "POINT_X"],
      code = as.factor(landgebruik),
      landgebruik = recode_factor(code,
        "1" = "Open natuur", "2" = "Bos",
        "3" = "Grasland", "4" = "Akker",
        "5" = "Urbaan", "6" = "Laag groen",
        "7" = "Hoog groen", "8" = "Water",
        "9" = "Overig"
      ),
      type = type,
      year = year
    ) %>%
    left_join(tbltrans[, -1], by = c("code" = "lucode")) %>%
    droplevels()
  rm(data)
  return(maps)
}

Cleanchangeareadata <- function(file, tbltrans, type) {
  maparea <- read_csv2(file = file)
  maparea %>%
    mutate(
      LG2013 = as.factor(LG2013),
      LG2016 = as.factor(LG2016),
      type = type
    ) %>%
    left_join(tbltrans[, c("lucode", "valid_eng")],
      by = c("LG2013" = "lucode")
    ) %>%
    dplyr::select(-LG2013) %>%
    rename(LG2013 = valid_eng) %>%
    left_join(tbltrans[, c("lucode", "valid_eng")],
      by = c("LG2016" = "lucode")
    ) %>%
    dplyr::select(-LG2016) %>%
    rename(LG2016 = valid_eng) %>%
    filter(LG2013 != "Water" &
      LG2016 != "Water") %>%
    mutate(
      changecat = as.factor(str_c(LG2013, LG2016, sep = "-")),
      changebool = as.factor(ifelse(LG2013 == LG2016,
        "No change", "Change"
      ))
    ) %>%
    group_by(changecat, changebool, type) %>%
    summarize(area = sum(Count)) %>%
    ungroup() %>%
    arrange(changecat)
}
# mapdata <- nara13$valid
# refdata <- as.factor(points_id$lu13oord)
# both arrays need to be factor variables with the same levels.
calculate_accuracy <- function(mapdata, refdata) {
  assert_that(length(mapdata) == length(refdata),
    msg = "Length of the arrays is not equal"
  )
  assert_that(
    nlevels(mapdata) == nlevels(refdata) &
      all(levels(mapdata) %in% levels(refdata)),
    msg = "The data are not factors or don't have the same levels."
  )
  conf <- confusionMatrix(
    data = mapdata,
    reference = refdata
  )
  # overall accuracy
  # sum(diag(conf$table))/sum(conf$table)
  return(conf)
}

confusion_matrix <- function(maparea, ma) {
  aoi <- sum(maparea) # calculate the area proportions for each map class
  propmaparea <- maparea / aoi

  # convert the absolute cross tab into a probability cross tab
  ni <- rowSums(ma) # number of reference points per map class
  propma <- as.matrix(ma / ni * propmaparea)
  propma[is.nan(propma)] <- 0 # for classes with ni. = 0
  return(propma)
}



calc_oa <- function(maparea, ma, propma = NULL) {
  if (is.null(propma)) {
    propma <- confusion_matrix(maparea = maparea, ma = ma)
  }
  # overall accuracy (Eq. 1 in Olofsson et al. 2014)
  oa <- sum(diag(propma))

  # variance of overall accuracy (Eq. 5 in Olofsson et al. 2014)
  ni <- rowSums(ma) # number of reference points per map class
  aoi <- sum(maparea) # calculate the area proportions for each map class
  propmaparea <- maparea / aoi
  ua <- diag(propma) / rowSums(propma)
  v_oa <- sum(propmaparea^2 * ua * (1 - ua) / (ni - 1), na.rm = TRUE)
  me_oa <- 1.96 * sqrt(v_oa)
  oa_low <- oa - me_oa
  oa_high <- oa + me_oa
  return(
    data.frame(oa_est = oa, oa_var = v_oa, oa_low = oa_low, oa_high = oa_high)
  )
}

calc_ua_pa <- function(maparea, ma, propma = NULL) {
  if (is.null(propma)) {
    propma <- confusion_matrix(maparea = maparea, ma = ma)
  }
  dyn <- rownames(ma)
  aoi <- sum(maparea) # calculate the area proportions for each map class
  propmaparea <- maparea / aoi
  ni <- rowSums(ma) # number of reference points per map class

  # estimate the accuracies
  # user's accuracy (Eq. 2 in Olofsson et al. 2014)
  # producer's accuracy (Eq. 3 in Olofsson et al. 2014)
  pa <- diag(propma) / colSums(propma)
  ua <- diag(propma) / rowSums(propma)

  # estimate confidence intervals for the accuracies
  # variance of user's accuracy (Eq. 6 in Olofsson et al. 2014)
  # variance of producer's accuracy (Eq. 7 in Olofsson et al. 2014)
  v_ua <- ua * (1 - ua) / (ni - 1)
  n_j <- vector(mode = "numeric", length = length(dyn))
  aftersumsign <- vector(mode = "numeric", length = length(dyn))
  for (cj in seq_len(length(dyn))) {
    n_j[cj] <- sum(maparea / ni * ma[, cj], na.rm = TRUE)
    aftersumsign[cj] <- sum(maparea[-cj]^2 * ma[-cj, cj] / ni[-cj] *
      (1 - ma[-cj, cj] / ni[-cj]) /
      (ni[-cj] - 1), na.rm = TRUE)
  }
  v_pa <- 1 / n_j^2 * (maparea^2 * (1 - pa)^2 * ua * (1 - ua) / (ni - 1) +
    pa^2 * aftersumsign)
  v_pa[is.nan(v_pa)] <- 0

  ua_me <- 1.96 * sqrt(v_ua)
  pa_me <- 1.96 * sqrt(v_pa)
  ua_low <- ua - ua_me
  ua_high <- ua + ua_me
  pa_low <- pa - pa_me
  pa_high <- pa + pa_me


  return(tibble::tibble(
    class = dyn,
    ua_est = ua,
    ua_var = v_ua,
    ua_low = ua_low,
    ua_high = ua_high,
    pa_est = pa,
    pa_var = v_pa,
    pa_low = pa_low,
    pa_high = pa_high
  ))
}

calc_areas <- function(maparea, ma, pixelsize = 0.01, propma = NULL) {
  if (is.null(propma)) {
    propma <- confusion_matrix(maparea = maparea, ma = ma)
  }
  # Estimate area
  # proportional area estimation
  # proportion of area (Eq. 8 and 9 in Olofsson et al. 2014)
  propareaest <- colSums(propma)

  # standard errors of the area estimation (Eq. 10 in Olofsson et al. 2014)
  dyn <- rownames(ma)
  ni <- rowSums(ma)
  aoi <- sum(maparea) # calculate the area proportions for each map class
  propmaparea <- maparea / aoi
  v_propareaest <- vector(mode = "numeric", length = length(dyn))
  for (cj in seq_len(length(dyn))) {
    v_propareaest[cj] <- sum((propmaparea * propma[, cj] - propma[, cj]^2) /
      (ni + 0.001 - 1))
    # + 0.001 voor klassen met maar 1 punt
  }
  v_propareaest[is.na(v_propareaest)] <- 0
  me_propareaest <- 1.96 * sqrt(v_propareaest)

  out <- tibble::tibble(
    class = dyn,
    n_points = ni,
    prop_est = propareaest,
    prop_var = v_propareaest
  ) |>
    mutate(
      prop_low = prop_est - me_propareaest,
      prop_high = prop_est + me_propareaest,
      prop_map_unadjusted = propmaparea,
      prop_map_bias = prop_map_unadjusted - prop_est,
      prop_map_rbias = prop_map_bias / prop_est,
      area_pixelcount_ha = maparea * pixelsize, # in ha
      area_est_ha = prop_est * aoi * pixelsize, # in ha
      area_low_ha = area_est_ha - me_propareaest * aoi * pixelsize, # in ha
      area_high_ha = area_est_ha + me_propareaest * aoi * pixelsize, # in ha
      area_bias_ha = area_pixelcount_ha - area_est_ha,
      area_rbias_ha = area_bias_ha / area_est_ha,
      area_rme = me_propareaest * aoi * pixelsize / area_est_ha,
      prop_mse_map = prop_map_bias^2,
      prop_mse_sample = v_propareaest
    )
  return(out)
}


# maparea is the surface area of each change class
# ma is the n_{ij} confusion matrix for the change classes
validation_uncertainty <- function(ma, maparea, pixelsize) {
  dyn <- rownames(ma)
  aoi <- sum(maparea) # calculate the area proportions for each map class
  propmaparea <- maparea / aoi
  ni <- rowSums(ma) # number of reference points per map class
  propma <- confusion_matrix(maparea = maparea, ma = ma)

  pa <- diag(propma) / colSums(propma)
  # estimate the accuracies
  oa <- sum(diag(propma))
  # overall accuracy (Eq. 1 in Olofsson et al. 2014)
  ua <- diag(propma) / rowSums(propma)
  # user's accuracy (Eq. 2 in Olofsson et al. 2014)
  # producer's accuracy (Eq. 3 in Olofsson et al. 2014)

  # estimate confidence intervals for the accuracies
  v_oa <- sum(propmaparea^2 * ua * (1 - ua) / (ni - 1), na.rm = TRUE)
  # variance of overall accuracy (Eq. 5 in Olofsson et al. 2014)

  v_ua <- ua * (1 - ua) / (rowSums(ma) - 1)
  # variance of user's accuracy (Eq. 6 in Olofsson et al. 2014)

  # variance of producer's accuracy (Eq. 7 in Olofsson et al. 2014)
  n_j <- array(0, dim = length(dyn))
  aftersumsign <- array(0, dim = length(dyn))
  for (cj in seq_len(length(dyn))) {
    n_j[cj] <- sum(maparea / ni * ma[, cj], na.rm = TRUE)
    aftersumsign[cj] <- sum(maparea[-cj]^2 * ma[-cj, cj] / ni[-cj] *
      (1 - ma[-cj, cj] / ni[-cj]) /
      (ni[-cj] - 1), na.rm = TRUE)
  }
  v_pa <- 1 / n_j^2 * (maparea^2 * (1 - pa)^2 * ua * (1 - ua) / (ni - 1) +
    pa^2 * aftersumsign)
  v_pa[is.nan(v_pa)] <- 0

  ### Estimate area

  # proportional area estimation
  propareaest <- colSums(propma)
  # proportion of area (Eq. 8 in Olofsson et al. 2014)

  # standard errors of the area estimation (Eq. 10 in Olofsson et al. 2014)
  v_propareaest <- array(0, dim = length(dyn))
  for (cj in seq_len(length(dyn))) {
    v_propareaest[cj] <- sum((propmaparea * propma[, cj] - propma[, cj]^2) /
      (rowSums(ma) + 0.001 - 1)) # + 0.001 voor klassen met maar 1 punt
  }
  v_propareaest[is.na(v_propareaest)] <- 0

  # produce the overview table
  ov <- as.data.frame(round(propma, 3))
  ov$class <- rownames(ov)
  ov <- dplyr::select(ov, class)
  ov$totpunt <- rowSums(ma)
  ov$area_ha <- round(maparea * pixelsize) # in ha
  ov$prop_area <- round(propmaparea, 3)
  ov$adj_proparea <- round(propareaest, 3)
  ov$ci_adj_proparea <- round(1.96 * sqrt(v_propareaest), 3)
  ov$adj_area <- round(ov$adj_proparea * aoi * pixelsize, 3)
  # in ha
  ov$ci_adj_area <- round(1.96 * sqrt(v_propareaest) * aoi * pixelsize, 3)
  # in ha
  ov$ua <- round(ua, 3)
  ov$ci_ua <- round(1.96 * sqrt(v_ua), 3)
  ov$pa <- round(pa, 3)
  ov$ci_pa <- round(1.96 * sqrt(v_pa), 3)
  rownames(ov) <- colnames(ma)
  ov$oa <- c(round(oa, 3), rep(NA, times = length(dyn) - 1))
  ov$ci_oa <- c(round(1.96 * sqrt(v_oa), 3), rep(NA, times = length(dyn) - 1))
  ov
}


plot_validation_data <- function(ov) {
  plot_val <- ov %>%
    dplyr::select(class, area_ha, adj_area, ci_adj_area, ua, pa) %>%
    mutate(conf.low = adj_area - ci_adj_area, conf.high = adj_area +
      ci_adj_area) %>%
    mutate(signif0 = ifelse(conf.low <= 0, "", "*")) %>%
    separate(class, c("lu13", "lu16"), sep = "-", remove = FALSE) %>%
    unite("classfull", lu13:lu16, sep = " > ", remove = FALSE) %>%
    mutate(paperc = scales::percent(pa, accuracy = 1)) %>%
    mutate(uaperc = scales::percent(ua, accuracy = 1))

  options(scipen = 999)

  bar <- ggplot(
    plot_val %>% filter(lu13 != lu16),
    aes(x = classfull, y = adj_area, text = paste(
      "PA:", paperc, " - UA:", uaperc,
      "\nValidated area:", round(adj_area), " ha",
      "\nArea on the map:", round(area_ha), " ha",
      "\nCI:", round(conf.low), " ha - ", round(conf.high), " ha"
    ))
  ) +
    geom_bar(aes(fill = lu13),
      stat = "identity", position = "dodge",
      width = 0.7
    ) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
      width = 0.2,
      colour = "black", position = position_dodge(width = 0.7)
    ) +
    geom_point(aes(y = area_ha), colour = "black") +
    labs(y = "Area (ha)", fill = "Class 2013") +
    geom_text(aes(y = -1000, label = signif0), colour = inbo_hoofd) +
    theme(
      axis.title.x = element_text(margin = margin(
        t = 5, r = 0, b = 0,
        l = 0
      ), hjust = 0),
      axis.line.x = element_line(color = "black"),
      axis.title.y = element_blank(),
      axis.line.y = element_line(color = "black"),
      panel.grid.major.x = element_line(colour = "grey", linetype = "dotted"),
      panel.grid.major.y = element_blank(),
      legend.key.size = unit(0.3, "cm")
    ) +
    coord_flip()

  return(bar)
}

validation_data <- function(data_root) {
  punten <- read_vc(
    "validatiepunten",
    root = file.path(data_root, "data")
  ) # Gevalideerde punten
  combine <- read_vc(
    "combine",
    root = file.path(data_root, "data")
  )
  # Combine van de validatieklassenkaart van 2013
  # en 2016 -> geeft de oppervlakte van de landgebruiksveranderingen en
  # van de stabiele klassen
  # Omzettingstabel landgebruiken
  lu <- c(
    "Open natuur", "Bos", "Grasland", "Akker", "Urbaan", "Laag groen",
    "Hoog groen", "Water", "Overig"
  )
  lucode <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
  valid <- c(
    "Open natuur", "Hoog groen", "Open natuur", "Akker", "Urbaan",
    "Open natuur", "Hoog groen", "Water", "Overig"
  )
  valid_eng <- c(
    "Open nature", "High green", "Open nature", "Field", "Urban",
    "Open nature", "High green", "Water", "Other"
  )
  validcode <- c(1, 2, 1, 4, 5, 1, 2, 8, 9)
  tbltrans2 <- data.frame(lu, lucode, valid, valid_eng, validcode)
  tbltrans <- tbltrans2 %>% mutate(
    valid = as.factor(valid),
    valid_eng = as.factor(valid_eng),
    lucode = as.factor(lucode),
    validcode = as.factor(validcode),
    lu = as.factor(lu)
  )

  combine <- combine %>%
    mutate(
      LG2013_ChangeCla = as.factor(LG2013_ChangeCla),
      LG2016_ChangeCla = as.factor(LG2016_ChangeCla)
    ) %>%
    left_join(tbltrans[, c("lucode", "valid_eng")],
      by = c("LG2013_ChangeCla" = "lucode")
    ) %>%
    rename(LG2013 = valid_eng) %>%
    left_join(tbltrans[, c("lucode", "valid_eng")],
      by = c("LG2016_ChangeCla" = "lucode")
    ) %>%
    rename(LG2016 = valid_eng) %>%
    filter(LG2013 != "Water" &
      LG2016 != "Water") %>%
    mutate(
      changecat = as.factor(str_c(LG2013, LG2016, sep = "-")),
      changebool = as.factor(ifelse(LG2013 == LG2016,
        "No change", "Change"
      ))
    ) %>%
    group_by(changecat, changebool) %>%
    summarize(area = sum(Count)) %>%
    ungroup()
  # OPMAAK DATASET -----------------------------------------------------
  # De dataset met validatiepunten bevat informatie over de verandering van een
  # cel (0/1) en de aard van de verandering (klasse A -> klasse B). Voor de
  # validatie werden achteraf een aantal moeilijk te onderscheiden klassen
  # samengenomen. De analyse kan dus op 3 niveaus uitgevoerd worden:
  #   (1) verandering - geen verandering, (2) originele landgebruiksklassen en
  # (3) de geaggregeerde klassen. In het onderstaande script worden deze
  # drie validatiesets apart aangemaakt.

  # 1. Puntenset -----
  points <- punten %>%
    mutate(
      across(where(is.factor), \(x) as.character(x) |> na_if("<Null>"))
    ) %>%
    # Alle <Null> omzetten in NA
    gather(key = klasse, value = oordeel, X2013:change.9) %>%
    # Van breed formaat naar lang formaat
    separate(col = klasse, into = c("klasse", "eval"), sep = "\\.") %>%
    replace_na(list(eval = "0")) %>%
    # Voeg evaluator toe
    drop_na() %>%
    # Alle NA laten vallen
    rename(
      x = POINT_X, y = POINT_Y, lu2013 = change2013, lu2016 = change2016,
      verandering = type
    ) %>%
    mutate(oordeel = recode(oordeel, "Nee" = "nochange", "Ja" = "change")) %>%
    mutate(klasse = gsub("\\..*", "", klasse)) %>%
    # Alles na "." weglaten -> \\.. definieert . en * betekent "alles na"
    mutate(klasse = recode(klasse,
      X2013 = "lu2013", X2016 = "lu2016",
      change = "verandering"
    )) %>%
    mutate(lu2013 = recode(lu2013,
      "1" = "Open natuur", "2" = "Bos", "3" = "Grasland", "4" = "Akker",
      "5" = "Urbaan",
      "6" = "Laag groen", "7" = "Hoog groen", "8" = "Water", "9" = "Overig"
    )) %>%
    # codes naar tekst
    mutate(lu2016 = recode(lu2016,
      "1" = "Open natuur", "2" = "Bos", "3" = "Grasland",
      "4" = "Akker", "5" = "Urbaan", "6" = "Laag groen",
      "7" = "Hoog groen", "8" = "Water", "9" = "Overig"
    )) %>%
    # codes naar tekst
    rowwise() %>%
    mutate(lu_c = ifelse(klasse == "lu2013" & lu2013 == oordeel,
      1, ifelse(klasse == "lu2016" & lu2016 == oordeel,
        1, 0
      )
    )) %>%
    # Check of de gevalideerde landgebruiken overeenkomen met de LG van de kaart
    mutate(change_c = ifelse(verandering == oordeel,
      1, 0
    )) %>%
    # Check of de beoordeling "change/nochange" overeenkomt met die van de
    # LG-kaart
    ungroup() %>%
    group_by(objectid) %>%
    mutate(n = n() / 3) %>%
    ungroup() %>%
    filter(n >= 1) %>%
    # alleen punten die gevalideerd zijn
    arrange(objectid) %>%
    left_join(dplyr::select(tbltrans2, valid, lu), by = c("lu2013" = "lu")) %>%
    rename(luval13 = valid) %>%
    # validatieklassen toevoegen (= aggregatie van oorspronkelijke lu-klassen)
    left_join(dplyr::select(tbltrans2, valid, lu), by = c("lu2016" = "lu")) %>%
    rename(luval16 = valid) %>%
    left_join(dplyr::select(tbltrans2, valid, lu), by = c("oordeel" = "lu")) %>%
    rename(oordeelval = valid) %>%
    group_by(objectid, eval) %>%
    mutate(oordeelval = ifelse(is.na(oordeelval),
      ifelse(identical(oordeelval[1], oordeelval[2]),
        "nochange", "change"
      ), oordeelval
    )) %>%
    # Aanpassen beoordeling "verandering" ->
    # als de validatieklasse 2 x hetzelfde
    # is per evaluator, dan "nochange"
    rowwise() %>%
    mutate(
      veranderingval =
        ifelse(luval13 == luval16, "nochange", "change")
    ) %>%
    mutate(
      luval_c =
        ifelse(klasse == "lu2013" & luval13 == oordeelval, 1,
          ifelse(klasse == "lu2016" & luval16 == oordeelval, 1, 0)
        )
    ) %>%
    # Check of de gevalideerde landgebruiken overeenkomen met de LG van de kaart
    mutate(changeval_c = ifelse(veranderingval == oordeelval, 1, 0)) %>%
    # Check of de beoordeling "change/nochange" overeenkomt met die van
    # de LG-kaart
    ungroup()

  points_id <- points %>%
    # puntenset met 1 waarde per objectid
    # filter(klasse != "verandering") %>%
    spread(klasse, oordeelval) %>%
    group_by(objectid, eval) %>%
    summarise(
      luval13 = first(luval13), luval16 = first(luval16),
      verand = first(veranderingval),
      lu13oord = first(na.omit(lu2013)), lu16oord = first(na.omit(lu2016)),
      verandoord = first(na.omit(verandering))
    ) %>%
    group_by(objectid) %>%
    sample_n(1) %>% # 1 random classificatie w gekozen bij conflict tss experten
    na.omit() %>%
    mutate(codeval13 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", luval13,
      perl = TRUE
    )) %>%
    mutate(codeval16 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", luval16,
      perl = TRUE
    )) %>%
    unite(changeclass, codeval13, codeval16, sep = "_") %>%
    mutate(codeval13 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", lu13oord,
      perl = TRUE
    )) %>%
    mutate(codeval16 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", lu16oord,
      perl = TRUE
    )) %>%
    unite(changeclassref, codeval13, codeval16, sep = "_") %>%
    as.data.frame() %>%
    mutate(
      changeclass = as.factor(changeclass),
      changeclassref = as.factor(changeclassref)
    ) %>%
    filter( # No water classes
      !changeclass %in% c("W_W", "ON_W", "O_W"),
      !changeclassref %in% c("W_W", "ON_W", "O_W")
    ) %>%
    droplevels() %>%
    left_join(punten[, c("objectid", "POINT_X", "POINT_Y")],
      by = c("objectid" = "objectid")
    ) %>%
    # left_join(unique(tbltrans[,c("valid", "valid_eng")]),
    #           by = c("lu13oord" = "valid")) %>%
    # rename(valid_eng = lu13oord) %>%
    left_join(unique(tbltrans2[, c("valid", "valid_eng")]),
      by = c("lu13oord" = "valid")
    ) %>%
    rename(lu13oord_eng = valid_eng) %>%
    left_join(unique(tbltrans2[, c("valid", "valid_eng")]),
      by = c("lu16oord" = "valid")
    ) %>%
    rename(lu16oord_eng = valid_eng)
  save(points_id, tbltrans, combine, file = "data/validation.Rdata")

  ############################# Get areas #################################
  lgarea <- combine %>%
    dplyr::select(LG2013_ChangeCla, LG2016_ChangeCla, Count) %>%
    rename(
      lu2013 = LG2013_ChangeCla, lu2016 = LG2016_ChangeCla,
      count = Count
    ) %>%
    mutate(
      lu2013 = as.factor(lu2013),
      lu2016 = as.factor(lu2016)
    ) %>%
    left_join(dplyr::select(tbltrans, valid, lucode),
      by = c("lu2013" = "lucode")
    ) %>%
    rename(val2013 = valid) %>%
    left_join(dplyr::select(tbltrans, valid, lucode),
      by = c("lu2016" = "lucode")
    ) %>%
    rename(val2016 = valid) %>%
    mutate(codeval13 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", val2013,
      perl = TRUE
    )) %>%
    mutate(codeval16 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", val2016,
      perl = TRUE
    )) %>%
    unite(class, codeval13, codeval16, sep = "_") %>%
    dplyr::select(class, count) %>%
    group_by(class) %>%
    summarise(count = sum(count)) %>%
    filter(!str_detect(class, "W")) %>%
    # alles met water in weglaten
    droplevels() %>%
    mutate(area = count / sum(count))

  # change-no change
  lgarea_change <- lgarea %>%
    mutate(change = ifelse(!class %in% c("A_A", "HG_HG", "O_O", "ON_ON", "U_U"),
      "change", "nochange"
    )) %>%
    group_by(change) %>%
    summarise(count = sum(count), area = sum(area))

  # 2.3. Originele landgebruiksklassen ----

  lgarea_lu <- combine %>%
    dplyr::select(LG2013_ChangeCla, LG2016_ChangeCla, Count) %>%
    rename(
      lu2013code = LG2013_ChangeCla, lu2016code = LG2016_ChangeCla,
      count = Count
    ) %>%
    mutate(
      lu2013code = as.factor(lu2013code),
      lu2016code = as.factor(lu2016code)
    ) %>%
    left_join(dplyr::select(tbltrans, lucode, lu),
      by = c("lu2013code" = "lucode")
    ) %>%
    rename(lu2013 = lu) %>%
    left_join(dplyr::select(tbltrans, lucode, lu),
      by = c("lu2016code" = "lucode")
    ) %>%
    rename(lu2016 = lu) %>%
    mutate(code13 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", lu2013,
      perl = TRUE
    )) %>%
    mutate(code16 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", lu2016,
      perl = TRUE
    )) %>%
    unite(class, code13, code16, sep = "_") %>%
    dplyr::select(class, count) %>%
    group_by(class) %>%
    summarise(count = sum(count)) %>%
    filter(!str_detect(class, "W")) %>%
    # alles met water in weglaten
    droplevels() %>%
    mutate(area = count / sum(count))
}
# nolint end
