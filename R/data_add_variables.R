lmessage <- function(...) {}

# lmessage <- function(...)
#   do.call(message, list(...))

add_variables <- function(dcd, dmg, rule_spec = NULL) {
  dat <- left_join(dcd, dmg, by = "champs_deid")

  if (is.null(rule_spec))
    rule_spec <- rules

  rgx <- paste(
    "immediate_cause_of_death_etiol",
    "morbid_condition_.._etiol",
    "underlying_cause_factor_",
    sep = "|"
  )
  nms <- names(dat)
  nms <- nms[grepl(rgx, nms)]
  nms <- nms[!grepl("_othr", nms)]

  rgx2 <- "ic_champs_group_desc|[0-9]_champs_group_desc"
  nms2 <- names(dat)
  nms2 <- nms2[grepl(rgx2, nms2)]

  for (rule in rule_spec) {
    if (rule$method == "etiol") {
      lmessage("Adding variable for ", rule$check,
        " [", rule$var, "]... ", appendLF = FALSE)
      tmp <- dat[, nms] == rule$check
      dat[[rule$var]] <- as.integer(apply(tmp, 1,
        function(x) sum(x, na.rm = TRUE)) > 0)
      lmessage("(", sum(dat[[rule$var]]), ")")
    } else if (rule$method == "champs_group") {
      var2 <- paste0(rule$var, "Sum")
      lmessage("Adding variables for ", rule$check,
        " [", rule$var, ", ", var2, "]... ", appendLF = FALSE)
      tmp <- dat[, nms2] == rule$check
      dat[[rule$var]] <- as.integer(apply(tmp, 1,
        function(x) sum(x, na.rm = TRUE)) > 0)
      dat[[var2]] <- as.integer(dat[[rule$var]]
        | dat$uc_champs_group_desc == rule$check)
      lmessage("(", sum(dat[[rule$var]]), ")")
    }
  }

  lmessage("Adding hosp_los_24h2, hosp_los_48h2, infectious_cause, ",
    "acquired, acquired2, pmirange")
  starts <- c(-Inf, 0, 4, 7, 10, 13, 16, 19, 22, 25, Inf)
  labels <- c("bad", "0 to 3", "4 to 6", "7 to 9", "10 to 12",
    "13 to 15", "16 to 18", "19 to 21", "22 to 24", "Over 24h")
  tmp <- cut(dat$calc_postmortem_hrs, starts - 0.5, labels)
  tmp[tmp == "bad"] <- NA
  tmp <- droplevels(tmp)

  # rng_lvls <- c("0 to 3", "4 to 6", "7 to 9", "10 to 12", "13 to 15",
  #   "16 to 18", "19 to 21", "22 to 24", "Over 24h")
  # pmirange <- rep(rng_lvls, ceiling(nrow(dat) / length(rng_lvls)))[1:nrow(dat)]

  dat %>%
    mutate(
      hosp_los_24h2 = ifelse(is.na(hosp_los_24h), 0, hosp_los_24h),
      hosp_los_48h2 = ifelse(is.na(hosp_los_48h), 0, hosp_los_48h),
      infectious_cause = as.integer(
        DiarrhealDSum | CISum | HIVSum | LRISum | MalariaSum
        | MeaslesSum | MESum | NeoSepsisSum | OtherInfectionSum
        | SepsisSum | SyphilisSum | TBSum | URISum),
      acquired = ifelse(
        case_type_desc == "Death in the first 24 hours" |
        location_of_death == "community" |
        hosp_los_24h2 == 1, "Community", "Facility"),
      acquired2 = ifelse(
        case_type_desc == "Stillbirth" |
        case_type_desc == "Death in the first 24 hours" |
        location_of_death == "community" |
        hosp_los_24h2 == 1 |
        hosp_los_48h2 == 1, "Community", "Facility"),
      pmirange = tmp
    )
}

add_variables_old <- function(dcd_res, dcd_meas, rule_spec = NULL) {
  dat <- dplyr::left_join(dcd_res, dcd_meas, by = "champsid")

  if (is.null(rule_spec))
    rule_spec <- rules

  rgx <- paste(
    "immediate_cause_of_death_etiol",
    "morbid_condition_.._etiol",
    "underlying_cause_factor_",
    sep = "|"
  )
  nms <- names(dat)
  nms <- nms[grepl(rgx, nms)]
  nms <- nms[!grepl("_othr", nms)]

  rgx2 <- "ic_champs_group_desc|[0-9]_champs_group_desc"
  nms2 <- names(dat)
  nms2 <- nms2[grepl(rgx2, nms2)]

  for (rule in rule_spec) {
    if (rule$method == "etiol") {
      lmessage("Adding variable for ", rule$check,
        " [", rule$var, "]... ", appendLF = FALSE)
      tmp <- dat[, nms] == rule$check
      dat[[rule$var]] <- as.integer(apply(tmp, 1,
        function(x) sum(x, na.rm = TRUE)) > 0)
      lmessage("(", sum(dat[[rule$var]]), ")")
    } else if (rule$method == "champs_group") {
      var2 <- paste0(rule$var, "Sum")
      lmessage("Adding variables for ", rule$check,
        " [", rule$var, ", ", var2, "]... ", appendLF = FALSE)
      tmp <- dat[, nms2] == rule$check
      dat[[rule$var]] <- as.integer(apply(tmp, 1,
        function(x) sum(x, na.rm = TRUE)) > 0)
      dat[[var2]] <- as.integer(dat[[rule$var]]
        | dat$uc_champs_group_desc == rule$check)
      lmessage("(", sum(dat[[rule$var]]), ")")
    }
  }

  lmessage("Adding hosp_los_24h2, hosp_los_48h2, infectious_cause, ",
    "acquired, acquired2, pmirange")
  starts <- c(-Inf, 0, 4, 7, 10, 13, 16, 19, 22, 25, Inf)
  labels <- c("bad", "0 to 3", "4 to 6", "7 to 9", "10 to 12",
    "13 to 15", "16 to 18", "19 to 21", "22 to 24", "Over 24h")
  tmp <- cut(dat$calc_postmortem_hrs, starts - 0.5, labels)
  tmp[tmp == "bad"] <- NA
  tmp <- droplevels(tmp)

  dat %>%
    mutate(
      hosp_los_24h2 = ifelse(is.na(hosp_los_24h), 0, hosp_los_24h),
      hosp_los_48h2 = ifelse(is.na(hosp_los_48h), 0, hosp_los_48h),
      infectious_cause = as.integer(
        DiarrhealDSum | CISum | HIVSum | LRISum | MalariaSum
        | MeaslesSum | MESum | NeoSepsisSum | OtherInfectionSum
        | SepsisSum | SyphilisSum | TBSum | URISum),
      acquired = ifelse(
        case_type == "CH01404" |
        calc_location == "community" |
        hosp_los_24h2 == 1, "Community", "Facility"),
      acquired2 = ifelse(
        case_type == "CH01404" |
        calc_location == "community" |
        hosp_los_24h2 == 1 |
        hosp_los_48h2 == 1, "Community", "Facility"),
      pmirange = pmirange
    )
}
