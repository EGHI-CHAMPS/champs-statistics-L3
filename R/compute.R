# TODO: if specimen_types is NULL, don't do any filtering
# TODO: add annotations to the data that provide information that will go into
# a caption
# TODO: allow multiple condition and pathogen specifications, using AND logic

#' Tabulate the number of cases of a given condition found in the causal chain
#'   by site and by case type, within the context of all cases
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param condition A string specifying the condition to tabulate causal chain
#'   presence for. To list all possibilities, see
#'   \code{\link{valid_conditions}}.
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @export
#' @examples
#' \dontrun{
#' d <- load_data(...)
#' calc_cc_allcases_by_site_age(d, condition = "GBS")
#' }
calc_cc_allcases_by_site_age <- function(d, condition, sites = NULL) {
  check_champs_data(d)
  condition <- check_valid_condition(condition)

  ridx <- sapply(rules, function(x) x$var)
  if (!condition %in% ridx)
    stop("'", condition, "' is not a valid condition.\n",
      "Call valid_conditions() to get a list of valid conditions.")

  dat <- d$dcd_join

  if (!is.null(sites)) {
    sites <- check_valid_sites(d, sites)
    dat <- dplyr::filter(dat, site %in% sites)
    exclude <- dplyr::setdiff(levels(dat$site), sites)
    dat$site <- forcats::fct_drop(dat$site, only = exclude)
  }

  denominator <- stats::xtabs(~ case_type_desc + site,
    data = dat) %>% stats::addmargins()
  numerator <- stats::xtabs(~ case_type_desc + site,
    data = dplyr::filter(dat, !!rlang::sym(condition) == 1)) %>%
    stats::addmargins()

  structure(list(
    condition = condition,
    numerator = numerator,
    denominator = denominator
  ), class = c("champs_computed", "cc_allcases_by_site_age"))
}

#' Tabulate the number of cases of a given condition found in the causal chain
#'   by site and by case type, within the context of cases where the pathogen(s)
#'   are detected in TAC results (dying from vs. dying with)
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param condition A string specifying the condition to tabulate causal chain
#'   presence for. To list all possibilities, see
#'   \code{\link{valid_conditions}}.
#' @param pathogen A string specifying the pathogen to count positive cases of
#'   in the TAC results.
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @param specimen_types An optional vector of specimen types to include in the
#'   calculation. If not provided, all specimen types will be used. See
#'   \code{\link{valid_specimen_types}}.
#' @export
#' @examples
#' \dontrun{
#' d <- load_data(...)
#' specimen_types <- c(
#'   "Cerebrospinal fluid sample",
#'   "Tissue specimen from lung",
#'   "Whole blood",
#'   "Rectal swab",
#'   "Plasma or spun blood specimen"
#' )
#'
#' calc_cc_detected_by_site_age(d,
#'   condition = "GBS",
#'   pathogen = "Group B Streptococcus",
#'   specimen_types = specimen_types)
#' }
calc_cc_detected_by_site_age <- function(d, condition, pathogen,
  sites = NULL, specimen_types = NULL) {

  check_champs_data(d)
  condition <- check_valid_condition(condition)
  pathogen <- check_valid_pathogen(d, pathogen)
  specimen_types <- check_valid_specimen_types(d, specimen_types)

  dd <- d$tac_pivot %>%
    dplyr::left_join(d$dmg, by = "champs_deid") %>%
    dplyr::filter(
      pathogen == !!pathogen
      & result == "Positive"
      & specimen_type %in% !!specimen_types) %>%
    dplyr::select(champs_deid, case_type_desc) %>%
    dplyr::distinct() %>%
    left_join(dplyr::select(d$dmg, champs_deid, site), by = "champs_deid")

  if (!is.null(sites)) {
    sites <- check_valid_sites(d, sites)
    dd <- dplyr::filter(dd, site %in% sites)
    exclude <- dplyr::setdiff(levels(dd$site), sites)
    dd$site <- forcats::fct_drop(dd$site, only = exclude)
  }

  nn <- d$dcd_join %>%
    dplyr::filter(!!rlang::sym(condition) == 1
     & champs_deid %in% dd$champs_deid)

  if (!is.null(sites)) {
    sites <- check_valid_sites(d, sites)
    nn <- dplyr::filter(nn, site %in% sites)
    exclude <- dplyr::setdiff(levels(nn$site), sites)
    nn$site <- forcats::fct_drop(nn$site, only = exclude)
  }

  denominator <- stats::xtabs(~ case_type_desc + site,
    data = dd) %>%
    stats::addmargins()
  numerator <- stats::xtabs(~ case_type_desc + site,
    data = nn) %>%
    stats::addmargins()

  structure(list(
    condition = condition,
    pathogen = pathogen,
    sites = sites,
    specimen_types = specimen_types,
    numerator = numerator,
    denominator = denominator
  ), class = c("champs_computed", "cc_detected_site_case"))
}

#' Tabulate the number of cases where the pathogen(s) are detected in TAC
#'   results in the context of all cases by site and by case type
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param condition A string specifying the condition to tabulate causal chain
#'   presence for. To list all possibilities, see \code{\link{valid_conditions}}.
#' @param pathogen A string specifying the pathogen to count positive cases of
#'   in the TAC results.
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @param specimen_types An optional vector of specimen types to include in the
#'   calculation. If not provided, all specimen types will be used. See
#'   \code{\link{valid_specimen_types}}.
#' @export
#' @examples
#' \dontrun{
#' d <- load_data(...)
#' specimen_types <- c(
#'   "Cerebrospinal fluid sample",
#'   "Tissue specimen from lung",
#'   "Whole blood",
#'   "Rectal swab",
#'   "Plasma or spun blood specimen"
#' )
#'
#' calc_detected_allcases_by_site_age(d,
#'   condition = "GBS",
#'   pathogen = "Group B Streptococcus",
#'   specimen_types = specimen_types)
#' }
calc_detected_allcases_by_site_age <- function(d,
  condition, pathogen, sites = NULL, specimen_types = NULL) {

  check_champs_data(d)
  condition <- check_valid_condition(condition)
  pathogen <- check_valid_pathogen(d, pathogen)
  specimen_types <- check_valid_specimen_types(d, specimen_types)
  # tac_variable <- check_valid_tac_variable(d, tac_variable)

  nn <- d$tac_pivot %>%
    dplyr::left_join(d$dmg, by = "champs_deid") %>%
    dplyr::filter(
      pathogen == !!pathogen
      & result == "Positive"
      & champs_deid %in% d$dcd_join$champs_deid
      & specimen_type %in% !!specimen_types) %>%
    dplyr::select(champs_deid, case_type_desc, site) %>%
    dplyr::distinct()

  if (!is.null(sites)) {
    sites <- check_valid_sites(d, sites)
    nn <- dplyr::filter(nn, site %in% sites)
    exclude <- dplyr::setdiff(levels(nn$site), sites)
    nn$site <- forcats::fct_drop(nn$site, only = exclude)
  }

  dd <- d$dmg

  # dd <- d$tac %>%
  #   select(champs_deid) %>%
  #   distinct() %>%
  #   filter(champs_deid %in% d$dmg$champs_deid) %>%
  #   left_join(d$dmg)
    # dplyr::filter(!!rlang::sym(tac_variable) %in%
    #   c("Negative", "Positive")) %>%

  if (!is.null(sites)) {
    sites <- check_valid_sites(d, sites)
    dd <- dplyr::filter(dd, site %in% sites)
    exclude <- dplyr::setdiff(levels(dd$site), sites)
    dd$site <- forcats::fct_drop(dd$site, only = exclude)
  }

  numerator <- stats::xtabs(~ case_type_desc + site, data = nn) %>%
    stats::addmargins()
  denominator <- stats::xtabs(~ case_type_desc + site, data = dd) %>%
    stats::addmargins()

  structure(list(
    condition = condition,
    pathogen = pathogen,
    sites = sites,
    specimen_types = specimen_types,
    numerator = numerator,
    denominator = denominator
  ), class = c("champs_computed", "detected_allcases_by_site_age"))
}

#' Tabulate the number of cases where the pathogen(s) are detected in TAC
#'   results by DeCoDe result (in causal chain or not, etc.) and by either
#'   site or case type
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param by One of either "site" or "casetype", indicating the second dimension
#'   of tabulation.
#' @param condition A string specifying the condition to tabulate causal chain
#'   presence for. To list all possibilities, see
#'   \code{\link{valid_conditions}}.
#' @param pathogen A string specifying the pathogen to count positive cases of
#'   in the TAC results.
#' @param icds A vector of ICD10 codes to check for the DeCoDe result of
#'   "Contributing (P2)".
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @param specimen_types An optional vector of specimen types to include in the
#'   calculation. If not provided, all specimen types will be used. See
#'  \code{\link{valid_specimen_types}}.
#' @export
#' @examples
#' \dontrun{
#' d <- load_data(...)
#' specimen_types <- c(
#'   "Cerebrospinal fluid sample",
#'   "Tissue specimen from lung",
#'   "Whole blood",
#'   "Rectal swab",
#'   "Plasma or spun blood specimen"
#' )
#'
#' calc_detected_by_decode(d,
#'   by = "site",
#'   condition = "GBS",
#'   pathogen = "Group B Streptococcus",
#'   icds = c("P36.0", "A40.1", "P23.3", "G00.2"),
#'   specimen_types = specimen_types)
#' }
calc_detected_by_decode <- function(d, by = "site", condition, pathogen,
  icds, sites = NULL, specimen_types = NULL) {

  check_champs_data(d)

  if (!by %in% c("site", "casetype"))
    stop("Argument 'by' must be one of 'site' or 'casetype'.")

  vars <- c(
    "champs_deid", condition, "infectious_cause", "site", "case_type_desc",
    paste0("other_significant_condition_0", 1:8)
  )

  dd <- d$tac_pivot %>%
    dplyr::filter(
      pathogen == !!pathogen
      & result == "Positive"
      & champs_deid %in% d$dcd_join$champs_deid
      & specimen_type %in% !!specimen_types) %>%
    dplyr::select(champs_deid) %>%
    dplyr::distinct() %>%
    dplyr::left_join(dplyr::select(d$dcd_join, tidyselect::one_of(vars)),
      by = "champs_deid")

  if (!is.null(sites)) {
    sites <- check_valid_sites(d, sites)
    dd <- dplyr::filter(dd, site %in% sites)
    exclude <- dplyr::setdiff(levels(dd$site), sites)
    dd$site <- forcats::fct_drop(dd$site, only = exclude)
  }

  tmp <- dplyr::select(dd, tidyselect::starts_with("other_significant"))
  dd$p2 <- apply(tmp, 1,
    function(x) sum(x %in% icds, na.rm = TRUE)) > 0

  dd$result <- ""
  dd$result[dd[[condition]] == 1] <- "In causal chain"
  dd$result[dd$result != "In causal chain" & dd$p2 == 1] <- "Contributing (P2)"
  dd$result[dd$result == "" & dd$infectious_cause == 1] <- "Other Infectious"
  dd$result[dd$result == ""] <- "Not in CC/no ID"
  dd$result <- forcats::fct_relevel(dd$result,
    c("In causal chain", "Contributing (P2)",
      "Other Infectious", "Not in CC/no ID"))

  if (by == "site") {
    numerator <- stats::xtabs(~ result + site, data = dd) %>%
      stats::addmargins()
  } else {
    numerator <- stats::xtabs(~ case_type_desc + result, data = dd) %>%
      stats::addmargins()
  }

  structure(list(
    condition = condition,
    pathogen = pathogen,
    icds = icds,
    sites = sites,
    specimen_types = specimen_types,
    numerator = numerator
  ), class = c("champs_computed", "detected_by_decode"))
}

#' Tabulate the number of positive specimens for each case where the pathogen(s)
#'   are detected in TAC results by postmortem interval range (PMI -
#'   time from death to MITS)
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param pathogen A string specifying the pathogen to count positive cases of
#'   in the TAC results.
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @export
#' @examples
#' \dontrun{
#' d <- load_data(...)
#' calc_nspecimen_by_pmi(d, pathogen = "Group B Streptococcus")
#' }
calc_nspecimen_by_pmi <- function(d, pathogen, sites = NULL) {

  message("L3 data currently does not contain postmortem interval...",
    " Using fake data for now to ensure function works with new L3 format.")

  check_champs_data(d)

  dd <- d$tac_pivot %>%
    dplyr::filter(
      pathogen == !!pathogen
      # & result == "Positive"
      & champs_deid %in% d$dcd_join$champs_deid) %>%
    dplyr::group_by(champs_deid) %>%
    dplyr::tally() %>%
    dplyr::right_join(dplyr::select(d$dcd_join, champs_deid, pmirange, site),
      by = "champs_deid") %>%
    dplyr::filter(!is.na(pmirange))

  if (!is.null(sites)) {
    sites <- check_valid_sites(d, sites)
    dd <- dplyr::filter(dd, site %in% sites)
    exclude <- dplyr::setdiff(levels(dd$site), sites)
    dd$site <- forcats::fct_drop(dd$site, only = exclude)
  }

  dd$n[is.na(dd$n)] <- 0
  dd$n <- factor(dd$n)

  numerator <- stats::xtabs(~ n + pmirange, data = dd) %>%
    stats::addmargins()

  structure(list(
    pathogen = pathogen,
    sites = sites,
    numerator = numerator
  ), class = c("champs_computed", "nspecimen_by_pmi"))
}

#' Calculate average postmortem interval (PMI - average time from death to MITS
#'   in hours) by specimen type of NP Only vs. Blood/CSF/Lung and by site
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param pathogen A string specifying the pathogen to count positive cases of
#'   in the TAC results.
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @export
#' @examples
#' \dontrun{
#' d <- load_data(...)
#' calc_pmi_by_specimen_site(d, pathogen = "Group B Streptococcus")
#' }
calc_pmi_by_specimen_site <- function(d, pathogen, sites = NULL) {

  message("L3 data currently does not contain postmortem interval...",
    " Using fake data for now to ensure function works with new L3 format.")

  check_champs_data(d)

  dd <- dplyr::filter(d$tac_pivot,
      pathogen == !!pathogen
      & result == "Positive"
      & champs_deid %in% d$dcd_join$champs_deid) %>%
    dplyr::group_by(champs_deid) %>%
    dplyr::summarise(
      has_np = any(specimen_type %in% "Nasopharyngeal and Oropharyngeal swab"),
      has_other = any(specimen_type %in% c(
        "Cerebrospinal fluid sample",
        "Tissue specimen from lung",
        "Whole blood",
        "Rectal swab",
        "Plasma or spun blood specimen"
      ))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      specimen_type2 = ifelse(has_np & !has_other,
        "NP Only", "Blood, CSF, or Lung")
    ) %>%
    dplyr::left_join(dplyr::select(d$dmg, champs_deid, calc_postmortem_hrs,
      site),
      by = "champs_deid") %>%
    dplyr::filter(!is.na(calc_postmortem_hrs) & calc_postmortem_hrs >= 0)

  if (!is.null(sites)) {
    sites <- check_valid_sites(d, sites)
    dd <- dplyr::filter(dd, site %in% sites)
    exclude <- dplyr::setdiff(levels(dd$site), sites)
    dd$site <- forcats::fct_drop(dd$site, only = exclude)
  }

  dd$specimen_type2 <- forcats::fct_relevel(dd$specimen_type2,
    c("NP Only", "Blood, CSF, or Lung"))

  numerator <- stats::xtabs(calc_postmortem_hrs ~ specimen_type2 +
    site, data = dd) %>%
    stats::addmargins()
  denominator <- stats::xtabs(~ specimen_type2 + site, data = dd) %>%
    stats::addmargins()

  structure(list(
    pathogen = pathogen,
    numerator = numerator,
    denominator = denominator
  ), class = c("champs_computed", "pmi_by_specimen_site"))
}

# internal
get_pmi_decode_data <- function(
  d, condition, pathogen, icds, specimen_types, sites
) {

  message("L3 data currently does not contain postmortem interval...",
    " Using fake data for now to ensure function works with new L3 format.")

  vars <- c(
    "champs_deid", condition, "infectious_cause", "calc_postmortem_hrs",
    paste0("other_significant_condition_0", 1:8), "site", "case_type_desc"
  )

  dd <- d$tac_pivot %>%
    dplyr::filter(
      pathogen == !!pathogen
      & result == "Positive"
      & champs_deid %in% d$dcd_join$champs_deid
      & specimen_type %in% !!specimen_types
    ) %>%
    dplyr::select(champs_deid) %>%
    dplyr::distinct() %>%
    dplyr::left_join(dplyr::select(d$dcd_join, tidyselect::one_of(vars)),
      by = "champs_deid") %>%
    dplyr::filter(!is.na(calc_postmortem_hrs) & calc_postmortem_hrs >= 0)
    # dplyr::filter(calc_postmortem_hrs >= 0 | is.na(calc_postmortem_hrs))

  if (!is.null(sites)) {
    sites <- check_valid_sites(d, sites)
    dd <- dplyr::filter(dd, site %in% sites)
    exclude <- dplyr::setdiff(levels(dd$site), sites)
    dd$site <- forcats::fct_drop(dd$site, only = exclude)
  }

  tmp <- dplyr::select(dd, tidyselect::starts_with("other_significant"))
  dd$p2 <- apply(tmp, 1,
    function(x) sum(x %in% icds, na.rm = TRUE)) > 0

  dd$result <- ""
  dd$result[dd[[condition]] == 1] <- "In causal chain"
  dd$result[dd$result != "In causal chain" & dd$p2 == 1] <- "Contributing (P2)"
  dd$result[dd$result == "" & dd$infectious_cause == 1] <- "Other Infectious"
  dd$result[dd$result == ""] <- "Not in CC/no ID"
  dd$result <- forcats::fct_relevel(dd$result,
    c("In causal chain", "Contributing (P2)",
      "Other Infectious", "Not in CC/no ID"))

  dd
}

#' Calculate average postmortem interval (PMI - average time from death to MITS
#'   in hours) by DeCoDe result and site
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param condition A string specifying the condition to tabulate causal chain
#'   presence for. To list all possibilities, see
#'   \code{\link{valid_conditions}}.
#' @param pathogen A string specifying the pathogen to count positive cases of
#'   in the TAC results.
#' @param icds A vector of ICD10 codes to check for the DeCoDe result of
#'   "Contributing (P2)".
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @param specimen_types An optional vector of specimen types to include in the
#'   calculation. If not provided, all specimen types will be used. See
#'   \code{\link{valid_specimen_types}}.
#' @export
#' @examples
#' \dontrun{
#' d <- load_data(...)
#' specimen_types <- c(
#'   "Cerebrospinal fluid sample",
#'   "Tissue specimen from lung",
#'   "Whole blood",
#'   "Rectal swab",
#'   "Plasma or spun blood specimen"
#' )
#'
#' gbs11_2 <- calc_pmi_by_decode_site(d,
#'   condition = "GBS",
#'   pathogen = "Group B Streptococcus",
#'   icds = c("P36.0", "A40.1", "P23.3", "G00.2"),
#'   specimen_types = specimen_types)
#' }
calc_pmi_by_decode_site <- function(d, condition, pathogen, icds,
  specimen_types = NULL, sites = NULL) {

  check_champs_data(d)

  dd <- get_pmi_decode_data(d, condition, pathogen, icds, specimen_types, sites)

  t1 <- stats::xtabs(calc_postmortem_hrs ~ result + site,
    data = dd) %>%
    stats::addmargins()
  t2 <- stats::xtabs(~ result + site, data = dd) %>%
    stats::addmargins()

  structure(list(
    condition = condition,
    pathogen = pathogen,
    icds = icds,
    specimen_types = specimen_types,
    sites = sites,
    numerator = t1,
    denominator = t2
  ), class = c("champs_computed", "pmi_by_decode_site"))
}

#' Calculate average postmortem interval (PMI - average time from death to MITS
#'   in hours) by case type and DeCoDe result
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param condition A string specifying the condition to tabulate causal chain
#'   presence for. To list all possibilities, see
#'   \code{\link{valid_conditions}}.
#' @param pathogen A string specifying the pathogen to count positive cases of
#'   in the TAC results.
#' @param icds A vector of ICD10 codes to check for the DeCoDe result of
#'   "Contributing (P2)".
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @param specimen_types An optional vector of specimen types to include in the
#'   calculation. If not provided, all specimen types will be used. See
#'   \code{\link{valid_specimen_types}}.
#' @export
#' @examples
#' \dontrun{
#' d <- load_data(...)
#' specimen_types <- c(
#'   "Cerebrospinal fluid sample",
#'   "Tissue specimen from lung",
#'   "Whole blood",
#'   "Rectal swab",
#'   "Plasma or spun blood specimen"
#' )
#' calc_pmi_by_age_decode(d,
#'   condition = "GBS",
#'   pathogen = "Group B Streptococcus",
#'   icds = c("P36.0", "A40.1", "P23.3", "G00.2"),
#'   specimen_types = specimen_types)
#' }
calc_pmi_by_age_decode <- function(d, condition, pathogen, icds,
  specimen_types = NULL, sites = NULL) {

  check_champs_data(d)

  dd <- get_pmi_decode_data(d, condition, pathogen, icds, specimen_types, sites)

  # dd$case_type_desc <- fct_other(dd$casetype_ofcl,
  #   drop = c(
  #     "Death in the first 24 hours",
  #     "Early Neonate (1 to 6 days)",
  #     "Late Neonate (7 to 27 days)"
  #   ),
  #   other_level = "Neonates"
  # )

  # dd$case_type_desc <- forcats::fct_relevel(dd$case_type_desc,
  #   c(
  #     "Stillbirth",
  #     "Neonates",
  #     "Infant (28 days to less than 12 months)",
  #     "Child (12 months to less than 60 Months)"
  #   )
  # )

  t1 <- stats::xtabs(calc_postmortem_hrs ~ case_type_desc + result,
    data = dd) %>%
    stats::addmargins()
  t2 <- stats::xtabs(~ case_type_desc + result, data = dd) %>%
    stats::addmargins()

  structure(list(
    condition = condition,
    pathogen = pathogen,
    icds = icds,
    specimen_types = specimen_types,
    sites = sites,
    numerator = t1,
    denominator = t2
  ), class = c("champs_computed", "pmi_by_age_decode"))
}

# take cases filtered positive TAC for GBS but GBS not equal 1
# (not in causal chain) and infectious cause = 1 and then see
# what pathogens they are

#' Tabulate top other pathogens found in the causal chain for cases that have
#'   positive TAC results for a given pathogen
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param condition A string specifying the condition to tabulate causal chain
#'   presence for. To list all possibilities, see
#'   \code{\link{valid_conditions}}.
#' @param pathogen A string specifying the pathogen to count positive cases of
#'   in the TAC results.
#' @param icds A vector of ICD10 codes to check for the DeCoDe result of
#'   "Contributing (P2)".
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @param specimen_types An optional vector of specimen types to include in the
#'   calculation. If not provided, all specimen types will be used. See
#'   \code{\link{valid_specimen_types}}.
#' @export
#' @examples
#' \dontrun{
#' specimen_types <- c(
#'   "Cerebrospinal fluid sample",
#'   "Tissue specimen from lung",
#'   "Whole blood",
#'   "Rectal swab",
#'   "Plasma or spun blood specimen"
#' )
#' d <- load_data(...)
#' calc_top_tac_pathogens(d,
#'   condition = "GBS",
#'   pathogen = "Group B Streptococcus",
#'   icds = c("P36.0", "A40.1", "P23.3", "G00.2"),
#'   specimen_types = specimen_types)
#' }
calc_top_tac_pathogens <- function(d, condition, pathogen, icds,
  specimen_types = NULL, sites = NULL) {

  check_champs_data(d)

  dd <- d$tac_pivot %>%
    dplyr::filter(
      pathogen == !!pathogen
      & result == "Positive"
      & champs_deid %in% d$dmg$champs_deid
      & specimen_type %in% !!specimen_types) %>%
    dplyr::select(champs_deid) %>%
    dplyr::distinct() %>%
    dplyr::left_join(d$dcd_join, by = "champs_deid") %>%
    dplyr::filter(!!rlang::sym(condition) != 1 & infectious_cause == 1)

  if (!is.null(sites)) {
    sites <- check_valid_sites(d, sites)
    dd <- dplyr::filter(dd, site %in% sites)
    exclude <- dplyr::setdiff(levels(dd$site), sites)
    dd$site <- forcats::fct_drop(dd$site, only = exclude)
  }

  tmp <- dplyr::select(dd, tidyselect::starts_with("other_significant"))
  dd$p2 <- apply(tmp, 1,
    function(x) sum(x %in% icds, na.rm = TRUE)) > 0
  dd <- dplyr::filter(dd, p2 == 0)

  rgx <- paste(
    "champs_deid",
    "immediate_cause_of_death_etiol",
    "morbid_condition_.._etiol",
    "underlying_cause_factor_",
    "ic_champs_group_desc",
    "[0-9]_champs_group_desc",
    "uc_champs_group_desc",
    sep = "|"
  )

  res <- dd %>%
    dplyr::select(tidyselect::matches(rgx)) %>%
    tidyr::pivot_longer(-champs_deid) %>%
    dplyr::filter(value %in% tac_table$assay) %>%
    dplyr::mutate(
      value = recode(value,
        `Candida glabrata` = "Candida spp",
        `Candida albicans` = "Candida spp",
        `Candida glabrata` = "Candida spp",
        `Candida albicans` = "Candida spp",
        `Candida` = "Candida spp",
        `Candida Krusei` = "Candida spp",
        `candida sp` = "Candida spp",
        `Candida spp.` = "Candida spp",
        `Candida tropicalis` = "Candida spp")
    ) %>%
    dplyr::select(-name) %>%
    dplyr::distinct() %>%
    dplyr::group_by(value) %>%
    dplyr::tally() %>%
    dplyr::arrange(-n)

  structure(list(
    condition = condition,
    pathogen = pathogen,
    icds = icds,
    specimen_types = specimen_types,
    sites = sites,
    df = res,
    n = nrow(dd)
  ), class = c("champs_computed", "top_pathogens"))
}

#' Tabulate top pathogens associated with a condition by site of acquisition
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param condition A string specifying the condition to tabulate causal chain
#'   presence for. To list all possibilities, see
#'   \code{\link{valid_conditions}}.
#' @param age_groups Optional vector of age groups.
#'   See \code{\link{valid_age_groups}}.
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @export
#' @examples
#' \dontrun{
#' calc_top_dcd_pathogens_by_acq(d,
#'   condition = "Lower respiratory infections",
#'   age_groups = c(
#'     "Infant (28 days to less than 12 months)",
#'     "Child (12 months to less than 60 Months)")
#' )
#' }
calc_top_dcd_pathogens_by_acq <- function(
  d, condition, age_groups = NULL, sites = NULL
) {
  check_champs_data(d)

  age_groups <- check_valid_age_groups(d, age_groups)

  dd <- d$dcd_pivot
  condition <- check_valid_condition(condition, form = "long")

  if (!is.null(sites)) {
    sites <- check_valid_sites(d, sites)
    dd <- dplyr::filter(dd, site %in% sites)
    exclude <- dplyr::setdiff(levels(dd$site), sites)
    dd$site <- forcats::fct_drop(dd$site, only = exclude)
  }

  dd %>%
    filter(champs_group_desc == condition &
      case_type_desc %in% age_groups) %>%
    group_by(etiol, acquired2) %>%
    summarise(n = length(unique(champs_deid))) %>%
    tidyr::pivot_wider(names_from = acquired2, values_from = n) %>%
    mutate(Total = Community + Facility) %>%
    arrange(-Total) %>%
    select(etiol, Community, Facility, Total) %>%
    filter(etiol != "Other Etiology/Agent")
}

#' Tabulate TAC top pathogens detected and causal chain
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param condition A string specifying the condition to tabulate causal chain
#'   presence for. To list all possibilities, see
#'   \code{\link{valid_conditions}}.
#' @param age_groups Optional vector of age groups.
#'   See \code{\link{valid_age_groups}}.
#' @param specimen_types An optional vector of specimen types to include as columns in the table.
#' @param specimen_abbrv An optional vector of abbreviations to use for the specified \code{specimen_types}.
#'   calculation. If not provided, all specimen types will be used. See
#'   \code{\link{valid_specimen_types}}.
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @export
#' @examples
#' \dontrun{
#' calc_top_tac_pathogens_cc(d,
#'   condition = "Lower respiratory infections",
#'   age_groups = c(
#'     "Infant (28 days to less than 12 months)",
#'     "Child (12 months to less than 60 Months)"),
#'   specimen_types = c(
#'     "Nasopharyngeal and Oropharyngeal swab",
#'     "Tissue specimen from lung"),
#'   specimen_abbrv = c("# NP+", "# Lung+")
#' )
#' }
calc_top_tac_pathogens_cc <- function(
  d, condition, specimen_types = NULL, specimen_abbrv = NULL,
  age_groups = NULL, sites = NULL
) {
  check_champs_data(d)
  condition <- check_valid_condition(condition, form = "long")

  specimen_types <- check_valid_specimen_types(d, specimen_types)
  age_groups <- check_valid_age_groups(d, age_groups)
  # TODO: handle age groups/subgroups appropriately

  if (is.null(specimen_abbrv))
    specimen_abbrv <- specimen_types

  nn <- left_join(d$tac, d$dmg, by = "champs_deid") %>%
    filter(case_type_desc %in% age_groups) %>%
    dplyr::pull(champs_deid) %>%
    unique() %>%
    length()

  nn2 <- d$dcd_pivot %>%
    filter(
      case_type_desc %in% age_groups &
      champs_group_desc == condition) %>%
    dplyr::pull(champs_deid) %>%
    unique() %>%
    length()

  specimen_abbrv <- paste0(specimen_abbrv, " (n=", nn, ")")

  rcd <- specimen_abbrv
  if (is.null(rcd))
    rcd <- specimen_types
  names(rcd) <- specimen_types

  # TODO: add filter by sites...
  # if (!is.null(sites)) {
  #   sites <- check_valid_sites(d, sites)
  #   dd <- dplyr::filter(dd, site %in% sites)
  #   exclude <- dplyr::setdiff(levels(dd$site), sites)
  #   dd$site <- forcats::fct_drop(dd$site, only = exclude)
  # }

  res1 <- d$tac_pivot %>%
    left_join(d$dmg, by = "champs_deid") %>%
    filter(
      case_type_desc %in% age_groups &
      !is.na(pathogen) &
      result == "Positive" &
      specimen_type %in% specimen_types) %>%
    group_by(specimen_type, pathogen) %>%
    summarise(n = length(unique(champs_deid))) %>%
    arrange(-n) %>%
    dplyr::ungroup() %>%
    mutate(specimen_type = recode(specimen_type, !!!rcd)) %>%
    tidyr::pivot_wider(names_from = specimen_type, values_from = n)

  nm <- paste0("# in CC (n=", nn, ")")
  res2 <- d$dcd_pivot %>%
    filter(
      case_type_desc %in% age_groups &
      etiol != "Other Etiology/Agent") %>%
    group_by(etiol) %>%
    summarise(n = length(unique(champs_deid))) %>%
    arrange(-n) %>%
    rename(pathogen = "etiol", !!nm := "n")

  nm <- paste0("# in CC as LRI (n=", nn2, ")")
  res3 <- d$dcd_pivot %>%
    filter(
      case_type_desc %in% age_groups &
      champs_group_desc == "Lower respiratory infections" &
      etiol != "Other Etiology/Agent") %>%
    group_by(etiol) %>%
    summarise(n = length(unique(champs_deid))) %>%
    arrange(-n) %>%
    rename(pathogen = "etiol", !!nm := "n")

  res <- res1 %>%
    left_join(res2, by = "pathogen") %>%
    left_join(res3, by = "pathogen")

  structure(list(
    table = res,
    n = nn,
    n_condition = nn2
  ), class = c("champs_computed", "calc_top_tac_pathogens_cc"))
}

#' Tabulate top etiologies by specified age groups
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param age_groups Optional vector of age groups.
#'   See \code{\link{valid_age_groups}}.
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @export
#' @examples
#' \dontrun{
#' calc_top_pathogens_by_age(d,
#'   age_groups = c(
#'     "Death in the first 24 hours",
#'     "Early Neonate (24-72 hours)",
#'     "Early Neonate (72+hrs to 6 days)",
#'     "Late Neonate (7 to 27 days)")
#' )
#' }
calc_top_etiol_by_age <- function(
  d, age_groups = NULL, sites = NULL
) {
  check_champs_data(d)

  age_groups <- check_valid_age_groups(d, age_groups)
  age_groups_var <- attr(age_groups, "group_var")

  # TODO: support site

  res <- d$dcd_pivot %>%
    dplyr::mutate(
      etiol = recode(etiol,
        `Candida glabrata` = "Candida spp",
        `Candida albicans` = "Candida spp",
        `Candida glabrata` = "Candida spp",
        `Candida albicans` = "Candida spp",
        `Candida` = "Candida spp",
        `Candida Krusei` = "Candida spp",
        `candida sp` = "Candida spp",
        `Candida spp.` = "Candida spp",
        `Candida tropicalis` = "Candida spp")
    ) %>%
    filter(!!rlang::sym(age_groups_var) %in% age_groups &
      infectious_cause == 1) %>%
    group_by(!!rlang::sym(age_groups_var), etiol) %>%
    summarise(n = length(unique(champs_deid))) %>%
    # filter(etiol %in% incl) %>%
    filter(etiol != "Other Etiology/Agent") %>%
    arrange(!!rlang::sym(age_groups_var), etiol)
    # mutate(
    #   etiol = recode(etiol, !!!rcd),
    #   etiol = factor(etiol, levels = incl_names))

  top <- res %>%
    group_by(etiol) %>%
    summarise(n = sum(n)) %>%
    filter(etiol != "Other Etiology/Agent") %>%
    arrange(-n)

  # denominators
  denoms <- d$dcd_pivot %>%
    filter(!!rlang::sym(age_groups_var) %in% age_groups) %>%
    group_by(!!rlang::sym(age_groups_var)) %>%
    summarise(n = length(unique(champs_deid)))

  # no etiology
  tmp <- d$dcd_pivot %>%
    filter(!!rlang::sym(age_groups_var) %in% age_groups &
      infectious_cause == 1 &
      is.na(etiol)) %>%
    group_by(!!rlang::sym(age_groups_var), champs_deid) %>%
    tally()
  noet <- tmp %>%
    filter(n == max(tmp$n)) %>%
    group_by(!!rlang::sym(age_groups_var)) %>%
    tally()

  structure(list(
    etiol_counts = res,
    top = top,
    denominators = denoms,
    no_etiol = noet
  ), class = c("champs_computed", "calc_top_etiol_by_age"))
}

# if (is.null(incl_names))
#   incl_names <- incl
# rcd <- incl_names
# if (is.null(rcd))
#   rcd <- incl
# names(rcd) <- incl

# incl <- c("Acinetobacter baumannii",
#   "Staphylococcus aureus",
#   "Candida spp",
#   "Klebsiella pneumoniae",
#   "Escherichia coli",
#   "Streptococcus agalactiae")

# incl_names <- c("A. baumanii", "S. aureus", "Candida spp",
#   "K. pneumoniae", "E. coli", "GBS")

#' Tabulate cases with condition in causal chain by age and site of acquisition
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param condition A string specifying the condition to tabulate causal chain
#'   presence for. To list all possibilities, see
#'   \code{\link{valid_conditions}}.
#' @param age_groups Optional vector of age groups.
#'   See \code{\link{valid_age_groups}}.
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @export
#' @examples
#' \dontrun{
#' calc_cc_allcases_by_age_acq(d,
#'   condition = "Klebsiella pneumoniae",
#'   age_groups = valid_age_subgroups(d)
#' )
#' }
calc_cc_allcases_by_age_acq <- function(
  d, condition, age_groups = NULL, sites = NULL
) {
  check_champs_data(d)
  condition <- check_valid_condition(condition, form = "long")

  age_groups <- check_valid_age_groups(d, age_groups)
  age_groups_var <- attr(age_groups, "group_var")

  res <- d$dcd_pivot %>%
    filter(etiol == condition) %>%
    group_by(!!rlang::sym(age_groups_var), acquired2) %>%
    summarise(n = length(unique(champs_deid))) %>%
    tidyr::pivot_wider(names_from = acquired2, values_from = n,
      values_fill = list(n = 0)) %>%
    mutate(Total = Community + Facility) %>%
    select(!!rlang::sym(age_groups_var), Community, Facility, Total) %>%
    arrange(!!rlang::sym(age_groups_var))

  # denominator
  denoms <- d$dcd_pivot %>%
    group_by(!!rlang::sym(age_groups_var)) %>%
    summarise(n = length(unique(champs_deid))) %>%
    arrange(!!rlang::sym(age_groups_var))

  structure(list(
    table = res,
    denominators = denoms
  ), class = c("champs_computed", "calc_cc_allcases_by_age_acq"))
}

#' Tabulate syndrome combinations for a specified condition
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param condition A string specifying the condition to tabulate causal chain
#'   presence for. To list all possibilities, see
#'   \code{\link{valid_conditions}}.
#' @param syndrome_names A vector specifying conditions to tabulate
#'   combinations of presence of for the given \code{condition}
#' @param syndrome_values An optional vector of the same length of
#'   \code{syndrome_names} that specifies the mapping from the official
#'   condition names to the values to use for reporting. See example.
#' @param specimen_types An optional vector of specimen types to include in the
#'   calculation. If not provided, all specimen types will be used. See
#'   \code{\link{valid_specimen_types}}.
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @export
#' @examples
#' \dontrun{
#' calc_syndrome_combinations(d,
#'   condition = "Streptococcus pneumoniae",
#'   conditions = c(
#'     "Lower respiratory infections",
#'     "Meningitis/Encephalitis",
#'     "Neonatal sepsis",
#'     "Congenital infection"),
#'   syndrome_values = c(
#'     "Pneumonia",
#'     "Meningitis",
#'     "Sepsis",
#'     "Sepsis"),
#'   specimen_types = c(
#'     "Cerebrospinal fluid sample",
#'     "Tissue specimen from lung",
#'     "Whole blood")
#' )
#' }
calc_syndrome_combinations <- function(
  d, condition, syndrome_names, syndrome_values = NULL,
  specimen_types = NULL, sites = NULL
) {
  check_champs_data(d)
  condition <- check_valid_condition(condition, form = "long")
  # TODO: check syndrome_names conditions as well...

  rcd <- syndrome_values
  if (is.null(rcd))
    rcd <- syndrome_names
  names(rcd) <- syndrome_names

  # TODO: sites

  specimen_types <- check_valid_specimen_types(d, specimen_types)

  tmp <- d$dcd_pivot %>%
    filter(etiol == condition & infectious_cause == 1) %>%
    dplyr::mutate(
      champs_group_desc = recode(champs_group_desc, !!!rcd)) %>%
    group_by(champs_deid) %>%
    mutate(syndrome = paste(sort(unique(champs_group_desc)), collapse = ", ")) %>%
    group_by(syndrome, champs_deid) %>%
    tally()

  res <- tmp %>%
    group_by(syndrome) %>%
    tally() %>%
    arrange(-n)

  age_breakdown <- d$dmg %>%
    dplyr::filter(champs_deid %in% tmp$champs_deid) %>%
    dplyr::group_by(case_type_desc) %>%
    dplyr::tally()

  tmp <- dplyr::filter(d$tac_pivot, champs_deid %in% d$dcd$champs_deid)

  tmp2 <- tmp %>%
    dplyr::filter(
      specimen_type %in% specimen_types &
      condition == !!condition &
      result == "Positive")

  unique(tmp2$specimen_type)
  length(unique(tmp2$champs_deid)) / length(unique(tmp$champs_deid))

  tac_age_breakdown <- tmp2 %>%
    left_join(d$dmg, by = "champs_deid") %>%
    group_by(case_type_desc) %>%
    tally() %>%
    mutate(pct = 100 * n / sum(n))

  structure(list(
    table = res,
    age_breakdown = age_breakdown,
    cc_leading_to_death = list(
      numerator = sum(res$n),
      denominator = nrow(d$dmg),
      pct = sum(res$n) / nrow(d$dmg) * 100
    ),
    tac_age_breakdown = tac_age_breakdown
  ), class = c("champs_computed", "calc_syndrome_combinations"))
}


# #' TODO
# #'
# #' @param d A data object obtained from \code{\link{load_data}}.
# #' @param pathogen A string specifying the pathogen to count cases of
# #'   in the causal chain.
# #' @param condition A string specifying the condition to tabulate causal chain
# #'   presence for. To list all possibilities, see
# #'   \code{\link{valid_conditions}}.
# #' @param syndrome_values An optional vector of the same length of
# #'   \code{conditions} that specifies the mapping from the official
# #'   condition names to the values to use for reporting. See example.
# #' @param specimen_types An optional vector of specimen types to include in the
# #'   calculation. If not provided, all specimen types will be used. See
# #'   \code{\link{valid_specimen_types}}.
# #' @param sites An optional vector of site names to include in the tabulation.
# #'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
# #' @export
# #' @examples
# #' \dontrun{
# # calc_syndrome_combinations(d, pathogen = "Streptococcus pneumoniae")
# #' }
# calc_syndrome_combinations <- function(
#   d, pathogen, conditions, syndrome_values = NULL,
#   specimen_types = NULL, sites = NULL
# ) {
#   check_champs_data(d)

#   # TODO: sites

#   res <- d$dcd_pivot %>%
#     filter(etiol == pathogen) %>%
#     group_by(case_type_desc, acquired2) %>%
#     summarise(n = length(unique(champs_deid))) %>%
#     dplyr::ungroup() %>%
#     tidyr::complete(case_type_desc, acquired2, fill = list(n = 0)) %>%
#     tidyr::pivot_wider(names_from = acquired2, values_from = n) %>%
#     mutate(Total = Community + Facility) %>%
#     select(case_type_desc, Community, Facility, Total) %>%
#     arrange(case_type_desc)

#   structure(list(
#     table = res
#   ), class = c("champs_computed", "calc_syndrome_combinations"))
# }


#' Tabulate cases for a condition by age and syndrome
#'
#' @param d A data object obtained from \code{\link{load_data}}.
#' @param condition A string specifying the condition to count cases of
#'   in the causal chain.
#' @param condition A string specifying the condition to tabulate causal chain
#'   presence for. To list all possibilities, see
#'   \code{\link{valid_conditions}}.
#' @param syndrome_names A vector specifying conditions to tabulate
#'   combinations of presence of for the given \code{condition}
#' @param syndrome_values An optional vector of the same length of
#'   \code{syndrome_names} that specifies the mapping from the official
#'   condition names to the values to use for reporting. See example.
#' @param specimen_types An optional vector of specimen types to include in the
#'   calculation. If not provided, all specimen types will be used. See
#'   \code{\link{valid_specimen_types}}.
#' @param sites An optional vector of site names to include in the tabulation.
#'   If not provided, all sites will be used. See \code{\link{valid_sites}}.
#' @export
#' @examples
#' \dontrun{
#' calc_cc_by_age_syndrome(d,
#'   condition = "Streptococcus pneumoniae",
#'   syndrome_names = c(
#'     "Lower respiratory infections",
#'     "Meningitis/Encephalitis",
#'     "Neonatal sepsis",
#'     "Congenital infection"),
#'   syndrome_values = c(
#'     "Pneumonia",
#'     "Meningitis",
#'     "Sepsis",
#'     "Sepsis")
#' )
#' }
calc_cc_by_age_syndrome <- function(
  d, condition, syndrome_names, syndrome_values = NULL,
  specimen_types = NULL, sites = NULL
) {
  check_champs_data(d)
  condition <- check_valid_condition(condition, form = "long")
  # TODO: check syndrome_names conditions as well...

  rcd <- syndrome_values
  if (is.null(rcd))
    rcd <- syndrome_names
  names(rcd) <- syndrome_names

  # TODO: sites

  res <- d$dcd_pivot %>%
    filter(etiol == condition & infectious_cause == 1) %>%
    dplyr::mutate(
      champs_group_desc = recode(champs_group_desc, !!!rcd)) %>%
    group_by(champs_deid) %>%
    mutate(syndrome = paste(sort(unique(champs_group_desc)), collapse = ", ")) %>%
    group_by(syndrome, case_type_desc, champs_deid) %>%
    tally() %>%
    dplyr::ungroup() %>%
    tidyr::complete(case_type_desc, fill = list(n = 0)) %>%
    group_by(syndrome, case_type_desc) %>%
    tally() %>%
    tidyr::pivot_wider(names_from = syndrome, values_from = n,
      values_fill = list(n = 0)) %>%
    arrange(match(case_type_desc, levels(d$dcd_pivot$case_type_desc)))

  structure(list(
    table = res
  ), class = c("champs_computed", "calc_cc_by_age_syndrome"))
}
