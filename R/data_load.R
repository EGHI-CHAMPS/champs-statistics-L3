#' @importFrom utils str
#' @export
print.champs_data <- function(x, ...) {
  utils::str(x, 1)
}

fix_names <- function(x) {
  nms <- names(x)
  nms <- tolower(nms)
  nms <- gsub(" ", "_", nms)
  names(x) <- nms
  x
}

#' List rules for creating new variables in \code{\link{load_data}}
#' @export
get_rules <- function() {
  rules
}

# list.files("../data_l3")

# load_data(
#   dmg_path = "../data_l3/CHAMPS_L3_basic_demographics_V3.csv",
#   lab_path = "../data_l3/CHAMPS_L3_lab_results_V3.csv",
#   dcd_path = "../data_l3/CHAMPS_L3_decode_results_V3.csv",
#   tac_path = "../data_l3/CHAMPS_L3_tac_results_V3.csv"
# )

sm <- suppressMessages

#' Load and pre-process data from relevant CHAMPS database tables
#'
#' @param data_dir Path to directory containing L3 csv files.
#'   Note that if this is provided, the additional `path`
#'   parameters are not required, but csv files must exist
#'   with the naming scheme as described for these files.
#' @param dmg_path Path to L3 demographics csv file (typically
#'   named "CHAMPS_L3_basic_demographics_V[x].csv").
#' @param lab_path Path to L3 lab results csv file (typically
#'   named "CHAMPS_L3_lab_results_V[x].csv").
#' @param dcd_path Path to L3 DeCoDe results csv file (typically
#'   named "CHAMPS_L3_decode_results_V[x].csv").
#' @param tac_path Path to L3 TAC results csv file (typically 
#'   namde "CHAMPS_L3_tac_results_V[x].csv").
#' @param rules Optional list of rules for calculating new variables
#'   (see \code{\link{get_rules}})
#' @note The excel files must be fresh exports from the database not have had
#'   any processing, such as additional sheets added to them. One exception is
#'   that it is assumed that the TAC results file (indicated by
#'   \code{tac_res_path}) is assumed to have indeterminate results removed,
#'   although the logic can easily be changed so that this is not required. This
#'   function and all functions in this package rely on the input data coming
#'   directly from the database with no variable names changed, etc.
#' @export
#' @importFrom readr read_csv
load_data <- function(
  data_dir = NULL,
  dmg_path = NULL,
  lab_path = NULL,
  dcd_path = NULL,
  tac_path = NULL,
  rules = NULL
) {
  get_path <- function(x, ff, data_dir) {
    idx <- which(grepl(x, ff))
    if (length(idx) > 0) {
      if (length(idx) > 1) {
        message("There were ", length(idx), " files that match ",
          x, "... using the first: ", ff[idx[1]])
        idx <- idx[1]
      }
      return(file.path(data_dir, ff[idx]))
    } else {
      stop("Couldn't find file with name matching ", x, " in ",
        "provided dat_dir: ", data_dir, ".", call. = FALSE)
    }
  }

  if (!is.null(data_dir)) {
    ff <- list.files(data_dir)
    dmg_path <- get_path("CHAMPS_L3_basic_demographics",
      ff, data_dir)
    lab_path <- get_path("CHAMPS_L3_lab_results",
      ff, data_dir)
    dcd_path <- get_path("CHAMPS_L3_decode_results",
      ff, data_dir)
    tac_path <- get_path("CHAMPS_L3_tac_results",
      ff, data_dir)
  } else {
    if (is.null(dmg_path) || is.null(lab_path) ||
      is.null(dcd_path) || is.null(tac_path))
    stop("Must provide paths for all files, or provide the ",
      "path that these files reside in, 'data_dir'.",
        call. = FALSE)
  }

  rd <- function(path) {
    tmp <- sm(readr::read_csv(path, guess_max = 100000)) %>%
      fix_names()
    attr(tmp, "spec") <- NULL
    tmp
  }

  message("Reading demographics data: ", dmg_path, "...")
  dmg <- rd(dmg_path)

  # temporarily add new case_type_desc variable (until subcat is in L3)
  # TODO: fix when subcat is in L3
  dmg$case_type_desc2 <- as.character(dmg$case_type_desc)
  idx <- which(dmg$case_type_desc == "Early Neonate (1 to 6 days)")
  dmg$case_type_desc2[idx] <- "Early Neonate (72+hrs to 6 days)"
  dmg$case_type_desc2[idx][dmg$age_days[idx] < 3] <- "Early Neonate (24-72 hours)"
  dmg$case_type_desc2 <- factor(dmg$case_type_desc2,
    levels = c(
      "Stillbirth",
      "Death in the first 24 hours",
      "Early Neonate (24-72 hours)",
      "Early Neonate (72+hrs to 6 days)",
      "Late Neonate (7 to 27 days)",
      "Infant (28 days to less than 12 months)",
      "Child (12 months to less than 60 Months)"))

  # make up fake PMI variable for now so we can test
  # TODO: fix when PMI is in L3
  dmg$calc_postmortem_hrs <- rep(0:40, ceiling(nrow(dmg) / 41))[1:nrow(dmg)]

  message("Reading lab data: ", lab_path, "...")
  lab <- rd(lab_path)

  message("Reading DeCoDe results: ", dcd_path, "...")
  dcd <- rd(dcd_path)

  message("Reading TAC results: ", tac_path, "...")
  tac <- rd(tac_path)
    # dplyr::filter(mits_flag != 0) %>%
    # dplyr::rename(champsid = "champs_id")

  # TODO: add checks here to make sure all variable names exist that will be
  # used in the code

  dmg <- dmg %>%
    mutate(case_type_desc = forcats::fct_relevel(case_type_desc, c(
      "Stillbirth",
      "Death in the first 24 hours",
      "Early Neonate (1 to 6 days)",
      "Late Neonate (7 to 27 days)",
      "Infant (28 days to less than 12 months)",
      "Child (12 months to less than 60 Months)"
    ))) %>%
    mutate(
      site = forcats::fct_relevel(recode(site_iso_code,
        BD = "Bangladesh",
        ET = "Ethiopia",
        KE = "Kenya",
        ML = "Mali",
        MZ = "Mozambique",
        SL = "Sierra Leone",
        ZA = "South Africa"
      ),
      c("Bangladesh", "Kenya", "Mali", "Mozambique",
        "South Africa", "Ethiopia", "Sierra Leone"))
    )

  # site_lookup <- dcd_meas %>%
  #   dplyr::select(champsid, site_name) %>%
  #   dplyr::rename(site_name_ofcl = "site_name") %>%
  #   dplyr::distinct() %>%
  #   dplyr::arrange(site_name_ofcl) %>%
  #   dplyr::mutate(site_name_ofcl = factor(site_name_ofcl))

  # extra processing
  dcd_join <- add_variables(dcd, dmg)

  message("Building TAC pivot...")
  # build TAC pivot
  tac_pivot <- build_tac_pivot(tac)


  # TODO: add etiol_tac column to dcd_pivot that maps...
  message("Building DeCoDe pivot...")
  # incl are variables to pull into pivot table from demographics
  incl <- c("case_type_desc", "case_type_desc2",
    "acquired2", "infectious_cause", "site")
  dcd_pivot <- build_dcd_pivot(dcd_join, incl)

  structure(list(
    dmg = dmg,
    lab = lab,
    dcd = dcd,
    tac = tac,
    dcd_pivot = dcd_pivot,
    tac_pivot = tac_pivot,
    dcd_join = dcd_join
  ), class = c("champs_data", "list"))
}

# #' Load and pre-process data from relevant CHAMPS database tables
# #'
# #' @param dcd_res_path Path to Excel export of the DeCoDe results database
# #'   table ("vw_decode_results_classification").
# #' @param dcd_meas_path Path to Excel export of the DeCoDe measures database
# #'   table ("lk_vw_dm_all.xlsx").
# #' @param tac_res_path Path to Excel export of the TAC results database table
# #'   ("vw_TACResult_By_ChampsID").
# #' @param tac_lab_path Path to Excel export of the TAC laboratory results
# #'   database pivot table ("TacLaboratoryResult_Pivot").
# #' @param rules Optional list of rules for calculating new variables
# #'   (see \code{\link{get_rules}})
# #' @note The excel files must be fresh exports from the database not have had
# #'   any processing, such as additional sheets added to them. One exception is
# #'   that it is assumed that the TAC results file (indicated by
# #'   \code{tac_res_path}) is assumed to have indeterminate results removed,
# #'   although the logic can easily be changed so that this is not required. This
# #'   function and all functions in this package rely on the input data coming
# #'   directly from the database with no variable names changed, etc.
# #' @export
# load_data_old <- function(
#   dcd_res_path,
#   dcd_meas_path,
#   tac_res_path,
#   tac_lab_path,
#   rules = NULL
# ) {
#   message("Reading TAC Laboratory Results: ", tac_lab_path, "...")
#   tac_lab <- readxl::read_xlsx(tac_lab_path, guess_max = 100000) %>%
#     fix_names() %>%
#     dplyr::rename(champsid = "champs_id")

#   message("Reading TAC Results: ", tac_res_path, "...")
#   tac_res <- readxl::read_xlsx(tac_res_path, guess_max = 100000) %>%
#     fix_names() %>%
#     dplyr::rename(champsid = "champs_id")

#   message("Reading DeCoDe results: ", dcd_res_path, "...")
#   dcd_res <- readxl::read_xlsx(dcd_res_path, guess_max = 100000) %>%
#     fix_names()

#   message("Reading DeCoDe measurements: ", dcd_meas_path, "...")
#   dcd_meas <- readxl::read_xlsx(dcd_meas_path, guess_max = 100000) %>%
#     fix_names() %>%
#     dplyr::filter(mits_flag != 0) %>%
#     dplyr::rename(champsid = "champs_id")

#   # TODO: add checks here to make sure all variable names exist that will be
#   # used in the code

#   case_lookup <- dcd_res %>%
#     dplyr::select(champsid, casetype) %>%
#     dplyr::rename(casetype_ofcl = "casetype") %>%
#     dplyr::mutate(
#       casetype_ofcl = dplyr::recode(casetype_ofcl,
#         Infant = "Infant (28 days to less than 12 months)",
#         Child = "Child (12 months to less than 60 Months)"
#       ),
#       casetype_ofcl = forcats::fct_relevel(casetype_ofcl, c(
#         "Stillbirth",
#         "Death in the first 24 hours",
#         "Early Neonate (1 to 6 days)",
#         "Late Neonate (7 to 27 days)",
#         "Infant (28 days to less than 12 months)",
#         "Child (12 months to less than 60 Months)"
#       ))
#     )

#   site_lookup <- dcd_meas %>%
#     dplyr::select(champsid, site_name) %>%
#     dplyr::rename(site_name_ofcl = "site_name") %>%
#     dplyr::distinct() %>%
#     dplyr::arrange(site_name_ofcl) %>%
#     dplyr::mutate(site_name_ofcl = factor(site_name_ofcl))

#   # extra processing
#   dcd_join <- add_variables_old(dcd_res, dcd_meas)

#   structure(list(
#     tac_lab = tac_lab,
#     tac_res = tac_res,
#     dcd_res = dcd_res,
#     dcd_meas = dcd_meas,
#     dcd_join = dcd_join,
#     case_lookup = case_lookup,
#     site_lookup = site_lookup
#   ), class = "champs_data")
# }
