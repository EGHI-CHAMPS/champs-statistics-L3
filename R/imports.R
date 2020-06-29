#' @importFrom dplyr filter left_join setdiff select mutate arrange rename
#'   recode distinct group_by tally right_join summarise
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect starts_with matches one_of
#' @importFrom forcats fct_drop fct_relevel
#' @importFrom stats xtabs addmargins
#' @importFrom tibble tibble
#' @importFrom rlang sym :=
NULL

#' Pipe
#'
#' Use the pipe function, \code{\%>\%} to turn function composition
#' into a series of imperative statements.
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @param lhs,rhs an object and a function to apply to it
NULL

#' tac_table
#'
"tac_table"

utils::globalVariables(c(
  "CISum", "DiarrhealDSum", "HIVSum", "LRISum", "MESum", "MalariaSum",
  "MeaslesSum", "NeoSepsisSum", "OtherInfectionSum", "SepsisSum",
  "SyphilisSum", "TBSum", "URISum", "calc_location", "calc_postmortem_hrs",
  "case_type", "case_type_desc", "casetype", "casetype_ofcl", "champsid",
  "has_np", "has_other", "hosp_los_24h", "hosp_los_24h2", "hosp_los_48h",
  "hosp_los_48h2", "infectious_cause", "mits_flag", "n", "name",
  "p2", "pct", "pmirange", "result", "site_name", "site_name_ofcl",
  "specimen_type", "tbl", "txt", "value", "var",
  "Community", "Facility", "Total", "acquired2", "assay",
  "case_type_desc2", "champs_deid", "champs_group_desc", "code",
  "count", "desc", "eshort", "etiol", "etiolgies",
  "location_of_death", "pathogen", "percent", "percent_label",
  "site", "site_iso_code", "syndrome", "tac_table", "xuse", "y"
))
