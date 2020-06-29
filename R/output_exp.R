#' Experimental html table
#' 
#' @param x Object of class "cc_allcases_by_site_casetype".
#' @importFrom grDevices colorRampPalette
#' @importFrom formattable formatter
#' @importFrom formattable style
#' @importFrom kableExtra kable
#' @importFrom kableExtra kable_styling
#' @importFrom kableExtra add_header_above
#' @importFrom kableExtra row_spec
#' @importFrom kableExtra column_spec
#' @importFrom kableExtra save_kable
#' @importFrom utils browseURL
#' @export
write_html_table_experimental <- function(x) {
  if (!inherits(x, "cc_allcases_by_site_casetype"))
    stop("Experimental table only works for cc_allcases_by_site_casetype")

  nn <- as.data.frame.matrix(x$numerator)
  dd <- as.data.frame.matrix(x$denominator)
  nms <- names(nn)
  nms[nms == "Sum"] <- "Total"
  nr <- nrow(nn)
  nc <- ncol(nn)

  a <- unlist(nn) / unlist(dd)
  a[is.nan(a)] <- 0
  rng <- range(a)
  cols <- grDevices::colorRampPalette(c("white", "red"))(101)

  tbl <- lapply(seq_len(nc), function(ii) {
    lapply(seq_len(nr), function(jj) {
      list(nn = nn[jj, ii], dd = dd[jj, ii])
    })
  })
  names(tbl) <- nms
  tbl <- dplyr::as_tibble(tbl)
  rnms <- rownames(nn)
  rnms[length(rnms)] <- "Total"
  tbl <- dplyr::bind_cols(dplyr::tibble(`Case Type` = rnms), tbl)

  fn <- function(x) {
    sapply(x, function(a) paste0(a$nn, "/", a$dd))
  }

  fmt <- formattable::formatter("span",
    style = function(x) {
      nns <- sapply(x, function(a) a$nn)
      dds <- sapply(x, function(a) a$dd)
      vals <- floor(100 * (nns / dds) / rng[2]) + 1
      vals[is.nan(vals)] <- 1
      formattable::style(
        display = "block",
        padding = "0 4px",
        color = "black !important",
        `text-align` = "right",
        `border-radius` = "2px",
        `background-color` = cols[vals])
    },
    `data-toggle` = "popover",
    `data-trigger` = "hover",
    `data-placement` = "right",
    `data-content` = function(x) {
      nns <- sapply(x, function(a) a$nn)
      dds <- sapply(x, function(a) a$dd)
      vals <- round(100 * (nns / dds), 1)
      vals[is.nan(vals)] <- 0
      paste0(vals, "%")
    },
    fn
  )

  tbl %>%
    dplyr::mutate_at(nms, fmt) %>%
    kableExtra::kable("html", escape = FALSE) %>%
    kableExtra::kable_styling(c("hover", "condensed"), full_width = FALSE) %>%
    kableExtra::column_spec(2:ncol(tbl), width = "120px") %>%
    kableExtra::add_header_above(c(" ", "Site" = ncol(tbl) - 2, " ")) %>%
    kableExtra::row_spec(0, align = "c") %>%
    kableExtra::row_spec(nrow(tbl), bold = TRUE,
      extra_css = "border-top: 1px solid rgb(119, 119, 119);") %>%
    kableExtra::column_spec(ncol(tbl), bold = TRUE, border_left = TRUE) %>%
    kableExtra::column_spec(1, border_right = TRUE) %>%
    print_and_browse()
}

print_and_browse <- function(x) {
  x <- paste0("<script>
$(document).ready(function(){
    $('[data-toggle=\"popover\"]').popover(); 
});
</script>
", x)
  tf <- tempfile(fileext = ".html")
  kableExtra::save_kable(x, file = tf, self_contained = TRUE)
  utils::browseURL(tf)
}
