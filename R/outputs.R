get_table <- function(
  x, digits = 1,
  margin_pct = TRUE,
  margin_denom = FALSE,
  inside_pct = FALSE,
  inside_denom = FALSE
) {
  if (is.null(x$numerator) || is.null(x$denominator))
    stop("Currently only tables can be output that have both", " a numerator and denominator.")

  nn <- x$numerator
  dd <- x$denominator

  nms <- colnames(nn)
  ncol <- length(nms)
  col_tots <- dd[nrow(dd), ]
  nms <- paste0(nms, " (n=", col_tots, ")")
  nms[ncol] <- "Total"

  margin_pct <- TRUE
  digits <- 1
  tbl_pars <- lapply(seq_along(nms), function(ii) {
    res <- unname(nn[, ii])
    if (nms[ii] == "Total") {
      if (margin_denom) {
        res <- paste0(res, "/", dd[,ii])
      }
      if (margin_pct) {
        pct <- round(100 * nn[, ii] / dd[, ii], digits)
        res <- paste0(res, " (", pct, "%)")
      }
    } else {
      n <- length(res)
      ss <- seq_len(n - 1)
      if (inside_denom) {
        res[ss] <- paste0(res[ss],
          ifelse(dd[ss, ii] == 0, "", "/"),
          ifelse(dd[ss, ii] == 0, "", dd[ss, ii]))
      }
      if (margin_denom) {
        res[n] <- paste0(res[n],
          ifelse(dd[n, ii] == 0, "", "/"),
          ifelse(dd[n, ii] == 0, "", dd[n, ii]))
      }
      if (margin_pct) {
        pct <- round(100 * nn[n, ii] / dd[n, ii], digits)
        res[n] <- paste0(res[n], " (", pct, "%)")
      }
      if (inside_pct) {
        pct <- round(100 * nn[ss, ii] / dd[ss, ii], digits)
        res[ss] <- paste0(res[ss], " (", pct, "%)")
      }
    }
    res
  })

  names(tbl_pars) <- nms

  rnms <- rownames(nn)
  nrow <- length(rnms)
  rnms[nrow] <- "Total"
  row_tots <- dd[, ncol(dd)]
  rnms <- paste0(rnms, " (n=", row_tots, ")")

  tbl_pars <- c(list(`Case Type` = rnms), tbl_pars)

  do.call(tibble::tibble, tbl_pars)
}

#' Write an html table of a computed CHAMPS object
#'
#' @param x A computed CHAMPS object.
#' @param margin_pct Should percentages be printed in the margins?
#' @param margin_denom Should denominators be printed in the margins?
#' @param inside_pct Should percentages be printed inside the table?
#' @param inside_denom Should denominators be printed inside the table?
#' @importFrom htmlTable htmlTable
#' @export
write_html_table <- function(x,
  margin_pct = TRUE,
  margin_denom = FALSE,
  inside_pct = FALSE,
  inside_denom = FALSE
) {
  check_champs_object(x)

  tbl <- get_table(x,
    margin_pct = margin_pct,
    margin_denom = margin_denom,
    inside_pct = inside_pct,
    inside_denom = inside_denom
  )
  htmlTable::htmlTable(tbl,
    rnames = FALSE,
    total = TRUE,
    align = c("l", rep("c", ncol(tbl) - 1)),
    css.table = "font-family: Lato",
    css.cell = "padding: 5px;"
  )
}

# doc <- officer::read_pptx()

#' Print a computed CHAMPS table to a PowerPoint presentation
#'
#' @param x A computed CHAMPS object.
#' @param doc A PowerPoint presentation object.
#' @param margin_pct Should percentages be printed in the margins?
#' @param margin_denom Should denominators be printed in the margins?
#' @param inside_pct Should percentages be printed inside the table?
#' @param inside_denom Should denominators be printed inside the table?
#' @importFrom officer add_slide ph_with ph_location_type
#' @export
write_ppt_slide <- function(x, doc,
  margin_pct = TRUE,
  margin_denom = FALSE,
  inside_pct = FALSE,
  inside_denom = FALSE
) {
  check_champs_object(x)

  tbl <- get_table(x,
    margin_pct = margin_pct,
    margin_denom = margin_denom,
    inside_pct = inside_pct,
    inside_denom = inside_denom
  )

  doc <- officer::add_slide(doc)
  officer::ph_with(x = doc, value = tbl,
    location = officer::ph_location_type(type = "body"))
}

#' Plot the margins (row and column totals) of a computed CHAMPS table
#'
#' @param x A computed CHAMPS object.
#' @importFrom ggplot2 ggplot aes geom_text theme_bw theme element_text labs ylim
#' @importFrom cowplot plot_grid
#' @importFrom utils tail
#' @export
plot_margins <- function(x) {
  check_champs_object(x)

  get_data <- function(nx, dx, nms, which) {
    tibble::tibble(
      pct = as.vector(unname(100 * nx / dx)),
      var = forcats::fct_relevel(nms, nms),
      txt = paste(nx, dx, sep = "/"),
      which = which
    )
  }

  nn <- x$numerator
  dd <- x$denominator

  nms <- colnames(dd)
  nms[length(nms)] <- "Total"
  nx <- utils::tail(nn, 1)
  dx <- utils::tail(dd, 1)

  d1 <- get_data(nx, dx, nms, "Site")

  nms <- rownames(dd)
  nms[length(nms)] <- "Total"
  nx <- nn[, ncol(nn)]
  dx <- dd[, ncol(dd)]

  d2 <- get_data(nx, dx, nms, "Case Type")

  ylims <- range(c(d1$pct, d2$pct, 0))
  ylims[2] <- ylims[2] + 0.3

  p1 <- ggplot2::ggplot(d1, ggplot2::aes(var, pct, label = txt)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(vjust = 0, nudge_y = 0.1) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(y = "Percent", x = NULL, title = "Site") +
    ggplot2::ylim(ylims)

  p2 <- ggplot2::ggplot(d2, ggplot2::aes(var, pct, label = txt)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(vjust = 0, nudge_y = 0.1) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(y = "Percent", x = NULL, title = "Case Type") +
    ggplot2::ylim(ylims)

  cowplot::plot_grid(p1, p2, align = "h", axis = "b")
}
