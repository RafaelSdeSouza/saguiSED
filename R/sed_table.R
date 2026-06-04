#' Convert SAGUI regional photometry to a saguiSED table
#'
#' @param x A `sagui::extract_region_sed()` result, a wide regional photometry
#'   data frame, or a path to a CSV file.
#' @param filters Optional filter column names. If omitted, columns ending in
#'   `_err` or `_n_eff` and bookkeeping columns are ignored. If `filter_set` is
#'   provided, its filter order is used by default.
#' @param filter_set Optional `sagui_filter_set` object created by
#'   [sed_filter_set()].
#' @param unit Flux unit in the table. Supported values are `"Jy"`, `"uJy"`,
#'   `"nJy"`, and `"10nJy"`.
#' @param redshift Optional redshift attached as metadata.
#' @param id_col Region identifier column.
#' @param n_pix_col Pixel-count column.
#' @return A `sagui_sed_table` object.
#' @export
as_sagui_sed_table <- function(x,
                               filters = NULL,
                               filter_set = NULL,
                               unit = c("Jy", "uJy", "nJy", "10nJy"),
                               redshift = NULL,
                               id_col = "region",
                               n_pix_col = "n_pix") {
  unit <- match.arg(unit)

  if (is.character(x) && length(x) == 1) {
    x <- utils::read.csv(x, check.names = FALSE)
  } else if (is.list(x) && !is.data.frame(x) && !is.null(x$flux_wide)) {
    x <- x$flux_wide
  }

  if (!is.data.frame(x)) {
    stop("`x` must be a SAGUI SED result, a data frame, or a CSV path.", call. = FALSE)
  }
  if (!id_col %in% names(x)) {
    stop("Region identifier column not found: ", id_col, call. = FALSE)
  }

  if (!is.null(filter_set) && !inherits(filter_set, "sagui_filter_set")) {
    stop("`filter_set` must be created with `sed_filter_set()`.", call. = FALSE)
  }

  if (is.null(filters) && !is.null(filter_set)) {
    filters <- filter_set$filters
  }

  if (is.null(filters)) {
    drop_cols <- c(id_col, n_pix_col)
    filters <- setdiff(
      names(x),
      c(
        drop_cols,
        grep("(^err_|_err$|_n_eff$)", names(x), value = TRUE),
        grep("^flux_", names(x), value = TRUE)
      )
    )
  }

  missing_filters <- setdiff(filters, names(x))
  if (length(missing_filters)) {
    prefixed_flux <- paste0("flux_", filters)
    if (all(prefixed_flux %in% names(x))) {
      for (i in seq_along(filters)) {
        x[[filters[[i]]]] <- x[[prefixed_flux[[i]]]]
      }
      missing_filters <- character()
    }
  }
  if (length(missing_filters)) {
    stop(
      "Missing filter columns: ", paste(missing_filters, collapse = ", "),
      ". SAGUI-style `flux_<filter>` columns are also accepted.",
      call. = FALSE
    )
  }

  missing_err <- setdiff(paste0(filters, "_err"), names(x))
  if (length(missing_err)) {
    prefixed_err <- paste0("err_", filters)
    if (all(prefixed_err %in% names(x))) {
      for (i in seq_along(filters)) {
        x[[paste0(filters[[i]], "_err")]] <- x[[prefixed_err[[i]]]]
      }
      missing_err <- character()
    }
  }
  if (length(missing_err)) {
    stop(
      "Missing error columns required by saguiSED: ",
      paste(missing_err, collapse = ", "),
      ". SAGUI-style `err_<filter>` columns are also accepted.",
      call. = FALSE
    )
  }

  if (n_pix_col %in% names(x) && !"n_pix" %in% names(x)) {
    x$n_pix <- x[[n_pix_col]]
  }

  structure(
    list(
      flux_wide = x,
      filters = filters,
      filter_set = filter_set,
      unit = unit,
      redshift = redshift,
      id_col = id_col,
      n_pix_col = n_pix_col
    ),
    class = "sagui_sed_table"
  )
}
