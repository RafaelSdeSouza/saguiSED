#' Define a photometric filter set for SED fitting
#'
#' @param throughput A data frame or CSV path containing filter transmission
#'   curves.
#' @param filters Optional filter names/order to use. If omitted, the order of
#'   first appearance in `throughput` is used.
#' @param name Human-readable filter-set name.
#' @param wavelength_unit Unit of the input wavelength column.
#' @param filter_col Column containing filter names.
#' @param wave_col Column containing wavelengths.
#' @param throughput_col Column containing transmission values.
#' @param lambda_eff_col Optional column containing effective/central
#'   wavelengths. If absent, throughput-weighted central wavelengths are
#'   computed.
#' @return A `sagui_filter_set` object.
#' @export
sed_filter_set <- function(throughput,
                           filters = NULL,
                           name = "custom",
                           wavelength_unit = c("micron", "angstrom", "nm"),
                           filter_col = "filter",
                           wave_col = "wavelength",
                           throughput_col = "throughput",
                           lambda_eff_col = NULL) {
  wavelength_unit <- match.arg(wavelength_unit)

  if (is.character(throughput) && length(throughput) == 1) {
    if (!file.exists(throughput)) {
      stop("Throughput file not found: ", throughput, call. = FALSE)
    }
    throughput <- utils::read.csv(throughput, check.names = FALSE)
  }
  if (!is.data.frame(throughput)) {
    stop("`throughput` must be a data frame or CSV path.", call. = FALSE)
  }

  required <- c(filter_col, wave_col, throughput_col)
  missing <- setdiff(required, names(throughput))
  if (length(missing)) {
    stop("Missing throughput columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  wave <- as.numeric(throughput[[wave_col]])
  wave_um <- switch(
    wavelength_unit,
    micron = wave,
    angstrom = wave / 1e4,
    nm = wave / 1e3
  )

  tbl <- data.frame(
    filter = as.character(throughput[[filter_col]]),
    wavelength_um = wave_um,
    throughput = as.numeric(throughput[[throughput_col]]),
    stringsAsFactors = FALSE
  )

  if (!is.null(lambda_eff_col)) {
    if (!lambda_eff_col %in% names(throughput)) {
      stop("`lambda_eff_col` not found: ", lambda_eff_col, call. = FALSE)
    }
    lambda_eff <- as.numeric(throughput[[lambda_eff_col]])
    tbl$lambda_eff_um <- switch(
      wavelength_unit,
      micron = lambda_eff,
      angstrom = lambda_eff / 1e4,
      nm = lambda_eff / 1e3
    )
  }

  tbl <- tbl[is.finite(tbl$wavelength_um) & is.finite(tbl$throughput), , drop = FALSE]
  if (!nrow(tbl)) {
    stop("No finite throughput rows remain after parsing.", call. = FALSE)
  }

  if (is.null(filters)) {
    filters <- unique(tbl$filter)
  } else {
    filters <- as.character(filters)
  }
  missing_filters <- setdiff(filters, unique(tbl$filter))
  if (length(missing_filters)) {
    stop("Filters missing from throughput table: ", paste(missing_filters, collapse = ", "), call. = FALSE)
  }

  tbl <- tbl[tbl$filter %in% filters, , drop = FALSE]
  tbl$filter <- factor(tbl$filter, levels = filters)
  tbl <- tbl[order(tbl$filter, tbl$wavelength_um), , drop = FALSE]
  tbl$filter <- as.character(tbl$filter)

  lambda_eff <- vapply(filters, function(filt) {
    curve <- tbl[tbl$filter == filt, , drop = FALSE]
    if ("lambda_eff_um" %in% names(curve) && any(is.finite(curve$lambda_eff_um))) {
      return(stats::median(curve$lambda_eff_um[is.finite(curve$lambda_eff_um)]))
    }
    weights <- pmax(curve$throughput, 0)
    if (sum(weights, na.rm = TRUE) > 0) {
      stats::weighted.mean(curve$wavelength_um, weights, na.rm = TRUE)
    } else {
      stats::median(curve$wavelength_um, na.rm = TRUE)
    }
  }, numeric(1))

  structure(
    list(
      name = name,
      filters = filters,
      throughput = tbl,
      lambda_eff_um = lambda_eff,
      wavelength_unit = "micron"
    ),
    class = "sagui_filter_set"
  )
}

#' JWST/NIRCam filter set bundled with saguiSED
#'
#' @param filters Optional NIRCam filters/order.
#' @return A `sagui_filter_set` object.
#' @export
jwst_nircam_filter_set <- function(filters = NULL) {
  path <- system.file("extdata", "throughput_nircam.csv", package = "saguiSED")
  if (!nzchar(path)) {
    stop(
      "Could not find bundled JWST/NIRCam throughput table. ",
      "Use `sed_filter_set()` with an explicit throughput table instead.",
      call. = FALSE
    )
  }
  sed_filter_set(
    throughput = path,
    filters = filters,
    name = "JWST/NIRCam",
    wavelength_unit = "micron",
    filter_col = "filter",
    wave_col = "wavelength",
    throughput_col = "throughput",
    lambda_eff_col = "lambda_c"
  )
}

#' Rubin/LSST filter set bundled with saguiSED
#'
#' The bundled LSST response curves are generated from the `speclite` filter
#' definitions using the current `lsst2023` Rubin throughput assumptions.
#'
#' @param filters Optional LSST filters/order. Defaults to `u,g,r,i,z,y`.
#' @return A `sagui_filter_set` object.
#' @export
lsst_filter_set <- function(filters = NULL) {
  path <- system.file("extdata", "throughput_lsst.csv", package = "saguiSED")
  if (!nzchar(path)) {
    stop(
      "Could not find bundled LSST throughput table. ",
      "Use `sed_filter_set()` with an explicit throughput table instead.",
      call. = FALSE
    )
  }
  if (is.null(filters)) {
    filters <- c("u", "g", "r", "i", "z", "y")
  }
  sed_filter_set(
    throughput = path,
    filters = filters,
    name = "Rubin/LSST",
    wavelength_unit = "micron",
    filter_col = "filter",
    wave_col = "wavelength",
    throughput_col = "throughput",
    lambda_eff_col = "lambda_c"
  )
}

write_filter_set_csv <- function(filter_set, path) {
  if (!inherits(filter_set, "sagui_filter_set")) {
    stop("`filter_set` must be created with `sed_filter_set()`.", call. = FALSE)
  }
  tbl <- filter_set$throughput
  lambda_tbl <- data.frame(
    filter = names(filter_set$lambda_eff_um),
    lambda_eff_um = unname(filter_set$lambda_eff_um),
    stringsAsFactors = FALSE
  )
  tbl <- merge(tbl[, c("filter", "wavelength_um", "throughput")], lambda_tbl, by = "filter", all.x = TRUE)
  tbl$filter <- factor(tbl$filter, levels = filter_set$filters)
  tbl <- tbl[order(tbl$filter, tbl$wavelength_um), , drop = FALSE]
  tbl$filter <- as.character(tbl$filter)
  utils::write.csv(tbl, path, row.names = FALSE)
  invisible(path)
}
