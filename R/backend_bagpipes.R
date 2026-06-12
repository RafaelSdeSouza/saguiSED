backend_script <- function() {
  system.file("python", "bagpipes_fit_regions.py", package = "saguiSED", mustWork = TRUE)
}

default_nircam_throughput <- function() {
  path <- system.file("extdata", "throughput_nircam.csv", package = "saguiSED")
  if (!nzchar(path)) {
    stop(
      "Could not find bundled JWST/NIRCam throughput table. ",
      "Pass `throughput_path` explicitly.",
      call. = FALSE
    )
  }
  path
}

#' Check whether the Bagpipes Python backend can be imported
#'
#' @param python Python executable.
#' @return Invisibly returns `TRUE` if Bagpipes can be imported.
#' @export
check_bagpipes <- function(python = Sys.which("python3")) {
  if (!nzchar(python)) {
    stop("Could not find `python3` on PATH.", call. = FALSE)
  }

  code <- "import bagpipes; print(getattr(bagpipes, '__version__', 'available'))"
  old <- Sys.getenv(c("NUMBA_DISABLE_JIT", "MPLCONFIGDIR"), unset = NA)
  on.exit({
    if (is.na(old[["NUMBA_DISABLE_JIT"]])) Sys.unsetenv("NUMBA_DISABLE_JIT") else Sys.setenv(NUMBA_DISABLE_JIT = old[["NUMBA_DISABLE_JIT"]])
    if (is.na(old[["MPLCONFIGDIR"]])) Sys.unsetenv("MPLCONFIGDIR") else Sys.setenv(MPLCONFIGDIR = old[["MPLCONFIGDIR"]])
  }, add = TRUE)
  Sys.setenv(NUMBA_DISABLE_JIT = "1", MPLCONFIGDIR = tempdir())

  status <- system2(python, c("-c", shQuote(code)), stdout = TRUE, stderr = TRUE)
  ok <- identical(attr(status, "status"), NULL)
  if (!ok) {
    stop("Bagpipes backend is not available:\n", paste(status, collapse = "\n"), call. = FALSE)
  }

  message("Bagpipes backend available: ", paste(status, collapse = " "))
  invisible(TRUE)
}

#' Fit regional SEDs with an optional backend
#'
#' @param sed A `sagui_sed_table` object returned by [saguiSED::as_sagui_sed_table()].
#' @param backend Backend name. Currently `"bagpipes"`.
#' @param model Model configuration, e.g. [saguiSED::bagpipes_model()].
#' @param redshift Redshift. Overrides the redshift stored in `sed`.
#' @param regions Optional region IDs to fit. If omitted, all regions are fit.
#' @param filter_set Optional `sagui_filter_set` object. If omitted, the filter
#'   set stored in `sed` is used.
#' @param throughput_path Deprecated compatibility shortcut. Prefer
#'   [sed_filter_set()] and `filter_set`.
#' @param out_dir Output directory.
#' @param python Python executable.
#' @param overwrite Logical; overwrite existing backend outputs.
#' @return A `sagui_sed_fit` object.
#' @export
fit_region_seds <- function(sed,
                            backend = c("bagpipes"),
                            model = bagpipes_model(),
                            redshift = NULL,
                            regions = NULL,
                            filter_set = NULL,
                            throughput_path = NULL,
                            out_dir = "sagui_sedfit",
                            python = Sys.which("python3"),
                            overwrite = TRUE) {
  backend <- match.arg(backend)
  if (!inherits(sed, "sagui_sed_table")) {
    stop("`sed` must be created with `as_sagui_sed_table()`.", call. = FALSE)
  }
  if (!inherits(model, "sagui_bagpipes_model")) {
    stop("`model` must be created with `bagpipes_model()` for the Bagpipes backend.", call. = FALSE)
  }
  if (!nzchar(python)) {
    stop("Could not find `python3` on PATH.", call. = FALSE)
  }

  z <- redshift
  if (is.null(z)) z <- sed$redshift
  if (is.null(z) || !is.numeric(z) || length(z) != 1 || !is.finite(z)) {
    stop("A finite scalar `redshift` is required.", call. = FALSE)
  }

  if (is.null(filter_set)) {
    filter_set <- sed$filter_set
  }
  if (is.null(filter_set) && !is.null(throughput_path)) {
    warning(
      "`throughput_path` is deprecated. Use `filter_set = sed_filter_set(...)` instead.",
      call. = FALSE
    )
    filter_set <- sed_filter_set(
      throughput = throughput_path,
      filters = sed$filters,
      name = "custom",
      wavelength_unit = "micron"
    )
  }
  if (!inherits(filter_set, "sagui_filter_set")) {
    stop(
      "A filter set is required. Pass `filter_set = sed_filter_set(...)`, ",
      "or attach one with `as_sagui_sed_table(..., filter_set = ...)`.",
      call. = FALSE
    )
  }
  missing_filters <- setdiff(sed$filters, filter_set$filters)
  if (length(missing_filters)) {
    stop("SED filters missing from filter set: ", paste(missing_filters, collapse = ", "), call. = FALSE)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  input_path <- file.path(out_dir, "regional_photometry_input.csv")
  filter_path <- file.path(out_dir, "filter_set.csv")
  utils::write.csv(sed$flux_wide, input_path, row.names = FALSE)
  write_filter_set_csv(filter_set, filter_path)

  args <- c(
    backend_script(),
    "--photometry-csv", input_path,
    "--throughput-csv", filter_path,
    "--out-dir", out_dir,
    "--filters", paste(sed$filters, collapse = ","),
    "--unit", sed$unit,
    "--redshift", as.character(z),
    "--systematic-frac", as.character(model$systematic_frac),
    "--sfh", model$sfh,
    "--dust", model$dust,
    "--metallicity", model$metallicity,
    "--fixed-metallicity", as.character(model$fixed_metallicity)
  )
  if (!is.null(regions)) {
    args <- c(args, "--regions", paste(regions, collapse = ","))
  }
  if (isTRUE(overwrite)) {
    args <- c(args, "--overwrite")
  }

  old <- Sys.getenv(c("NUMBA_DISABLE_JIT", "MPLCONFIGDIR"), unset = NA)
  on.exit({
    if (is.na(old[["NUMBA_DISABLE_JIT"]])) Sys.unsetenv("NUMBA_DISABLE_JIT") else Sys.setenv(NUMBA_DISABLE_JIT = old[["NUMBA_DISABLE_JIT"]])
    if (is.na(old[["MPLCONFIGDIR"]])) Sys.unsetenv("MPLCONFIGDIR") else Sys.setenv(MPLCONFIGDIR = old[["MPLCONFIGDIR"]])
  }, add = TRUE)
  Sys.setenv(NUMBA_DISABLE_JIT = "1", MPLCONFIGDIR = tempdir())

  output <- system2(python, args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    stop("Bagpipes backend failed:\n", paste(output, collapse = "\n"), call. = FALSE)
  }

  summary_path <- file.path(out_dir, "region_fit_summary.csv")
  photom_path <- file.path(out_dir, "model_photometry.csv")
  spectrum_path <- file.path(out_dir, "model_spectrum.csv")
  if (!file.exists(summary_path)) {
    stop("Backend completed but did not write: ", summary_path, call. = FALSE)
  }

  summary <- utils::read.csv(summary_path, check.names = FALSE)
  summary <- add_property_aliases(summary)
  model_photometry <- if (file.exists(photom_path)) {
    utils::read.csv(photom_path, check.names = FALSE)
  } else {
    data.frame()
  }
  model_spectrum <- if (file.exists(spectrum_path)) {
    utils::read.csv(spectrum_path, check.names = FALSE)
  } else {
    data.frame()
  }

  structure(
    list(
      backend = backend,
      model = model,
      redshift = z,
      filter_set = filter_set,
      out_dir = normalizePath(out_dir, mustWork = FALSE),
      summary = summary,
      model_photometry = model_photometry,
      model_spectrum = model_spectrum,
      log = output
    ),
    class = "sagui_sed_fit"
  )
}
