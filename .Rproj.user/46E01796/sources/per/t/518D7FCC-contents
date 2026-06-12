#' Paint region-level SED-fit properties back onto an image grid
#'
#' @param labels Region-label matrix or a SAGUI segmentation object containing
#'   `cluster_map`.
#' @param fit A `sagui_sed_fit` object.
#' @param properties Fit-summary columns to paint.
#' @return A named list of matrices.
#' @export
paint_sed_properties <- function(labels,
                                 fit,
                                 properties = c(
                                   "logMformed",
                                   "SFR",
                                   "logsSFRformed",
                                   "age_gyr",
                                   "logZ_Zsun",
                                   "Av",
                                   "fit_rms_sigma"
                                 )) {
  if (is.list(labels) && !is.null(labels$cluster_map)) {
    labels <- labels$cluster_map
  }
  if (!is.matrix(labels)) {
    stop("`labels` must be a matrix or a SAGUI segmentation object.", call. = FALSE)
  }
  if (!inherits(fit, "sagui_sed_fit")) {
    stop("`fit` must be returned by `fit_region_seds()`.", call. = FALSE)
  }

  fit$summary <- add_property_aliases(fit$summary)
  available <- intersect(properties, names(fit$summary))
  missing <- setdiff(properties, available)
  if (length(missing)) {
    warning("Skipping missing fit properties: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  region <- as.character(fit$summary$region)
  out <- lapply(available, function(property) {
    values <- fit$summary[[property]]
    names(values) <- region
    mat <- matrix(NA_real_, nrow(labels), ncol(labels))
    ok <- !is.na(labels) & labels > 0
    mat[ok] <- values[as.character(labels[ok])]
    mat
  })
  names(out) <- available
  out
}

add_property_aliases <- function(summary) {
  if (!is.data.frame(summary)) return(summary)
  if ("log_massformed" %in% names(summary) && !"logMformed" %in% names(summary)) {
    summary$logMformed <- summary$log_massformed
  }
  has_exp_sfh <- all(c("logMformed", "age_gyr", "tau_gyr") %in% names(summary))
  if (has_exp_sfh && !"SFR" %in% names(summary)) {
    massformed <- 10^summary$logMformed
    age_yr <- summary$age_gyr * 1e9
    tau_yr <- summary$tau_gyr * 1e9
    denom <- tau_yr * (1 - exp(-age_yr / tau_yr))
    sfr <- massformed * exp(-age_yr / tau_yr) / denom
    sfr[!is.finite(sfr) | sfr < 0] <- NA_real_
    summary$SFR <- sfr
  }
  if ("SFR" %in% names(summary) && "logMformed" %in% names(summary) && !"sSFRformed" %in% names(summary)) {
    summary$sSFRformed <- summary$SFR / 10^summary$logMformed
    summary$sSFRformed[!is.finite(summary$sSFRformed) | summary$sSFRformed < 0] <- NA_real_
  }
  if ("SFR" %in% names(summary) && "logMstar" %in% names(summary) && !"sSFR" %in% names(summary)) {
    summary$sSFR <- summary$SFR / 10^summary$logMstar
    summary$sSFR[!is.finite(summary$sSFR) | summary$sSFR < 0] <- NA_real_
  }
  if ("SFR" %in% names(summary) && !"logSFR" %in% names(summary)) {
    summary$logSFR <- log10(summary$SFR)
  }
  if ("sSFRformed" %in% names(summary) && !"logsSFRformed" %in% names(summary)) {
    summary$logsSFRformed <- log10(summary$sSFRformed)
  }
  if ("sSFR" %in% names(summary) && !"logsSFR" %in% names(summary)) {
    summary$logsSFR <- log10(summary$sSFR)
  }
  summary
}

#' Write painted property maps to FITS files
#'
#' @param maps Named list of matrices returned by [saguiSED::paint_sed_properties()].
#' @param out_dir Output directory.
#' @param prefix File prefix.
#' @return Invisibly returns written file paths.
#' @export
write_property_maps <- function(maps,
                                out_dir = "sagui_sedfit/property_maps",
                                prefix = "sagui_sed") {
  if (!requireNamespace("FITSio", quietly = TRUE)) {
    stop("Install FITSio to write property maps.", call. = FALSE)
  }
  if (!is.list(maps) || is.null(names(maps))) {
    stop("`maps` must be a named list of matrices.", call. = FALSE)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- character(length(maps))
  for (i in seq_along(maps)) {
    property <- names(maps)[[i]]
    path <- file.path(out_dir, paste0(prefix, "_", property, ".fits"))
    FITSio::writeFITSim(maps[[i]], file = path, type = "double")
    paths[[i]] <- path
  }
  invisible(paths)
}

#' Write painted SED-fit property maps as a FITS cube
#'
#' @param maps Named list of matrices returned by [saguiSED::paint_sed_properties()].
#' @param path Output FITS cube path.
#' @param manifest_path Optional CSV manifest giving the plane/property mapping.
#' @param units Optional named character vector of property units.
#' @param properties Optional subset/order of properties to write.
#' @return Invisibly returns a list with the FITS path and manifest path.
#' @export
write_property_cube <- function(maps,
                                path = "sagui_sedfit/property_cube.fits",
                                manifest_path = sub("\\.fits$", "_planes.csv", path),
                                units = NULL,
                                properties = names(maps)) {
  if (!requireNamespace("FITSio", quietly = TRUE)) {
    stop("Install FITSio to write property cubes.", call. = FALSE)
  }
  if (!is.list(maps) || is.null(names(maps))) {
    stop("`maps` must be a named list of matrices.", call. = FALSE)
  }
  properties <- intersect(properties, names(maps))
  if (!length(properties)) {
    stop("No requested properties found in `maps`.", call. = FALSE)
  }

  dims <- lapply(maps[properties], dim)
  same_dims <- vapply(dims, function(x) identical(x, dims[[1]]), logical(1))
  if (!all(same_dims)) {
    stop("All property maps must have the same dimensions.", call. = FALSE)
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  cube <- array(NA_real_, dim = c(dims[[1]], length(properties)))
  for (i in seq_along(properties)) {
    cube[, , i] <- maps[[properties[[i]]]]
  }
  FITSio::writeFITSim(cube, file = path, type = "double")

  if (!is.null(manifest_path)) {
    dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
    unit_values <- rep(NA_character_, length(properties))
    if (!is.null(units)) {
      unit_values <- unname(units[properties])
    }
    manifest <- data.frame(
      plane = seq_along(properties),
      property = properties,
      unit = unit_values,
      stringsAsFactors = FALSE
    )
    utils::write.csv(manifest, manifest_path, row.names = FALSE)
  }

  invisible(list(path = path, manifest_path = manifest_path))
}
