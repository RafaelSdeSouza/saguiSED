#' Smooth SED-fit property maps on the SAGUI region graph
#'
#' This is a post-processing step. It does not refit the photometry and does not
#' alter the independent SED-fit summaries returned by [fit_region_seds()].
#' Instead, it regularizes region-level physical quantities over the adjacency
#' graph of the segmentation using `sagui::smooth_region_field_laplacian()`.
#'
#' @param labels Region-label matrix or a SAGUI segmentation object containing
#'   `cluster_map`.
#' @param fit A `sagui_sed_fit` object returned by [fit_region_seds()].
#' @param properties Fit-summary columns to smooth and paint.
#' @param lambda Smoothing strength passed to
#'   `sagui::smooth_region_field_laplacian()`.
#' @param adjacency Neighborhood definition, either `"queen"` or `"rook"`.
#' @param unsmoothed_properties Properties to paint directly without smoothing.
#'   Diagnostics such as `fit_rms_sigma` should usually remain unsmoothed.
#' @param keep_na_outside Keep non-region pixels as `NA`.
#' @return A named list of matrices.
#' @export
smooth_sed_property_maps <- function(labels,
                                     fit,
                                     properties = c(
                                       "logMformed",
                                       "SFR",
                                       "logSFR",
                                       "logsSFRformed",
                                       "age_gyr",
                                       "logZ_Zsun",
                                       "Av",
                                       "fit_rms_sigma"
                                     ),
                                     lambda = 3,
                                     adjacency = c("queen", "rook"),
                                     unsmoothed_properties = "fit_rms_sigma",
                                     keep_na_outside = TRUE) {
  adjacency <- match.arg(adjacency)

  if (is.list(labels) && !is.null(labels$cluster_map)) {
    labels <- labels$cluster_map
  }
  if (!is.matrix(labels)) {
    stop("`labels` must be a matrix or a SAGUI segmentation object.", call. = FALSE)
  }
  if (!inherits(fit, "sagui_sed_fit")) {
    stop("`fit` must be returned by `fit_region_seds()`.", call. = FALSE)
  }
  if (!requireNamespace("sagui", quietly = TRUE)) {
    stop("Install `sagui` to use graph-Laplacian property smoothing.", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("`lambda` must be a finite non-negative scalar.", call. = FALSE)
  }

  summary <- add_property_aliases(fit$summary)
  properties <- intersect(properties, names(summary))
  if (!length(properties)) {
    stop("No requested properties found in `fit$summary`.", call. = FALSE)
  }

  raw_maps <- paint_sed_properties(labels, fit, properties = properties)
  region <- as.integer(summary$region)
  out <- lapply(properties, function(property) {
    if (property %in% unsmoothed_properties) {
      return(raw_maps[[property]])
    }

    values <- as.numeric(summary[[property]])
    names(values) <- region
    values <- values[is.finite(values)]
    if (length(values) < 2L) {
      return(raw_maps[[property]])
    }

    sagui::smooth_region_field_laplacian(
      region_id_mat = labels,
      region_values = values,
      adjacency = adjacency,
      lambda = lambda,
      keep_na_outside = keep_na_outside
    )$interpolated_matrix
  })
  names(out) <- properties
  attr(out, "lambda") <- lambda
  attr(out, "adjacency") <- adjacency
  attr(out, "unsmoothed_properties") <- unsmoothed_properties
  class(out) <- c("sagui_sed_property_maps", "list")
  out
}
