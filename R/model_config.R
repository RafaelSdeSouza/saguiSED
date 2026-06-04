#' Define a Bagpipes model configuration
#'
#' @param sfh Star-formation history family. The prototype currently supports
#'   `"exponential"`.
#' @param dust Dust attenuation curve. The prototype currently supports
#'   `"Calzetti"`.
#' @param metallicity `"free"` or `"fixed"`.
#' @param fixed_metallicity Metallicity in solar units used when
#'   `metallicity = "fixed"`.
#' @param systematic_frac Fractional uncertainty added in quadrature to every
#'   photometric point.
#' @param quick Logical; use a fast optimisation fit rather than a full
#'   posterior sampler. This is intended for diagnostics and website examples.
#' @return A `sagui_bagpipes_model` object.
#' @export
bagpipes_model <- function(sfh = c("exponential"),
                           dust = c("Calzetti"),
                           metallicity = c("free", "fixed"),
                           fixed_metallicity = 1.0,
                           systematic_frac = 0.06,
                           quick = TRUE) {
  sfh <- match.arg(sfh)
  dust <- match.arg(dust)
  metallicity <- match.arg(metallicity)

  if (!is.numeric(systematic_frac) || length(systematic_frac) != 1 || systematic_frac < 0) {
    stop("`systematic_frac` must be a non-negative scalar.", call. = FALSE)
  }

  structure(
    list(
      backend = "bagpipes",
      sfh = sfh,
      dust = dust,
      metallicity = metallicity,
      fixed_metallicity = fixed_metallicity,
      systematic_frac = systematic_frac,
      quick = isTRUE(quick)
    ),
    class = "sagui_bagpipes_model"
  )
}

