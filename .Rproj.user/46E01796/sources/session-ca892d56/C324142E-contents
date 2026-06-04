.sagui_filter_colors <- function(wavelength_um) {
  stopifnot(is.numeric(wavelength_um))
  anchors <- c(0.35, 0.48, 0.62, 0.76, 0.87, 0.98, 1.2, 1.6, 2.0, 3.0, 4.5, 5.0)
  colors <- c(
    "#5A45B8",
    "#2F80ED",
    "#36A96D",
    "#E6B84A",
    "#E68C3A",
    "#C44E3B",
    "#9D3A2D",
    "#20B2AA",
    "#8FD744",
    "#F4A64A",
    "#E56B3A",
    "#9D3A2D"
  )
  ramp <- grDevices::colorRamp(colors, space = "Lab")
  x <- pmin(pmax((wavelength_um - min(anchors)) / diff(range(anchors)), 0), 1)
  rgb <- ramp(x)
  grDevices::rgb(rgb[, 1], rgb[, 2], rgb[, 3], maxColorValue = 255)
}

.fit_photometry_for_plot <- function(fit,
                                     regions = NULL,
                                     normalize = c("none", "median")) {
  normalize <- match.arg(normalize)
  if (!inherits(fit, "sagui_sed_fit")) {
    stop("`fit` must be returned by `fit_region_seds()`.", call. = FALSE)
  }
  phot <- fit$model_photometry
  if (!is.data.frame(phot) || !nrow(phot)) {
    stop("`fit$model_photometry` is empty.", call. = FALSE)
  }
  required <- c("region", "filter", "wave_um", "obs_ujy", "err_ujy", "model_ujy")
  missing <- setdiff(required, names(phot))
  if (length(missing)) {
    stop("Missing model-photometry columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!is.null(regions)) {
    phot <- phot[phot$region %in% regions, , drop = FALSE]
  }
  if (!nrow(phot)) {
    stop("No photometry rows remain after applying `regions`.", call. = FALSE)
  }

  phot$region <- as.integer(phot$region)
  phot <- phot[order(phot$region, phot$wave_um), , drop = FALSE]
  scale_tbl <- stats::aggregate(obs_ujy ~ region, data = phot, FUN = function(x) {
    s <- stats::median(x[is.finite(x) & x > 0], na.rm = TRUE)
    if (!is.finite(s) || s <= 0) 1 else s
  })
  names(scale_tbl)[names(scale_tbl) == "obs_ujy"] <- "scale"
  if (identical(normalize, "median")) {
    scale <- scale_tbl$scale[match(phot$region, scale_tbl$region)]
    phot$obs_plot <- phot$obs_ujy / scale
    phot$err_plot <- phot$err_ujy / scale
    phot$model_plot <- phot$model_ujy / scale
    y_label <- "Relative flux"
  } else {
    scale_tbl$scale <- 1
    phot$obs_plot <- phot$obs_ujy
    phot$err_plot <- phot$err_ujy
    phot$model_plot <- phot$model_ujy
    y_label <- expression(F[nu]~"["*mu*"Jy]")
  }
  phot$region_label <- paste0("Region ", phot$region)
  list(data = phot, scale = scale_tbl, y_label = y_label)
}

.fit_spectrum_for_plot <- function(fit, prep) {
  spectrum <- fit$model_spectrum
  if (!is.data.frame(spectrum) || !nrow(spectrum)) {
    return(NULL)
  }
  required <- c("region", "wave_um", "model_ujy")
  if (!all(required %in% names(spectrum))) {
    return(NULL)
  }
  regions <- unique(prep$data$region)
  spectrum <- spectrum[spectrum$region %in% regions, , drop = FALSE]
  if (!nrow(spectrum)) return(NULL)
  spectrum$region <- as.integer(spectrum$region)
  spectrum <- merge(spectrum, prep$scale, by = "region", all.x = TRUE, sort = FALSE)
  spectrum$scale[!is.finite(spectrum$scale) | spectrum$scale <= 0] <- 1
  spectrum$model_plot <- spectrum$model_ujy / spectrum$scale
  spectrum$region_label <- paste0("Region ", spectrum$region)
  spectrum <- spectrum[is.finite(spectrum$wave_um) & is.finite(spectrum$model_plot), , drop = FALSE]
  spectrum[order(spectrum$region, spectrum$wave_um), , drop = FALSE]
}

.transmission_overlay <- function(phot,
                                  filter_set,
                                  transmission_height = 0.45) {
  if (!inherits(filter_set, "sagui_filter_set") || transmission_height <= 0) {
    return(NULL)
  }
  trans <- filter_set$throughput
  if (!all(c("filter", "wavelength_um", "throughput") %in% names(trans))) {
    return(NULL)
  }
  trans <- trans[trans$filter %in% unique(phot$filter), , drop = FALSE]
  if (!nrow(trans)) return(NULL)

  ranges <- stats::aggregate(
    cbind(obs_plot, model_plot) ~ region + region_label,
    data = phot,
    FUN = function(x) max(x[is.finite(x)], na.rm = TRUE)
  )
  ranges$ymax <- pmax(ranges$obs_plot, ranges$model_plot)
  ranges$ymax[!is.finite(ranges$ymax) | ranges$ymax <= 0] <- 1
  ranges$ybase <- 0.03 * ranges$ymax
  ranges$yheight <- transmission_height * ranges$ymax
  ranges <- ranges[, c("region", "region_label", "ymax", "ybase", "yheight"), drop = FALSE]

  out <- merge(trans, ranges, by = NULL)
  max_t <- stats::ave(out$throughput, out$filter, FUN = function(x) max(x, na.rm = TRUE))
  max_t[!is.finite(max_t) | max_t <= 0] <- 1
  out$trans_ymin <- out$ybase
  out$trans_ymax <- out$ybase + (out$throughput / max_t) * out$yheight
  out <- out[order(out$region, out$filter, out$wavelength_um), , drop = FALSE]
  out
}

#' Plot a mosaic of regional SED fits
#'
#' @param fit A `sagui_sed_fit` object returned by [fit_region_seds()].
#' @param regions Optional region IDs to show.
#' @param filter_set Optional `sagui_filter_set`. If omitted, the filter set
#'   stored in `fit` is used.
#' @param normalize Plot absolute fluxes (`"none"`) or divide each region by
#'   its median observed flux (`"median"`).
#' @param ncol Number of facet columns.
#' @param transmission_height Fractional height used for filter-transmission
#'   curves. Set to `0` to hide them.
#' @param point_size Observed-photometry point size.
#' @param errorbar_width Error-bar width in microns.
#' @param base_size Base font size.
#' @return A `ggplot2` object.
#' @export
plot_sed_fit_mosaic <- function(fit,
                                regions = NULL,
                                filter_set = NULL,
                                normalize = c("none", "median"),
                                ncol = NULL,
                                transmission_height = 0.20,
                                point_size = 2.8,
                                errorbar_width = 0.035,
                                base_size = 12) {
  normalize <- match.arg(normalize)
  prep <- .fit_photometry_for_plot(fit, regions = regions, normalize = normalize)
  phot <- prep$data
  spectrum <- .fit_spectrum_for_plot(fit, prep)
  if (is.null(filter_set)) {
    filter_set <- fit$filter_set
  }

  lambda_by_filter <- stats::aggregate(wave_um ~ filter, data = phot, FUN = stats::median)
  filter_cols <- .sagui_filter_colors(lambda_by_filter$wave_um)
  names(filter_cols) <- lambda_by_filter$filter

  trans <- .transmission_overlay(
    phot,
    filter_set = filter_set,
    transmission_height = transmission_height
  )
  labels <- unique(phot[, c("region", "region_label"), drop = FALSE])
  if ("fit_rms_sigma" %in% names(fit$summary)) {
    labels$rms <- fit$summary$fit_rms_sigma[match(labels$region, fit$summary$region)]
    labels$label <- ifelse(
      is.finite(labels$rms),
      paste0("region ", labels$region, "\nRMS = ", sprintf("%.2f", labels$rms), " sigma"),
      paste0("region ", labels$region)
    )
  } else {
    labels$label <- paste0("region ", labels$region)
  }
  xlim <- c(
    0.92 * min(phot$wave_um, na.rm = TRUE),
    1.07 * max(phot$wave_um, na.rm = TRUE)
  )

  p <- ggplot2::ggplot(phot, ggplot2::aes(x = wave_um))
  if (!is.null(trans)) {
    p <- p +
      ggplot2::geom_ribbon(
        data = trans,
        ggplot2::aes(
          x = wavelength_um,
          ymin = trans_ymin,
          ymax = trans_ymax,
          group = filter,
          fill = filter
        ),
        inherit.aes = FALSE,
        alpha = 0.23,
        color = NA
      )
  }

  if (!is.null(spectrum)) {
    p <- p +
      ggplot2::geom_line(
        data = spectrum,
        ggplot2::aes(x = wave_um, y = model_plot, group = region),
        inherit.aes = FALSE,
        color = "#213E60",
        linewidth = 1.05,
        lineend = "round"
      )
  }

  p +
    ggplot2::geom_line(
      ggplot2::aes(y = model_plot),
      color = "#213E60",
      linewidth = 0.55,
      alpha = 0.24,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = model_plot),
      color = "#213E60",
      size = point_size * 0.42,
      alpha = 0.75
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(obs_plot - err_plot, 0),
        ymax = obs_plot + err_plot,
        color = filter
      ),
      width = errorbar_width,
      linewidth = 0.62,
      alpha = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = obs_plot, color = filter),
      shape = 21,
      fill = "white",
      stroke = 1.05,
      size = point_size
    ) +
    ggplot2::geom_text(
      data = labels,
      ggplot2::aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.08,
      vjust = 1.18,
      color = "#213E60",
      fontface = "bold",
      size = base_size / 3.25,
      lineheight = 0.88
    ) +
    ggplot2::facet_wrap(ggplot2::vars(region_label), ncol = ncol, scales = "free_y") +
    ggplot2::scale_color_manual(values = filter_cols, guide = "none") +
    ggplot2::scale_fill_manual(values = filter_cols, guide = "none") +
    ggplot2::coord_cartesian(xlim = xlim, clip = "on") +
    ggplot2::labs(x = "Observed wavelength [micron]", y = prep$y_label) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(color = "#101820"),
      axis.text = ggplot2::element_text(color = "#101820"),
      panel.border = ggplot2::element_rect(fill = NA, color = "#222222", linewidth = 0.45),
      panel.spacing = grid::unit(0.38, "lines"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}
