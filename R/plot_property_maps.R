.map_to_df <- function(mat) {
  if (!is.matrix(mat)) stop("Property map must be a matrix.", call. = FALSE)
  grid <- expand.grid(
    row = seq_len(nrow(mat)),
    col = seq_len(ncol(mat))
  )
  grid$value <- as.vector(mat)
  # Match sagui::plot_region_map(): x is matrix row, y is matrix column.
  # This avoids a silent transpose relative to the segmentation map.
  grid$x <- grid$row
  grid$y <- grid$col
  grid
}

#' Plot painted SED-fit property maps
#'
#' @param maps Named list of matrices returned by [paint_sed_properties()].
#' @param properties Optional subset of properties to plot.
#' @param palette Continuous color palette.
#' @param background_color Color used for non-finite pixels.
#' @param base_size Base font size.
#' @return A named list of `ggplot2` objects.
#' @export
plot_sed_property_maps <- function(maps,
                                   properties = names(maps),
                                   palette = c("#213E60", "#94B6EF", "#F4F2EF", "#E68C3A"),
                                   background_color = "white",
                                   base_size = 13) {
  if (!is.list(maps) || is.null(names(maps))) {
    stop("`maps` must be a named list of matrices.", call. = FALSE)
  }
  properties <- intersect(properties, names(maps))
  if (!length(properties)) {
    stop("No requested properties found in `maps`.", call. = FALSE)
  }

  labels <- c(
    logMformed = expression(log[10]~M["formed"]~"[M"["sun"]*"]"),
    log_massformed = expression(log[10]~M["formed"]~"[M"["sun"]*"]"),
    SFR = expression(SFR~"[M"["sun"]~yr^{-1}*"]"),
    logSFR = expression(log[10]~SFR~"[M"["sun"]~yr^{-1}*"]"),
    sSFRformed = expression(SFR/M["formed"]~"[yr"^{-1}*"]"),
    logsSFRformed = expression(log[10]~"(SFR/M"["formed"]*")"),
    sSFR = expression(sSFR~"[yr"^{-1}*"]"),
    logsSFR = expression(log[10]~sSFR~"[yr"^{-1}*"]"),
    age_gyr = expression("<t>"[M]~"[Gyr]"),
    tau_gyr = expression(tau~"[Gyr]"),
    logZ_Zsun = expression(log[10]~"(Z/Z"["sun"]*")"),
    Av = expression(A[V]~"[mag]"),
    fit_rms_sigma = expression(RMS~"["*sigma*"]")
  )

  out <- lapply(properties, function(property) {
    df <- .map_to_df(maps[[property]])
    legend_label <- if (property %in% names(labels)) labels[[property]] else property
    finite_values <- df$value[is.finite(df$value)]
    breaks <- if (length(finite_values) && diff(range(finite_values)) > 0) {
      seq(min(finite_values), max(finite_values), length.out = 3)
    } else {
      ggplot2::waiver()
    }
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::geom_raster(na.rm = FALSE) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::scale_fill_gradientn(
        colours = palette,
        na.value = background_color,
        name = legend_label,
        breaks = breaks,
        labels = function(x) format(signif(x, 3), trim = TRUE, scientific = FALSE),
        guide = ggplot2::guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          barwidth = grid::unit(4.6, "cm"),
          barheight = grid::unit(0.2, "cm")
        )
      ) +
      ggplot2::coord_equal() +
      ggplot2::theme_void(base_size = base_size) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.title = ggplot2::element_text(color = "#101820", size = base_size * 0.82),
        legend.text = ggplot2::element_text(color = "#101820", size = base_size * 0.72),
        legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
        legend.box.margin = ggplot2::margin(t = -4, r = 0, b = 0, l = 0),
        plot.margin = ggplot2::margin(t = 4, r = 18, b = 8, l = 18),
        plot.background = ggplot2::element_rect(fill = background_color, color = NA),
        panel.background = ggplot2::element_rect(fill = background_color, color = NA)
      )
  })
  names(out) <- properties
  out
}

#' Save painted SED-fit property maps as PNG files
#'
#' @param maps Named list of matrices returned by [paint_sed_properties()].
#' @param out_dir Output directory.
#' @param prefix File prefix.
#' @param properties Optional subset of properties to save.
#' @param width,height Plot dimensions in inches.
#' @param dpi Output resolution.
#' @param ... Additional arguments passed to [plot_sed_property_maps()].
#' @return Invisibly returns written file paths.
#' @export
save_sed_property_map_pngs <- function(maps,
                                       out_dir = "sagui_sedfit/property_maps_png",
                                       prefix = "sagui_sed",
                                       properties = names(maps),
                                       width = 5.2,
                                       height = 5.2,
                                       dpi = 220,
                                       ...) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plots <- plot_sed_property_maps(maps, properties = properties, ...)
  paths <- character(length(plots))
  for (i in seq_along(plots)) {
    property <- names(plots)[[i]]
    path <- file.path(out_dir, paste0(prefix, "_", property, ".png"))
    ggplot2::ggsave(path, plots[[i]], width = width, height = height, dpi = dpi)
    paths[[i]] <- path
  }
  invisible(paths)
}
