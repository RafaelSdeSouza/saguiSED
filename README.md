# saguiSED

`saguiSED` is the SED-fitting companion package for `sagui` regional
photometry. The core `sagui` package segments broadband or medium-band image
cubes and exports flux-conserving regional SEDs. `saguiSED` starts from those
regional SED tables, calls an optional SED-fitting backend, and paints fitted
physical quantities back onto the segmentation map.

The package is intentionally narrow. It does not perform segmentation and it
does not implement IFU spectral fitting. Spectroscopy-focused fitting belongs in
`capivara` or a capivara-specific extension.

## Installation

```r
install.packages("remotes")
remotes::install_github("RafaelSdeSouza/saguiSED")
library(saguiSED)
```

The current backend uses [Bagpipes](https://bagpipes.readthedocs.io/en/latest/).
Install Bagpipes separately following the official documentation.

```r
check_bagpipes()
```

## Minimal workflow

```r
suppressPackageStartupMessages({
  library(sagui)
  library(saguiSED)
})

filters <- jwst_nircam_filter_set()

sed <- as_sagui_sed_table(
  "region_seds_wide.csv",
  filter_set = filters,
  unit = "10nJy",
  redshift = 1.10,
  n_pix_col = "n_pix"
)

fit <- fit_region_seds(
  sed,
  backend = "bagpipes",
  model = bagpipes_model(
    sfh = "exponential",
    dust = "calzetti",
    metallicity = "free",
    systematic_frac = 0.10
  ),
  out_dir = "bagpipes_region_fits"
)
```

Plot the fitted regional SEDs:

```r
plot_sed_fit_mosaic(fit, normalize = "none")
```

Paint the fitted properties back onto the segmentation map:

```r
maps <- paint_sed_properties(seg$cluster_map, fit)

save_sed_property_map_pngs(
  maps,
  out_dir = "property_maps_png",
  prefix = "example"
)

write_property_cube(
  maps,
  path = "example_property_cube.fits"
)
```

Optionally compute graph-Laplacian-smoothed maps as a post-processing step:

```r
maps_smooth <- smooth_sed_property_maps(
  seg$cluster_map,
  fit,
  lambda = 3,
  adjacency = "queen"
)
```

The independent SED-fit summaries remain unchanged. The smoothed maps are
regularized map-level estimates for visualization and gradient-level
interpretation.

## Filter sets

`saguiSED` supports custom filter curves through `sed_filter_set()` and bundles
Rubin/LSST curves through `lsst_filter_set()`. JWST/NIRCam curves can be loaded
from an installed `sagui` package through `jwst_nircam_filter_set()`.

## Credits

The current backend uses Bagpipes. If you use this backend in a paper, cite
Bagpipes following the guidance in the official Bagpipes documentation:
<https://bagpipes.readthedocs.io/en/latest/>.
