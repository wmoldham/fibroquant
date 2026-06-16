
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fibroquant

`fibroquant` quantifies lung fibrosis from Olympus `.vsi` whole-slide
images (bleomycin mouse model; Masson trichrome and H&E) with no manual
annotation. It reads a scan, splits the two lung sections side by side
on each slide, and scores each section by unsupervised color clustering.

## Reading and splitting a slide

``` r
# library(fibroquant)
devtools::load_all()
#> ℹ Loading fibroquant

# vsi <- "/path/to/vsi/image"
vsi <- "/Volumes/Will/Mouse lung 6.10.26/Image_3470.vsi"
has_vsi <- file.exists(vsi)
```

The series-and-resolution table the file exposes:

``` r
fq_slide_info(vsi)
```

Read the scan series at a working resolution near 4 µm/px:

``` r
slide <-
  fq_read_slide(
    vsi,
    target_um_px = 4
  )
slide
```

Split the two lung sections:

``` r
sections <-
  fq_split_sections(
    slide,
    n = 2
  )
```

Visualization and the k-means analyzer slot in below as we build them.
