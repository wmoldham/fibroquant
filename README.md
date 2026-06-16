
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fibroquant

`fibroquant` quantifies lung fibrosis from histology whole-slide images
without manual annotation. It reads a slide, optionally splits a
multi-section slide into its separate sections, and scores tissue with a
choice of unsupervised analyzers: color clustering, collagen
proportionate area, tissue density, and texture.

## Reading and splitting a slide

``` r
# library(fibroquant)
devtools::load_all()
#> ℹ Loading fibroquant

# vsi <- "/path/to/vsi/image"
vsi <- "/Volumes/Will/Mouse lung 6.10.26/Image_3470.vsi"
has_vsi <- file.exists(vsi)
```

List the file’s image series and resolution levels with their pixel
dimensions and effective µm/px — the table `fq_read_slide()` reads to
pick the scan series and a working resolution:

``` r
fq_slide_info(vsi)
#> # A tibble: 24 × 5
#>    series   res size_x size_y um_px
#>     <int> <int>  <int>  <int> <dbl>
#>  1      1     1   8021   9366 0.274
#>  2      1     2   4011   4683 0.548
#>  3      1     3   2006   2342 1.09 
#>  4      1     4   1003   1171 2.19 
#>  5      1     5    502    586 4.37 
#>  6      1     6    251    293 8.75 
#>  7      2     1  18032   9148 0.274
#>  8      2     2   9016   4574 0.548
#>  9      2     3   4508   2287 1.10 
#> 10      2     4   2254   1144 2.19 
#> # ℹ 14 more rows
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
