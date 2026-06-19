
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

A slide file exposes several image series, each at a handful of
resolution levels. `fq_info()` lists them with pixel dimensions and
effective µm/px — the table `fq_read()` consults to pick the scan series
and a working resolution:

``` r
fq_info(vsi)
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

`fq_read()` returns an `fq_slide`: an RGB array in `[0, 1]` with its
physical scale and provenance. Printing one reports its dimensions,
resolution, and source:

``` r
slide <-
  fq_read(
    vsi,
    target_um_px = 4
  )
slide
#> <fq_slide> 8777 × 1349 px · 4.38 µm/px · Image_3470.vsi
```

A slide often carries more than one tissue section. `fq_split()` finds
them and returns an `fq_sections` collection — a list of `fq_section`s
that prints a one-line summary per section:

``` r
sections <-
  fq_split(
    slide,
    n = 2
  )
sections
#> <fq_sections> 2 section(s)
#>   <fq_section A> 2300 × 1145 px · 4.38 µm/px · 76% tissue · Image_3470.vsi
#>   <fq_section B> 2300 × 1143 px · 4.38 µm/px · 77% tissue · Image_3470.vsi
```

## Visualizing

`plot()` dispatches on each object. A slide plots as the whole scan:

``` r
plot(slide)
```

<img src="man/figures/README-plot-slide-1.png" alt="" width="100%" />

Plotting the collection lays the sections out as a contact sheet:

``` r
plot(sections)
```

<img src="man/figures/README-plot-sheet-1.png" alt="" width="100%" />

Handing the sections back to the slide draws each crop rectangle on the
parent — the quickest check that the split caught every section and
skipped streaks and debris:

``` r
plot(slide, sections = sections)
```

<img src="man/figures/README-plot-boxes-1.png" alt="" width="100%" />

A single section plots on its own, cropped tight to its tissue:

``` r
plot(sections[[1]])
```

<img src="man/figures/README-plot-section-1.png" alt="" width="100%" />

## Scoring fibrosis

Each section is scored by one or more unsupervised analyzers — color
clustering, collagen proportionate area, tissue density, and texture.
The workflow fits an analyzer once on the pooled sections, then reads
back a per-section score table and a severity map for each section. That
layer is the next piece under construction; its chunks slot in here as
the analyzers land.
