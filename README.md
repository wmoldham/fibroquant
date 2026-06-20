
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
# vsi <- "/Volumes/Will/Mouse lung 6.10.26/Image_3470.vsi"
vsi <- "/run/media/will/Will/Mouse lung 6.10.26/Image_3470.vsi"
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
#>   <fq_section A> 2300 × 1145 px · 4.38 µm/px · 18% tissue · Image_3470.vsi
#>   <fq_section B> 2300 × 1143 px · 4.38 µm/px · 18% tissue · Image_3470.vsi
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

### k-means colour analyzer

The first analyzer reproduces the LungDamage colour approach —
unsupervised CIELAB clustering of tissue colour, clusters ranked
dark-to-light into three severity grades (severe / moderate / healthy) —
with two changes. It clusters only masked tissue pixels: the mask does
the job LungDamage’s brightest, background cluster did, so three
clusters suffice where MATLAB used four. And it fits one shared cluster
basis across all sections at once (fit-once / apply-many) instead of
re-clustering each image.

``` r
slide    <- fq_read(vsi)
sections <- fq_split(slide, n = 2)
sec      <- sections[[1]]
```

#### 1. The spec

`fq_kmeans()` is just the recipe — how many grades, which colour
channels, the blur, the number of k-means restarts. It holds no fitted
state.

``` r
spec <- 
  fq_kmeans(
    k = 3,
    channels = c("a", "b"),
    smooth_sigma = 2,
    nstart = 3
  )
spec
#> <fibroquant::fq_kmeans>
#>  @ k           : num 3
#>  @ channels    : chr [1:2] "a" "b"
#>  @ smooth_sigma: num 2
#>  @ nstart      : num 3
#>  @ max_px      : num 1e+05
```

#### 2. Smoothing

A Gaussian blur (sigma in pixels) damps single-pixel stain speckle
before clustering. Same H×W×3 array, slightly softened.

``` r
sm <- fibroquant:::.smooth(sec@rgb, sigma = 0)
EBImage::display(EBImage::Image(sm, colormode = "Color"), method = "raster")
```

<img src="man/figures/README-kmeans-smooth-1.png" alt="" width="100%" />

#### 3. CIELAB — a and b channels

`rgb2lab`, then keep `a*` and `b*`: the colour axes, with lightness
(`L*`) set aside. The two channels are shown below.

``` r
lab <- fibroquant:::.lab(sm)
EBImage::display(EBImage::normalize(lab[, , 2]), method = "raster") # a*
```

<img src="man/figures/README-kmeans-lab-1.png" alt="" width="100%" />

``` r
EBImage::display(EBImage::normalize(lab[, , 3]), method = "raster") # b*
```

<img src="man/figures/README-kmeans-lab-2.png" alt="" width="100%" />

#### 4. The clustering input

Pull the `a*`/`b*` values at masked tissue pixels only — the clustering
input, and the fibroquant change from LungDamage. An N×2 matrix, one row
per tissue pixel.

``` r
feat <- fibroquant:::.features(lab, sec@mask, channels = c("a", "b"))
dim(feat)
#> [1] 463327      2
head(feat)
#>              a         b
#> [1,]  2.808943  -1.24876
#> [2,] 32.243008 -21.36482
#> [3,] 32.228354 -22.48672
#> [4,] 28.242389 -21.39894
#> [5,] 22.065641 -15.48909
#> [6,] 26.420241 -20.57937
```

#### 5. Fitting the basis

`fq_fit` pools tissue pixels across every section, runs k-means once,
then ranks the clusters by ascending tissue luminance (`L*`) — darkest
grade = most fibrotic. The fit carries the cluster centres and that
severity ordering; nothing is stored back on the sections.

``` r
fit <- fq_fit(spec, sections)
fit@centers  # 3 x 2 in a*b* space
#>          a         b
#> 1 23.65692 -16.48945
#> 3 31.97573 -21.89368
#> 2 48.22035 -26.16158
fit@luminance # cluster -> rank: mean L*, severity 1..3
#> [1] 69.98659 66.82713 64.28321
```

#### 6. Per-section score

`fq_score` applies the fitted basis to one section and measures area per
grade. One row: `severity_index` (area-weighted, 0–10) plus `frac_sev_1`
… `frac_sev_3`.

``` r
fq_score(fit, sec)
#> # A tibble: 1 × 4
#>   severity_index frac_sev_1 frac_sev_2 frac_sev_3
#>            <dbl>      <dbl>      <dbl>      <dbl>
#> 1           3.81      0.418      0.402      0.180
```

#### 7. The severity map

`fq_render` labels each tissue pixel with its grade, returning an
`fq_render` object (an H×W field, NA off tissue). Plotting it
pseudocolours the grades (severe / moderate / healthy → R/G/B); that
plot method lands with the plot functions.

``` r
rendered <- fq_render(fit, sec)
plot(rendered)
```

<img src="man/figures/README-kmeans-render-1.png" alt="" width="100%" />

## Batch processing

Scoring generalises from one slide to a folder of them. `fq_manifest()`
discovers the slides and gives each a stable `slide_id` for joining
scores back to files and naming saved maps:

``` r
manifest <- fq_manifest(dirname(vsi))
manifest
#> # A tibble: 173 × 2
#>    slide_id      path                                                     
#>    <chr>         <chr>                                                    
#>  1 Image_3470_01 /run/media/will/Will/Mouse lung 6.10.26/Image_3470_01.vsi
#>  2 Image_3470    /run/media/will/Will/Mouse lung 6.10.26/Image_3470.vsi   
#>  3 Image_3472_01 /run/media/will/Will/Mouse lung 6.10.26/Image_3472_01.vsi
#>  4 Image_3472    /run/media/will/Will/Mouse lung 6.10.26/Image_3472.vsi   
#>  5 Image_3473_01 /run/media/will/Will/Mouse lung 6.10.26/Image_3473_01.vsi
#>  6 Image_3473    /run/media/will/Will/Mouse lung 6.10.26/Image_3473.vsi   
#>  7 Image_3474_01 /run/media/will/Will/Mouse lung 6.10.26/Image_3474_01.vsi
#>  8 Image_3474    /run/media/will/Will/Mouse lung 6.10.26/Image_3474.vsi   
#>  9 Image_3475_01 /run/media/will/Will/Mouse lung 6.10.26/Image_3475_01.vsi
#> 10 Image_3475    /run/media/will/Will/Mouse lung 6.10.26/Image_3475.vsi   
#> # ℹ 163 more rows
```

For a batch the cluster basis is fit **once** and then frozen, so a
grade means the same thing on every slide. `fq_read_sections()` reads
and splits a chosen subsample into one pooled list, which `fq_fit()`
clusters into a single batch-wide severity scale. Pass a subsample —
e.g. `dplyr::slice_sample(manifest, n = 8)` — for large batches:

``` r
manifest <- 
  fq_manifest(dirname(vsi)) |> 
  dplyr::filter(!stringr::str_detect(slide_id, "_01")) |> 
  dplyr::slice_head(n = 10)
manifest$treatment <- "bleo"   # attach any covariates as columns

result <- fq_run(manifest, fq_kmeans(k = 3))
result$scores
#> # A tibble: 20 × 8
#>    slide_id   section n_sections severity_index frac_sev_1 frac_sev_2 frac_sev_3
#>    <chr>      <chr>        <int>          <dbl>      <dbl>      <dbl>      <dbl>
#>  1 Image_3470 A                2           3.94      0.404      0.403      0.193
#>  2 Image_3470 B                2           4.29      0.346      0.449      0.205
#>  3 Image_3472 A                2           2.65      0.589      0.292      0.119
#>  4 Image_3472 B                2           3.02      0.512      0.371      0.117
#>  5 Image_3473 A                2           3.55      0.468      0.354      0.178
#>  6 Image_3473 B                2           3.45      0.492      0.327      0.181
#>  7 Image_3474 A                2           3.99      0.378      0.447      0.175
#>  8 Image_3474 B                2           4.51      0.290      0.519      0.192
#>  9 Image_3475 A                2           3.64      0.519      0.234      0.246
#> 10 Image_3475 B                2           3.82      0.486      0.264      0.250
#> 11 Image_3477 A                2           5.31      0.256      0.427      0.317
#> 12 Image_3477 B                2           5.43      0.238      0.437      0.324
#> 13 Image_3478 A                2           4.61      0.344      0.390      0.266
#> 14 Image_3478 B                2           4.53      0.349      0.396      0.255
#> 15 Image_3479 A                2           5.39      0.219      0.484      0.297
#> 16 Image_3479 B                2           5.43      0.214      0.486      0.300
#> 17 Image_3480 A                2           5.09      0.286      0.411      0.303
#> 18 Image_3480 B                2           5.33      0.264      0.405      0.331
#> 19 Image_3481 A                2           6.09      0.179      0.426      0.396
#> 20 Image_3481 B                2           6.23      0.165      0.423      0.412
#> # ℹ 1 more variable: treatment <chr>
```

That frozen `batch_fit` scores and renders every section across the
folder on a shared scale. The driver that maps it over the manifest and
aggregates the per-section rows into one results tibble — in parallel
via `future`/`furrr` — is the next piece.
