# kmeans.R

#' k-means colour-segmentation spec
#'
#' [spec][fq_spec] for the k-means analyzer: cluster masked tissue pixels in
#' CIELAB colour space, then rank the clusters by luminance into severity
#' grades. Pass the result to [fq_fit()].
#'
#' @param k Number of severity grades (clusters). Default 3 reproduces the three
#'   tissue grades of the MATLAB LungDamage pipeline, once the tissue mask
#'   excludes background.
#' @param channels CIELAB channels to cluster on, a subset of `"L"`, `"a"`,
#'   `"b"`. Default `c("a", "b")` clusters on chroma alone, leaving lightness to
#'   the severity ordering.
#' @param smooth_sigma Gaussian-blur sigma in pixels applied before clustering,
#'   to damp single-pixel stain speckle. The blur is a mask-aware normalised
#'   convolution, so airspace carries no weight into tissue and cannot bleed
#'   into the rim. `0` disables smoothing.
#' @param nstart Number of k-means restarts; the lowest within-cluster
#'   sum-of-squares fit is kept.
#' @param max_px Maximum tissue pixels sampled per section when fitting. Caps
#'   each section's contribution so slides are weighted equally and the fit
#'   stays tractable on large batches; `Inf` uses every pixel. Sampling is
#'   random, so set a seed for a reproducible fit.
#' @return An `fq_kmeans` spec.
#' @export
fq_kmeans <-
  S7::new_class(
    "fq_kmeans",
    parent = fq_spec,
    properties = list(
      k = S7::new_property(S7::class_numeric, default = 3),
      channels = S7::new_property(S7::class_character, default = c("a", "b")),
      smooth_sigma = S7::new_property(S7::class_numeric, default = 2),
      nstart = S7::new_property(S7::class_numeric, default = 3),
      max_px = S7::new_property(S7::class_numeric, default = 1e5)
    ),
    validator = function(self) {
      if (length(self@k) != 1L || is.na(self@k) ||
          self@k < 2 || self@k %% 1 != 0) {
        return("`k` must be a single whole number >= 2")
      }
      if (length(self@channels) < 1L ||
          !all(self@channels %in% c("L", "a", "b"))) {
        return("`channels` must be a non-empty subset of 'L', 'a', 'b'")
      }
      if (length(self@smooth_sigma) != 1L || is.na(self@smooth_sigma) ||
          self@smooth_sigma < 0) {
        return("`smooth_sigma` must be a single number >= 0")
      }
      if (length(self@nstart) != 1L || is.na(self@nstart) ||
          self@nstart < 1 || self@nstart %% 1 != 0) {
        return("`nstart` must be a single whole number >= 1")
      }
      if (length(self@max_px) != 1L || is.na(self@max_px) ||
          self@max_px < 1) {
        return("`max_px` must be a single number >= 1 (Inf disables the cap)")
      }
      NULL
    }
  )

#' Fitted k-means analyzer
#'
#' The basis learned by [fq_fit()] from an [fq_kmeans()] spec: cluster centres
#' reordered into severity grades (row 1 = mildest, row k = most severe) and
#' each grade's mean lightness. Consumed by [fq_score()] and [fq_render()].
#'
#' @param spec The [fq_kmeans()] spec this basis was fit from.
#' @param centers Cluster centres in the clustered channels, one row per grade,
#'   ordered mildest (row 1) to most severe (row k).
#' @param luminance Mean L* per grade, in the same order as the rows of
#'   `centers`.
#' @return An `fq_kmeans_analyzer` object.
#' @export
fq_kmeans_analyzer <-
  S7::new_class(
    "fq_kmeans_analyzer",
    parent = fq_analyzer,
    properties = list(
      spec = fq_kmeans,
      centers = S7::class_numeric,
      luminance = S7::class_numeric
    )
  )

#' Per-pixel severity field
#'
#' The spatial output of [fq_render()]: an H x W grid carrying each tissue
#' pixel's severity grade (1 to `k`), `NA` off tissue. The plot functions
#' pseudocolour it.
#'
#' @param values Numeric H x W grid of severity grades (1 to `k`), `NA` off
#'   tissue.
#' @param k Number of severity grades.
#' @return An `fq_field` object.
#' @export
fq_field <-
  S7::new_class(
    "fq_field",
    properties = list(
      values = S7::class_numeric,
      k = S7::class_integer
    )
  )

# Gaussian-blur each RGB channel to damp single-pixel stain speckle. sigma is in
# pixels; 0 leaves the image unchanged. With a mask, blur as a normalised
# convolution so off-mask pixels (airspace) carry no weight into masked tissue.
.smooth <- function(rgb, sigma, mask = NULL) {
  if (sigma <= 0) {
    return(rgb)
  }
  if (is.null(mask)) {
    return(
      EBImage::gblur(
        rgb,
        sigma = sigma
      )
    )
  }

  m <- mask * 1
  weight <- EBImage::gblur(
    m,
    sigma = sigma
  )
  out <- rgb
  for (ch in seq_len(dim(rgb)[3L])) {
    blurred <- EBImage::gblur(
      rgb[, , ch] * m,
      sigma = sigma
    )
    smoothed <- blurred / weight
    smoothed[!mask] <- rgb[, , ch][!mask]  # off-tissue is unused; avoid 0/0
    out[, , ch] <- smoothed
  }
  out
}

# Convert sRGB to CIELAB pixelwise; returns an H x W x 3 array of L*, a*, b*.
.lab <- function(rgb) {
  flat <- matrix(rgb, ncol = 3)
  lab <- grDevices::convertColor(flat, from = "sRGB", to = "Lab")
  array(
    lab,
    dim = dim(rgb)
  )
}

# Select the chosen CIELAB channels at masked tissue pixels into an N x m
# matrix: rows are tissue pixels, columns are the requested channels.
.features <- function(lab, mask, channels) {
  idx <- match(channels, c("L", "a", "b"))
  flat <- matrix(lab[, , idx], ncol = length(idx))
  out <- flat[as.vector(mask), , drop = FALSE]
  colnames(out) <- channels
  out
}

# Severity grade of each masked tissue pixel: the nearest centre row, which is
# already in severity order. Returns a grade per pixel in column-major order.
.assign <- function(section, fit) {
  spec <- fit@spec
  lab <- .lab(.smooth(section@rgb, spec@smooth_sigma, section@mask))
  x <- .features(lab, section@mask, spec@channels)
  centers <- fit@centers
  d <- vapply(
    seq_len(nrow(centers)),
    function(g) rowSums(sweep(x, 2, centers[g, ])^2),
    numeric(nrow(x))
  )
  max.col(-d, ties.method = "first")
}

# Pool masked tissue pixels across sections, k-means on the chosen channels,
# then order clusters by descending mean L* so centre row i is severity grade i.
# Each section is capped to spec@max_px pixels before pooling, so a large
# section cannot outweigh a small one and the fit stays tractable.
S7::method(fq_fit, fq_kmeans) <- function(spec, sections, ...) {
  prepped <- lapply(sections, function(s) {
    lab <- .lab(.smooth(s@rgb, spec@smooth_sigma, s@mask))
    x <- .features(lab, s@mask, spec@channels)
    lum <- .features(lab, s@mask, "L")[, 1]
    if (nrow(x) > spec@max_px) {
      keep <- sample.int(nrow(x), as.integer(spec@max_px))
      x <- x[keep, , drop = FALSE]
      lum <- lum[keep]
    }
    list(x = x, lum = lum)
  })
  x <- do.call(rbind, lapply(prepped, function(p) p$x))
  lum <- unlist(lapply(prepped, function(p) p$lum))

  k <- as.integer(spec@k)
  if (nrow(x) < k) {
    stop(
      sprintf("Too few tissue pixels (%d) to fit %d clusters.", nrow(x), k),
      call. = FALSE
    )
  }

  km <- stats::kmeans(x, centers = k, nstart = spec@nstart)
  cluster_lum <- tapply(lum, km$cluster, mean)
  ord <- order(cluster_lum, decreasing = TRUE)

  fq_kmeans_analyzer(
    spec = spec,
    centers = km$centers[ord, , drop = FALSE],
    luminance = as.numeric(cluster_lum[ord])
  )
}

# Per-section metrics: area fraction in each grade, plus an area-weighted
# severity index on 0 (all mildest) to 10 (all most severe).
S7::method(fq_score, fq_kmeans_analyzer) <- function(fit, section, ...) {
  grade <- .assign(section, fit)
  k <- nrow(fit@centers)
  frac <- tabulate(grade, nbins = k) / length(grade)
  cols <- c(
    list(severity_index = sum(frac * (seq_len(k) - 1) / (k - 1)) * 10),
    setNames(as.list(frac), paste0("frac_sev_", seq_len(k)))
  )
  tibble::as_tibble(cols)
}

# Per-pixel severity field: each tissue pixel's grade in place, NA off tissue.
S7::method(fq_render, fq_kmeans_analyzer) <- function(fit, section, ...) {
  grade <- .assign(section, fit)
  field <- array(NA_real_, dim = dim(section@mask))
  field[section@mask] <- grade
  fq_field(values = field, k = nrow(fit@centers))
}
