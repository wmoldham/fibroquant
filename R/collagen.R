# collagen.R

# Classes ----------------------------------------------------------------------

#' Trichrome collagen (CPA) analyzer spec
#'
#' [spec][fq_spec] for the collagen analyzer: quantify Collagen Proportionate
#' Area (CPA), the fraction of tissue stained for collagen, on Masson's
#' trichrome sections. Collagen is isolated by colour deconvolution (Ruifrok &
#' Johnston) rather than by thresholding RGB directly, because the trichrome
#' dyes co-localize. The collagen (aniline-blue) density channel is then
#' thresholded and its area taken as a fraction of tissue. Pass the result to
#' [fq_fit()].
#'
#' This analyzer is **trichrome-specific**: it assumes a blue collagen stain
#' against a red counterstain. It is not meaningful on H&E, which does not
#' separate collagen by colour.
#'
#' @param collagen_quantile Quantile of the pooled collagen-density distribution
#'   used as the positive threshold, learned once on the fit batch. Higher is
#'   stricter. The default `0.98` was the most discriminating of several rules
#'   tested on a bleomycin series; collagen in lung is a density continuum
#'   rather than two separable populations, so a strict cutoff isolates
#'   pathological deposition from the faint basement-membrane signal present
#'   everywhere.
#' @param beta Optical-density floor, in OD units, below which a pixel is
#'   treated as too weakly stained to inform stain-vector estimation.
#' @param alpha Robust percentile (per side) for the stain-vector angle extremes
#'   in the colour-deconvolution estimate. Small values resist outliers.
#' @param max_px Maximum tissue pixels sampled per section when fitting, so one
#'   large section cannot dominate the stain estimate. `Inf` uses every pixel.
#' @return An `fq_collagen` spec.
#' @export
fq_collagen <-
  S7::new_class(
    "fq_collagen",
    parent = fq_spec,
    properties = list(
      collagen_quantile = S7::new_property(S7::class_numeric, default = 0.98),
      beta = S7::new_property(S7::class_numeric, default = 0.15),
      alpha = S7::new_property(S7::class_numeric, default = 1),
      max_px = S7::new_property(S7::class_numeric, default = 1e5)
    ),
    validator = function(self) {
      if (length(self@collagen_quantile) != 1L ||
          is.na(self@collagen_quantile) ||
          self@collagen_quantile <= 0 || self@collagen_quantile >= 1) {
        return("`collagen_quantile` must be a single number in (0, 1)")
      }
      if (length(self@beta) != 1L || is.na(self@beta) || self@beta < 0) {
        return("`beta` must be a single number >= 0")
      }
      if (length(self@alpha) != 1L || is.na(self@alpha) ||
          self@alpha <= 0 || self@alpha >= 50) {
        return("`alpha` must be a single number in (0, 50)")
      }
      if (length(self@max_px) != 1L || is.na(self@max_px) || self@max_px < 1) {
        return("`max_px` must be a single number >= 1 (Inf disables the cap)")
      }
      NULL
    }
  )

#' Fitted trichrome collagen analyzer
#'
#' The basis learned by [fq_fit()] from an [fq_collagen()] spec: the 3x3 stain
#' matrix estimated from the batch (rows are the unit optical-density vectors of
#' the collagen stain, the counterstain, and a residual), and the collagen
#' density threshold. Consumed by [fq_score()] and [fq_render()].
#'
#' @param spec The [fq_collagen()] spec this basis was fit from.
#' @param stain_matrix A 3x3 matrix whose rows are the unit optical-density
#'   vectors of the collagen stain, the counterstain, and their residual.
#' @param threshold Collagen-density threshold above which a pixel counts as
#'   collagen, learned as a quantile of the fit batch.
#' @return An `fq_collagen_analyzer` object.
#' @export
fq_collagen_analyzer <-
  S7::new_class(
    "fq_collagen_analyzer",
    parent = fq_analyzer,
    properties = list(
      spec = fq_collagen,
      stain_matrix = S7::class_numeric,
      threshold = S7::class_numeric
    )
  )

# Internal helpers -------------------------------------------------------------

# Reference collagen direction in optical-density space, the OD vector of a
# representative aniline-blue pixel. Used only to decide which of the two
# estimated stain vectors is collagen (the bluer one).
.collagen_ref <- local({
  v <- -log10(pmax(c(0.35, 0.45, 0.70), 1e-6))
  v / sqrt(sum(v^2))
})

# RGB in [0, 1] (array or matrix) to optical density, clamped to avoid log(0).
.rgb_to_od <- function(x) -log10(pmax(x, 1e-6))

# Masked tissue pixels of a section as an N x 3 optical-density matrix, sampled
# to at most max_px rows so each section contributes comparably to a fit.
.od_pixels <- function(section, max_px) {
  od <- matrix(.rgb_to_od(section@rgb), ncol = 3)
  od <- od[as.vector(section@mask), , drop = FALSE]
  if (nrow(od) > max_px) {
    od <- od[sample.int(nrow(od), as.integer(max_px)), , drop = FALSE]
  }
  od
}

# Macenko stain-vector estimation: the two stains are the robust extremes of the
# optical-density angle distribution in the plane of greatest variance. Returns
# a list of two unit OD vectors.
.macenko <- function(od, beta, alpha) {
  odhat <- od[!apply(od < beta, 1, any), , drop = FALSE]
  if (nrow(odhat) < 50L) {
    stop(
      "Too few strongly stained tissue pixels to estimate stain vectors.",
      call. = FALSE
    )
  }
  e <- eigen(stats::cov(odhat), symmetric = TRUE)
  v <- e$vectors[, 1:2]
  proj <- odhat %*% v
  phi <- atan2(proj[, 2], proj[, 1])
  lo <- stats::quantile(phi, alpha / 100)
  hi <- stats::quantile(phi, 1 - alpha / 100)
  unit <- function(x) {
    x <- x * sign(sum(x))
    x / sqrt(sum(x^2))
  }
  list(
    unit(as.numeric(v %*% c(cos(lo), sin(lo)))),
    unit(as.numeric(v %*% c(cos(hi), sin(hi))))
  )
}

# Build the 3x3 stain matrix from pooled tissue optical density. Rows are the
# unit OD vectors of collagen (the estimated vector nearest the blue reference),
# the counterstain, and their normalised cross product.
.stain_matrix <- function(od, beta, alpha) {
  vs <- .macenko(od, beta, alpha)
  d1 <- sum(vs[[1]] * .collagen_ref)
  d2 <- sum(vs[[2]] * .collagen_ref)
  collagen <- if (d1 >= d2) vs[[1]] else vs[[2]]
  counter <- if (d1 >= d2) vs[[2]] else vs[[1]]
  cx <- c(
    collagen[2] * counter[3] - collagen[3] * counter[2],
    collagen[3] * counter[1] - collagen[1] * counter[3],
    collagen[1] * counter[2] - collagen[2] * counter[1]
  )
  residual <- cx / sqrt(sum(cx^2))
  rbind(collagen = collagen, counter = counter, residual = residual)
}

# Collagen optical density per pixel via colour deconvolution: the collagen
# column of the inverse stain matrix applied to the OD image. Returns the values
# over masked tissue, or the full H x W grid (NA off tissue) when grid = TRUE.
.collagen_density <- function(section, stain, grid = FALSE) {
  inv_collagen <- solve(stain)[, 1]
  od <- matrix(.rgb_to_od(section@rgb), ncol = 3)
  dens <- as.numeric(od %*% inv_collagen)
  if (!grid) {
    return(dens[as.vector(section@mask)])
  }
  out <- array(NA_real_, dim = dim(section@mask))
  out[section@mask] <- dens[as.vector(section@mask)]
  out
}

# Methods ----------------------------------------------------------------------

# Pool masked tissue pixels across sections, estimate the stain matrix once, and
# set the collagen threshold as a quantile of the pooled collagen density. Cap
# each section to spec@max_px pixels first so one large section cannot dominate.
S7::method(fq_fit, fq_collagen) <- function(spec, sections, ...) {
  od <- do.call(rbind, lapply(sections, .od_pixels, max_px = spec@max_px))
  stain <- .stain_matrix(od, spec@beta, spec@alpha)
  dens <- as.numeric(od %*% solve(stain)[, 1])
  threshold <- as.numeric(stats::quantile(dens, spec@collagen_quantile))
  fq_collagen_analyzer(
    spec = spec,
    stain_matrix = stain,
    threshold = threshold
  )
}

# CPA for one section: the percentage of tissue pixels whose collagen density
# exceeds the fitted threshold.
S7::method(fq_score, fq_collagen_analyzer) <- function(fit, section, ...) {
  dens <- .collagen_density(section, fit@stain_matrix)
  cpa <- 100 * sum(dens > fit@threshold) / length(dens)
  tibble::tibble(severity_index = cpa)
}

# Collagen field. With density = FALSE the field is binary: grade 2 where tissue
# is collagen-positive, grade 1 elsewhere in tissue, NA off tissue. This is the
# pixels CPA counts. With density = TRUE the field grades the collagen channel
# itself into 32 levels, so plot() shows collagen density as a heatmap.
S7::method(fq_render, fq_collagen_analyzer) <-
  function(fit, section, density = FALSE, ...) {
    grid <- .collagen_density(section, fit@stain_matrix, grid = TRUE)
    if (density) {
      # Sequential blues, so the collagen channel reads as the isolated blue
      # stain. Clip the top 1% so a few strong pixels do not flatten the range.
      blues <- c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
      levels <- 32L
      hi <- stats::quantile(grid[section@mask], 0.99)
      scaled <- pmin(pmax(grid, 0) / hi, 1)
      field <- array(NA_real_, dim = dim(section@mask))
      field[section@mask] <- 1 + floor(scaled[section@mask] * (levels - 1))
      return(fq_field(values = field, k = levels, palette = blues))
    }
    # Binary mask on the default severity ramp: collagen red, tissue blue.
    field <- array(NA_real_, dim = dim(section@mask))
    field[section@mask] <- 1
    field[section@mask & grid > fit@threshold] <- 2
    fq_field(values = field, k = 2L)
  }
