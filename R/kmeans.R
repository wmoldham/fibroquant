# kmeans.R

#' k-means colour-segmentation spec
#'
#' [Spec][fq_spec] for the k-means analyzer: cluster masked tissue pixels in
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
#'   to damp single-pixel stain speckle. `0` disables smoothing.
#' @param nstart Number of k-means restarts; the lowest within-cluster
#'   sum-of-squares fit is kept.
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
      nstart = S7::new_property(S7::class_numeric, default = 3)
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
      NULL
    }
  )

# Gaussian-blur each RGB channel to damp single-pixel stain speckle. sigma is in
# pixels; 0 leaves the image unchanged.
.smooth <- function(rgb, sigma) {
  if (sigma <= 0) {
    return(rgb)
  }
  EBImage::gblur(
    rgb,
    sigma = sigma
  )
}
