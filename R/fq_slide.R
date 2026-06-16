# fq_slide.R

#' A whole-slide scan
#'
#' The core `fibroquant` object: a normalised RGB scan, its physical scale, and
#' the provenance needed to reproduce it. Every analyzer consumes an `fq_slide`
#' (or an [fq_section]).
#'
#' @section Properties:
#' - `rgb`: pixel data, a height x width x 3 numeric array in `[0, 1]`.
#' - `um_per_px`: physical scale in microns per pixel; `NA` if uncalibrated.
#' - `source`: a list of provenance (path, format, series, resolution, native
#'   microns/pixel, dims) describing where this image came from.
#'
#' @export
fq_slide <-
  S7::new_class(
    "fq_slide",
    properties = list(
      rgb = S7::class_double,
      um_per_px = S7::class_numeric,
      source = S7::class_list
    ),
    validator = function(self) {
      if (length(dim(self@rgb)) != 3L || dim(self@rgb)[3L] != 3L) {
        "@rgb must be a height x width x 3 array"
      } else if (length(self@um_per_px) != 1L) {
        "@um_per_px must be a single value"
      } else if (!is.na(self@um_per_px) && self@um_per_px <= 0) {
        "@um_per_px must be positive, or NA if uncalibrated"
      }
    }
  )

#' A single tissue section cropped from a slide
#'
#' One lung section extracted by [fq_split_sections()]. Inherits everything from
#' [fq_slide] and adds the tissue footprint and the crop provenance, so a section
#' carries both its pixels and where it sat in the parent scan.
#'
#' @section Additional properties:
#' - `mask`: logical tissue footprint, the same height x width as `rgb`; the
#'   denominator for density and the pixel pool for clustering.
#' - `bbox`: the crop box in parent coordinates, a list with `rows` and `cols`.
#' - `section`: the section label, e.g. `"A"` or `"B"`.
#'
#' @export
fq_section <-
  S7::new_class(
    "fq_section",
    parent = fq_slide,
    properties = list(
      mask = S7::class_logical,
      bbox = S7::class_list,
      section = S7::class_character
    ),
    validator = function(self) {
      if (length(dim(self@mask)) != 2L) {
        "@mask must be a 2-D logical matrix"
      } else if (!identical(dim(self@mask), dim(self@rgb)[1:2])) {
        "@mask must have the same height and width as @rgb"
      } else if (length(self@section) != 1L) {
        "@section must be a single label"
      }
    }
  )
