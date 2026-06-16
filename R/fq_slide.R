# fq_slide.R

#' A whole-slide scan
#'
#' The core `fibroquant` object: a normalised RGB scan, its physical scale, and
#' the provenance needed to reproduce it. Every analyzer consumes an `fq_slide`
#' (or an [fq_section]).
#'
#' @param rgb Pixel data, a height x width x 3 numeric array in `[0, 1]`.
#' @param um_per_px Physical scale in microns per pixel; `NA` if uncalibrated.
#' @param source A list of provenance (path, format, series, resolution, native
#'   microns/pixel, dims) describing where this image came from.
#' @return An `fq_slide` object.
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
#' @inheritParams fq_slide
#' @param mask Logical tissue footprint, the same height x width as `rgb`; the
#'   denominator for density and the pixel pool for clustering.
#' @param bbox The crop box in parent coordinates, a list with `rows` and `cols`.
#' @param section The section label, e.g. `"A"` or `"B"`.
#' @return An `fq_section` object.
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

# One-line summary of an fq_slide or fq_section.
.fq_summary <- function(x) {
  dims <- dim(x@rgb)
  size <- paste0(dims[1], " \u00d7 ", dims[2], " px")

  scale <-
    if (is.na(x@um_per_px)) {
      "uncalibrated"
    } else {
      paste0(signif(x@um_per_px, 3), " \u00b5m/px")
    }

  path <- x@source$path
  src <-
    if (is.null(path) || !nzchar(path)) {
      "unknown source"
    } else {
      basename(path)
    }

  if (S7::S7_inherits(x, fq_section)) {
    tag <- paste0("<fq_section ", x@section, ">")
    tissue <- paste0(round(100 * mean(x@mask)), "% tissue")
    details <- c(size, scale, tissue, src)
  } else {
    tag <- "<fq_slide>"
    details <- c(size, scale, src)
  }

  joined <-
    paste(
      details,
      collapse = " \u00b7 "
    )
  paste(tag, joined)
}

S7::method(print, fq_slide) <- function(x, ...) {
  cat(
    .fq_summary(x),
    "\n",
    sep = ""
  )
  invisible(x)
}
