# fq_slide.R

# Classes ----------------------------------------------------------------------

#' A whole-slide scan
#'
#' The core `fibroquant` object. It holds a normalised RGB scan, its physical
#' scale, and the provenance needed to reproduce it. Every analyzer operates on
#' an `fq_slide` or an [fq_section].
#'
#' @param rgb Pixel data, a height x width x 3 numeric array in `[0, 1]`.
#' @param um_per_px Physical scale in microns per pixel. `NA` if uncalibrated.
#' @param source A named list of provenance describing where the image came
#'   from, with entries such as path, format, series, resolution, native microns
#'   per pixel, and dimensions.
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
#' One lung section extracted by [fq_split()]. It inherits everything from
#' [fq_slide] and adds the analysis mask, the filled silhouette, and the crop
#' location, so a section carries both its pixels and where it sat in the parent
#' scan.
#'
#' @inheritParams fq_slide
#' @param mask Logical analysis mask, the same height and width as `rgb`. `TRUE`
#'   marks foreground tissue and excludes the airway and other bright lumen.
#'   This is the denominator for tissue density and the pool of pixels for
#'   clustering.
#' @param footprint Logical section silhouette, the same height and width as
#'   `rgb`. This is the filled connected component that sets the crop box and
#'   the drawn outline. It defaults to `mask` when a section has no distinct
#'   silhouette.
#' @param bbox The crop box in the parent scan's coordinates, a list with `rows`
#'   and `cols`, each a length 2 range.
#' @param section The section label, e.g. `"A"` or `"B"`.
#' @return An `fq_section` object.
#' @export
fq_section <-
  S7::new_class(
    "fq_section",
    parent = fq_slide,
    properties = list(
      mask = S7::class_logical,
      footprint = S7::class_logical,
      bbox = S7::class_list,
      section = S7::class_character
    ),
    constructor = function(rgb,
                           um_per_px = NA_real_,
                           source = list(),
                           mask,
                           bbox = list(),
                           section = NA_character_,
                           footprint = mask) {
      S7::new_object(
        fq_slide(
          rgb = rgb,
          um_per_px = um_per_px,
          source = source
        ),
        mask = mask,
        footprint = footprint,
        bbox = bbox,
        section = section
      )
    },
    validator = function(self) {
      if (length(dim(self@mask)) != 2L) {
        "@mask must be a 2-D logical matrix"
      } else if (!identical(dim(self@mask), dim(self@rgb)[1:2])) {
        "@mask must have the same height and width as @rgb"
      } else if (!identical(dim(self@footprint), dim(self@rgb)[1:2])) {
        "@footprint must have the same height and width as @rgb"
      } else if (length(self@section) != 1L) {
        "@section must be a single label"
      }
    }
  )

# Summary and printing ---------------------------------------------------------

# Summarise an fq_slide or fq_section in one line.
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

# fq_sections ------------------------------------------------------------------

#' A collection of tissue sections from one slide
#'
#' The return type of [fq_split()]. It extends a list, so it indexes, lengths,
#' and iterates like one, while carrying its own [plot()] method (a contact
#' sheet) and [print()] method.
#'
#' @param sections A list of [fq_section] objects.
#' @return An `fq_sections` object.
#' @export
fq_sections <-
  S7::new_class(
    "fq_sections",
    parent = S7::class_list,
    constructor = function(sections = list()) {
      S7::new_object(sections)
    }
  )

S7::method(print, fq_sections) <- function(x, ...) {
  cat("<fq_sections>", length(x), "section(s)\n")
  for (section in x) {
    cat(paste0("  ", .fq_summary(section), "\n"))
  }
  invisible(x)
}
