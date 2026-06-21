# sections.R

# fq_split ---------------------------------------------------------------------

#' Split a multi-section slide into one fq_section per tissue section
#'
#' Lung slides routinely carry several sections side by side. This finds them by
#' connected components -- consolidating each porous section into a solid blob so
#' alveolar airspace does not fragment it -- and returns one cropped [fq_section]
#' per section, ordered along the first spatial axis (left to right for a wide
#' strip).
#'
#' @param slide An [fq_slide] from [fq_read()].
#' @param n Number of sections to keep (largest by area).
#' @param close_um Morphological closing radius in microns; bridges airspace
#'   within a section without bridging the gap between sections.
#' @param min_area_frac Drop components smaller than this fraction of the
#'   largest, so debris and scan artifacts are ignored.
#' @return An [fq_sections] collection (a list of [fq_section]s, ordered left to
#'   right).
#' @export
fq_split <- function(slide,
                     n = 2,
                     close_um = 50,
                     min_area_frac = 0.05) {
  foreground <- .tissue_mask(slide@rgb)

  # Close each porous section into a solid blob, then fill alveolar holes; the
  # wide inter-section gap is too large to bridge and stays open.
  close_px <- .microns_to_px(close_um, slide@um_per_px)
  brush <-
    EBImage::makeBrush(
      2L * close_px + 1L,
      shape = "disc"
    )
  consolidated <- (foreground * 1) |>
    EBImage::closing(brush) |>
    EBImage::fillHull()

  labels <- EBImage::bwlabel(consolidated)
  keep <- .rank_components(labels, n, min_area_frac)
  if (length(keep) < n) {
    warning(
      "Found ", length(keep), " section(s), expected ", n, ".",
      call. = FALSE
    )
  }

  sections <-
    lapply(
      seq_along(keep),
      function(i) .extract_section(slide, labels, keep[i], i)
    )
  fq_sections(sections)
}

# Internal helpers -------------------------------------------------------------

# Labels of the components to keep: the n largest above min_area_frac of the
# biggest, reordered along the first spatial axis (left to right).
.rank_components <- function(labels, n, min_area_frac) {
  idx <- which(labels > 0)
  if (length(idx) == 0) {
    return(integer(0))
  }
  comp <- labels[idx]
  area <- tabulate(comp)
  big_enough <- which(area >= min_area_frac * max(area))
  ranked <- big_enough[order(area[big_enough], decreasing = TRUE)]
  keep <- ranked[seq_len(min(n, length(ranked)))]

  positions <- arrayInd(idx, dim(labels))[, 1]
  centroid_x <- tapply(positions, comp, mean)
  keep[order(centroid_x[as.character(keep)])]
}

# Crop the slide to a section's bounding box and build an fq_section.
.extract_section <- function(slide, labels, label, index) {
  region <- labels == label
  rows <- range(which(rowSums(region) > 0))
  cols <- range(which(colSums(region) > 0))

  rgb <- slide@rgb[rows[1]:rows[2], cols[1]:cols[2], , drop = FALSE]
  footprint <- region[rows[1]:rows[2], cols[1]:cols[2], drop = FALSE]

  # Re-threshold the cropped section to drop the airway and whitespace lumen the
  # split filled over; intersect with the footprint to hold the section boundary.
  mask <- .tissue_mask(rgb) & footprint

  fq_section(
    rgb = rgb,
    um_per_px = slide@um_per_px,
    source = list(
      path = slide@source$path,
      format = slide@source$format,
      series = slide@source$series,
      resolution = slide@source$resolution,
      native_um_px = slide@source$native_um_px,
      parent_dims = slide@source$dims,
      dims = dim(rgb)[1:2]
    ),
    mask = mask,
    footprint = footprint,
    bbox = list(rows = rows, cols = cols),
    section = LETTERS[index]
  )
}
