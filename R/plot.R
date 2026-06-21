# plot.R

# Drawing helpers --------------------------------------------------------------

# Build a displayable colour Image from a slide's pixel array.
.fq_image <- function(rgb) {
  EBImage::Image(
    rgb,
    colormode = "Color"
  )
}

# Paint a rectangle border into an RGB array in image space, so it stays visible
# at any display scale.
.draw_box <- function(rgb, rows, cols, color, thickness) {
  r1 <- rows[1]
  r2 <- rows[2]
  c1 <- cols[1]
  c2 <- cols[2]
  for (ch in 1:3) {
    rgb[r1:min(r1 + thickness - 1L, r2), c1:c2, ch] <- color[ch]
    rgb[max(r2 - thickness + 1L, r1):r2, c1:c2, ch] <- color[ch]
    rgb[r1:r2, c1:min(c1 + thickness - 1L, c2), ch] <- color[ch]
    rgb[r1:r2, max(c2 - thickness + 1L, c1):c2, ch] <- color[ch]
  }
  rgb
}

# Overlay each section's crop rectangle on the parent slide's pixels.
.draw_section_boxes <- function(rgb, sections, col) {
  thickness <- max(3L, as.integer(round(max(dim(rgb)[1:2]) / 300)))
  rgb_col <- grDevices::col2rgb(col)[, 1] / 255
  for (section in sections) {
    bbox <- section@bbox
    rgb <- .draw_box(rgb, bbox$rows, bbox$cols, rgb_col, thickness)
  }
  rgb
}

# plot methods -----------------------------------------------------------------

S7::method(plot, fq_section) <- function(x, ...) {
  EBImage::display(
    .fq_image(x@rgb),
    method = "raster"
  )
  invisible(x)
}

# plot() on a slide shows the scan. Supplying `sections` (from fq_split())
# overlays each section's crop rectangle on the parent instead -- a quick check
# that the split caught every section and ignored streaks or debris.
S7::method(plot, fq_slide) <- function(x, sections = NULL, col = "red", ...) {
  if (is.null(sections)) {
    EBImage::display(
      .fq_image(x@rgb),
      method = "raster"
    )
    return(invisible(x))
  }

  boxed <- .draw_section_boxes(x@rgb, sections, col)
  EBImage::display(
    .fq_image(boxed),
    method = "raster"
  )
  invisible(x)
}

# plot() on an fq_sections collection is a contact sheet: each section in its own
# panel, labelled by section.
S7::method(plot, fq_sections) <- function(x, ...) {
  old <- graphics::par(mfrow = c(1, length(x)))
  on.exit(graphics::par(old), add = TRUE)
  for (section in x) {
    plot(section)
    graphics::title(section@section)
  }
  invisible(x)
}

# Map a severity field to an H x W x 3 RGB array. Each grade takes a colour
# sampled along a blue -> orange -> red severity ramp; off-tissue pixels are
# white.
.pseudocolor <- function(field) {
  values <- field@values
  pal <- grDevices::colorRampPalette(c("#1F74E0", "#F2891C", "#BE1622"))
  grade_rgb <- t(grDevices::col2rgb(pal(field@k)) / 255)
  out <- array(1, dim = c(dim(values), 3))
  tissue <- !is.na(values)
  g <- as.integer(values[tissue])
  for (ch in 1:3) {
    layer <- out[, , ch]
    layer[tissue] <- grade_rgb[g, ch]
    out[, , ch] <- layer
  }
  out
}

S7::method(plot, fq_field) <- function(x, ...) {
  EBImage::display(.fq_image(.pseudocolor(x)), method = "raster")
  invisible(x)
}
