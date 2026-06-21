# utils.R

# Tissue mask: pixels darker than an Otsu threshold. Tissue absorbs light.
# Glass and airspace stay bright.
.tissue_mask <- function(rgb) {
  luminance <- (rgb[, , 1] + rgb[, , 2] + rgb[, , 3]) / 3
  threshold <-
    EBImage::otsu(
      luminance,
      range = c(0, 1)
    )
  luminance < threshold
}

.microns_to_px <- function(um, um_per_px) {
  if (is.na(um_per_px)) {
    return(10L) # uncalibrated fallback
  }
  max(1L, as.integer(round(um / um_per_px)))
}
