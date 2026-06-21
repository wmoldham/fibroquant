# helper-fixtures.R
#
# Synthetic fixtures shared across test files. testthat sources helper files
# before any test, so each test file does not need its own copy.

# Random RGB pixel array, height x width x 3 in [0, 1].
synth_rgb <- function(h = 6, w = 8) {
  array(runif(h * w * 3), dim = c(h, w, 3))
}

# A small fq_section: a dark tissue block on a white field, with the mask
# covering the block.
synth_section <- function(h = 12, w = 10) {
  rgb <- array(1, dim = c(h, w, 3))
  rgb[3:9, 3:7, ] <- 0.2
  mask <- matrix(FALSE, nrow = h, ncol = w)
  mask[3:9, 3:7] <- TRUE

  fq_section(
    rgb = rgb,
    um_per_px = 4,
    source = list(path = "synthetic"),
    mask = mask,
    bbox = list(rows = c(1L, h), cols = c(1L, w)),
    section = "A"
  )
}

# A minimal stub for fq_section. The kmeans methods only read @rgb and @mask.
section_stub <- S7::new_class(
  "section_stub",
  properties = list(
    rgb = S7::class_numeric,
    mask = S7::class_logical
  )
)

# Light cyan over dark red: two clusters separable in a*b*, cyan the brighter.
two_colour <- function(mask) {
  rgb <- array(0, dim = c(20, 20, 3))
  rgb[1:10, , 1] <- 0.6
  rgb[1:10, , 2] <- 0.9
  rgb[1:10, , 3] <- 0.9
  rgb[11:20, , 1] <- 0.4
  rgb[11:20, , 2] <- 0.1
  rgb[11:20, , 3] <- 0.1
  section_stub(rgb = rgb, mask = mask)
}

# Send graphics to a throwaway device so plotting tests stay headless.
with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(code)
}
