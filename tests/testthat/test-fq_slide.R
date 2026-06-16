# test-fq_slide.R

# Self-contained synthetic pixels so each test carries its own fixture.
synth_rgb <- function(h = 6, w = 8) {
  array(
    runif(h * w * 3),
    dim = c(h, w, 3)
  )
}

test_that("fq_slide constructs and exposes its properties", {
  rgb <- synth_rgb()
  s <-
    fq_slide(
      rgb = rgb,
      um_per_px = 4,
      source = list(path = "synthetic")
    )
  expect_equal(s@rgb, rgb)
  expect_equal(s@um_per_px, 4)
  expect_equal(s@source$path, "synthetic")
})

test_that("fq_slide allows NA scale for uncalibrated images", {
  expect_no_error(
    fq_slide(
      rgb = synth_rgb(),
      um_per_px = NA_real_,
      source = list()
    )
  )
})

test_that("fq_slide rejects a non-RGB pixel array", {
  flat <- matrix(runif(6 * 8), 6, 8)
  expect_error(
    fq_slide(rgb = flat, um_per_px = 4, source = list()),
    "height x width x 3"
  )

  four_ch <- array(runif(6 * 8 * 4), dim = c(6, 8, 4))
  expect_error(
    fq_slide(rgb = four_ch, um_per_px = 4, source = list()),
    "height x width x 3"
  )
})

test_that("fq_slide rejects a non-positive scale", {
  expect_error(
    fq_slide(rgb = synth_rgb(), um_per_px = 0, source = list()),
    "positive"
  )
})

test_that("fq_section inherits fq_slide and adds section properties", {
  rgb <- synth_rgb(6, 8)
  mask <- matrix(
    FALSE,
    nrow = 6,
    ncol = 8
  )
  mask[2:5, 2:7] <- TRUE

  sec <- fq_section(
    rgb = rgb,
    um_per_px = 4,
    source = list(path = "synthetic"),
    mask = mask,
    bbox = list(rows = c(10L, 35L), cols = c(5L, 28L)),
    section = "A"
  )

  expect_true(S7::S7_inherits(sec, fq_slide)) # a section is a slide
  expect_equal(sec@mask, mask)
  expect_equal(sec@section, "A")
  expect_equal(sec@um_per_px, 4) # inherited property
})

test_that("fq_section rejects a mask that does not match the pixels", {
  rgb <- synth_rgb(6, 8)
  wrong <- matrix(
    FALSE,
    nrow = 5,
    ncol = 8
  )
  expect_error(
    fq_section(
      rgb = rgb,
      um_per_px = 4,
      source = list(),
      mask = wrong,
      bbox = list(),
      section = "A"
    ),
    "same height and width"
  )
})

# .fq_summary -----------------------------------------------------------------

test_that(".fq_summary describes a calibrated slide", {
  slide <-
    fq_slide(
      rgb = synth_rgb(6, 8),
      um_per_px = 4,
      source = list(path = "/data/Image_3470.vsi")
    )
  expect_equal(
    .fq_summary(slide),
    "<fq_slide> 6 × 8 px · 4 µm/px · Image_3470.vsi"
  )
})

test_that(".fq_summary marks an uncalibrated slide and an unknown source", {
  slide <-
    fq_slide(
      rgb = synth_rgb(6, 8),
      um_per_px = NA_real_,
      source = list()
    )
  expect_equal(
    .fq_summary(slide),
    "<fq_slide> 6 × 8 px · uncalibrated · unknown source"
  )
})

test_that(".fq_summary adds the label and tissue fraction for a section", {
  mask <-
    matrix(
      FALSE,
      nrow = 6,
      ncol = 8
    )
  mask[, 1:4] <- TRUE # half the pixels are tissue

  section <-
    fq_section(
      rgb = synth_rgb(6, 8),
      um_per_px = 4,
      source = list(path = "/data/Image_3470.vsi"),
      mask = mask,
      bbox = list(rows = c(1L, 6L), cols = c(1L, 8L)),
      section = "A"
    )
  expect_equal(
    .fq_summary(section),
    "<fq_section A> 6 × 8 px · 4 µm/px · 50% tissue · Image_3470.vsi"
  )
})

test_that("print.fq_slide writes the summary and returns its input invisibly", {
  slide <-
    fq_slide(
      rgb = synth_rgb(6, 8),
      um_per_px = 4,
      source = list(path = "x.png")
    )
  expect_output(print(slide), "<fq_slide> 6 × 8 px")

  result <- withVisible(print(slide))
  expect_false(result$visible)
  expect_identical(result$value, slide)
})
