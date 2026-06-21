# test-plot.R
# Fixtures: synth_section(), with_null_device() in helper-fixtures.R

# .draw_box --------------------------------------------------------------------

test_that(".draw_box paints a coloured border and leaves the rest untouched", {
  rgb <-
    array(
      0,
      dim = c(20, 20, 3)
    )
  out <- .draw_box(rgb, c(5, 15), c(5, 15), c(1, 0, 0), thickness = 2L)

  expect_true(all(out[5:6, 5:15, 1] == 1)) # top border is red
  expect_true(all(out[5:6, 5:15, 2] == 0))
  expect_true(all(out[8:12, 8:12, ] == 0)) # interior untouched
  expect_true(all(out[1:4, , ] == 0))      # exterior untouched
})

test_that(".draw_box respects thickness and clamps when it exceeds the box", {
  rgb <-
    array(
      0,
      dim = c(20, 20, 3)
    )

  thick <- .draw_box(rgb, c(5, 15), c(5, 15), c(1, 1, 1), thickness = 3L)
  expect_true(all(thick[5:7, 5:15, 1] == 1))
  expect_equal(thick[8, 8, 1], 0)

  clamped <- .draw_box(rgb, c(5, 8), c(5, 8), c(1, 0, 0), thickness = 10L)
  expect_true(all(clamped[1:4, , ] == 0))
  expect_true(all(clamped[9:20, , ] == 0))
})

test_that(".draw_box handles a box flush against the edge", {
  rgb <-
    array(
      0,
      dim = c(20, 20, 3)
    )
  expect_no_error(
    .draw_box(rgb, c(1, 20), c(1, 20), c(1, 0, 0), thickness = 2L)
  )
})

# .draw_section_boxes ----------------------------------------------------------

test_that(".draw_section_boxes paints a box at each section's bbox", {
  rgb <-
    array(
      0,
      dim = c(40, 30, 3)
    )
  out <- .draw_section_boxes(rgb, list(synth_section()), "red")

  expect_true(all(out[1, 1:10, 1] == 1)) # top edge of the box is red
  expect_true(all(out[1, 1:10, 2] == 0))
  expect_equal(out[20, 20, ], c(0, 0, 0)) # parent interior untouched
})

# .fq_image --------------------------------------------------------------------

test_that(".fq_image builds a colour Image", {
  rgb <-
    array(
      runif(12 * 10 * 3),
      dim = c(12, 10, 3)
    )
  img <- .fq_image(rgb)
  expect_s4_class(img, "Image")
  expect_equal(dim(img)[1:2], c(12L, 10L))
})

# .pseudocolor -----------------------------------------------------------------

test_that(".pseudocolor returns an RGB array with white off tissue", {
  fld <- fq_field(values = matrix(c(1, 2, NA, 3), 2, 2), k = 3L)
  px <- .pseudocolor(fld)

  expect_equal(dim(px), c(2L, 2L, 3L))
  expect_true(all(px >= 0 & px <= 1))
  expect_equal(px[1, 2, ], c(1, 1, 1)) # NA pixel -> white
})

test_that(".pseudocolor gives each grade a distinct colour", {
  fld <- fq_field(values = matrix(1:3, 1, 3), k = 3L)
  px <- .pseudocolor(fld)
  cols <- rbind(px[1, 1, ], px[1, 2, ], px[1, 3, ])
  expect_equal(nrow(unique(cols)), 3L)
})

# plot -------------------------------------------------------------------------

test_that("plot renders a section and a slide, returning input invisibly", {
  section <- synth_section()
  slide <-
    fq_slide(
      rgb = section@rgb,
      um_per_px = 4,
      source = list(path = "synthetic")
    )

  with_null_device({
    sec_result <- withVisible(plot(section))
    slide_result <- withVisible(plot(slide))
  })
  expect_false(sec_result$visible)
  expect_identical(sec_result$value, section)
  expect_false(slide_result$visible)
  expect_identical(slide_result$value, slide)
})

test_that("plot overlays section crop boxes when given sections", {
  slide <-
    fq_slide(
      rgb = array(1, dim = c(40, 30, 3)),
      um_per_px = 4,
      source = list(path = "synthetic")
    )
  sections <- fq_sections(list(synth_section(), synth_section()))

  with_null_device(expect_no_error(plot(slide, sections = sections)))
})

test_that("plot renders an fq_sections contact sheet, returning it invisibly", {
  sheet <- fq_sections(list(synth_section(), synth_section()))

  with_null_device({
    result <- withVisible(plot(sheet))
  })
  expect_false(result$visible)
  expect_identical(result$value, sheet)
})

test_that("plot(fq_field) renders without error", {
  fld <- fq_field(values = matrix(c(1, 2, NA, 3), 2, 2), k = 3L)
  with_null_device(expect_no_error(plot(fld)))
})
