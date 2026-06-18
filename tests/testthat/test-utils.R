# test-utils.R

test_that(".tissue_mask marks dark tissue and leaves bright background out", {
  rgb <-
    array(
      0.9,
      dim = c(10, 10, 3)
    )
  rgb[3:7, 3:7, ] <- 0.1 # a dark tissue block

  mask <- .tissue_mask(rgb)
  expect_equal(dim(mask), c(10L, 10L))
  expect_true(is.logical(mask))
  expect_true(all(mask[3:7, 3:7])) # the dark block is tissue
  expect_false(any(mask[1:2, ]))   # the bright margin is not
})

test_that(".microns_to_px converts, floors at 1, and falls back when uncalibrated", {
  expect_equal(.microns_to_px(50, 4), 12L)
  expect_equal(.microns_to_px(1, 4), 1L)
  expect_equal(.microns_to_px(50, NA_real_), 10L)
})
