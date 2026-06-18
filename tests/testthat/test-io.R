# test-io.R

# Synthetic RBioFormats-style metadata (0-based indices, to exercise the
# normalisation). One entry per (series, level); the physical pixel size sits on
# the first entry's global metadata, as it does in a real .vsi.
fake_meta <-
  function() {
    entry <-
      function(series, res, sx, sy, gm = NULL) {
        list(
          coreMetadata = list(
            series = series,
            resolutionLevel = res,
            sizeX = sx,
            sizeY = sy
          ),
          globalMetadata = gm
        )
      }
    list(
      entry(0, 0, 1000, 800, list("Physical pixel size" = "(0.5, 0.5)")),
      entry(0, 1, 500, 400),
      entry(1, 0, 4000, 2000),
      entry(1, 1, 2000, 1000)
    )
  }

test_that(".physical_um_per_px parses the size, scanning past series without it", {
  ok <- list(list(globalMetadata = list("Physical pixel size" = "(0.2738, 0.2739)")))
  expect_equal(.physical_um_per_px(ok), 0.2738)

  later <-
    list(
      list(globalMetadata = NULL),
      list(globalMetadata = list("Other" = "x")),
      list(globalMetadata = list("Physical Pixel Size" = "(1.5, 1.5)"))
    )
  expect_equal(.physical_um_per_px(later), 1.5)

  missing <- list(list(globalMetadata = list("Nope" = "1")))
  expect_true(is.na(.physical_um_per_px(missing)))
})

test_that(".slide_info_table normalises indices and computes um/px per level", {
  info <- .slide_info_table(fake_meta())

  expect_equal(sort(unique(info$series)), c(1L, 2L)) # 0-based -> 1-based
  expect_equal(min(info$res), 1L)
  expect_equal(min(info$um_px), 0.5) # finest level = native scale

  s2 <- info[info$series == 2, ]
  s2 <- s2[order(s2$res), ]
  expect_equal(s2$um_px, c(0.5, 1.0)) # pixels double in size each level down
})

test_that("fq_slide_info reads a real .vsi (integration)", {
  skip_if_no_vsi()

  info <- fq_slide_info(vsi_path())
  expect_true(all(c("series", "res", "size_x", "size_y", "um_px") %in% names(info)))
  expect_gt(nrow(info), 0)
  expect_true(all(info$um_px > 0))
})


test_that(".pick_scan_series picks the largest series despite integer overflow", {
  info <- tibble::tibble(
    series = c(1L, 1L, 2L, 2L, 3L, 3L),
    res = c(1L, 2L, 1L, 2L, 1L, 2L),
    size_x = c(8021L, 4011L, 18032L, 9016L, 140428L, 70214L),
    size_y = c(9366L, 4683L, 9148L, 4574L, 21575L, 10788L)
  )
  # series 3 full-res area ~3.03e9 overflows 32-bit int; it must still win
  expect_equal(.pick_scan_series(info), 3L)
})

test_that(".pick_resolution chooses the level closest to the target um/px", {
  info <- tibble::tibble(
    series = rep(3L, 4),
    res = 1:4,
    size_x = c(140428L, 70214L, 35107L, 17554L),
    size_y = c(21575L, 10788L, 5394L, 2697L),
    um_px = c(0.274, 0.548, 1.10, 2.19)
  )
  expect_equal(.pick_resolution(info, series = 3, target_um_px = 1), 3L)
  expect_equal(.pick_resolution(info, series = 3, target_um_px = 2), 4L)
})

test_that(".as_rgb_array returns a height x width x 3 array in [0, 1]", {
  arr <- array(
    runif(10 * 12 * 3),
    dim = c(10, 12, 3)
  )
  img <- EBImage::Image(
    arr,
    colormode = "Color"
  )
  out <- .as_rgb_array(img)
  expect_equal(dim(out), c(10L, 12L, 3L))
  expect_true(all(out >= 0 & out <= 1))

  gray <- EBImage::Image(matrix(runif(10 * 12), 10, 12))
  expect_equal(dim(.as_rgb_array(gray)), c(10L, 12L, 3L)) # grayscale -> 3 channels
})

test_that("fq_read reads a standard raster into an fq_slide", {
  arr <-
    array(
      runif(8 * 6 * 3),
      dim = c(8, 6, 3)
    )
  img <-
    EBImage::Image(
      arr,
      colormode = "Color"
    )
  path <- withr::local_tempfile(fileext = ".png")
  EBImage::writeImage(img, path)

  s <- fq_read(path)
  expect_true(S7::S7_inherits(s, fq_slide))
  expect_equal(dim(s@rgb), c(8L, 6L, 3L))
  expect_equal(s@source$format, "ebimage")
  expect_true(is.na(s@um_per_px))
  expect_true(all(s@rgb >= 0 & s@rgb <= 1))
})

test_that("fq_read reads a real .vsi (integration)", {
  skip_if_no_vsi()

  s <- fq_read(vsi_path())
  expect_true(S7::S7_inherits(s, fq_slide))
  expect_equal(length(dim(s@rgb)), 3L)
  expect_gt(s@um_per_px, 0)
  expect_equal(s@source$format, "bioformats")
})
