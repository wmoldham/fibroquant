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

