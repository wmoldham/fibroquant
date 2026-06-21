# test-batch.R

# Synthetic slide with two coloured tissue blocks, written to a raster file.
# The colour variation gives the clustering distinct points to work with.
write_coloured_two_section_png <- function(path) {
  rgb <- array(1, dim = c(120, 40, 3))
  rgb[10:40, 8:32, ] <- array(runif(31 * 25 * 3, 0.1, 0.6), dim = c(31, 25, 3))
  rgb[80:110, 8:32, ] <- array(runif(31 * 25 * 3, 0.1, 0.6), dim = c(31, 25, 3))
  EBImage::writeImage(
    EBImage::Image(rgb, colormode = "Color"),
    path
  )
}

# fq_manifest ------------------------------------------------------------------

test_that("fq_manifest lists matching files in a directory", {
  dir <- withr::local_tempdir()
  file.create(file.path(dir, c("Image_1.vsi", "Image_2.tif", "notes.txt")))

  man <- fq_manifest(dir)
  expect_equal(nrow(man), 2L)
  expect_equal(man$slide_id, c("Image_1", "Image_2"))
  expect_true(all(file.exists(man$path)))
})

test_that("fq_manifest accepts an explicit vector of paths", {
  dir <- withr::local_tempdir()
  paths <- file.path(dir, c("a.png", "b.png"))
  file.create(paths)

  expect_equal(fq_manifest(paths)$slide_id, c("a", "b"))
})

test_that("fq_manifest errors on an empty result and a missing file", {
  expect_error(fq_manifest(character(0)), "No slide files")
  expect_error(fq_manifest("/no/such/file.vsi"), "No such file")
})

test_that("fq_manifest disambiguates duplicate slide ids", {
  dir <- withr::local_tempdir()
  sub <- file.path(dir, "sub")
  dir.create(sub)
  file.create(file.path(dir, "Image_1.vsi"))
  file.create(file.path(sub, "Image_1.tif"))

  expect_warning(
    man <- fq_manifest(dir, recursive = TRUE),
    "Duplicate slide_id"
  )
  expect_equal(anyDuplicated(man$slide_id), 0L)
})

# fq_run -----------------------------------------------------------------------

test_that("fq_run scores every section and joins covariates", {
  dir <- withr::local_tempdir()
  ids <- c("s1", "s2", "s3")
  paths <- file.path(dir, paste0(ids, ".png"))
  for (p in paths) {
    write_coloured_two_section_png(p)
  }

  manifest <- fq_manifest(dir)
  manifest$treatment <- c("saline", "bleo", "bleo")

  out <- fq_run(manifest, fq_kmeans(k = 3), n = 2, n_ref = 2)

  expect_named(out, c("scores", "fit"))
  expect_true(S7::S7_inherits(out$fit, fq_kmeans_analyzer))
  expect_equal(nrow(out$scores), 6L) # 3 slides x 2 sections
  expect_true(
    all(
      c("slide_id", "section", "treatment", "severity_index") %in%
        names(out$scores)
    )
  )
  expect_setequal(unique(out$scores$slide_id), ids)
  expect_true(
    all(out$scores$severity_index >= 0 & out$scores$severity_index <= 10)
  )
})

test_that("fq_run errors on a manifest missing required columns", {
  bad <- tibble::tibble(slide_id = "x")
  expect_error(fq_run(bad, fq_kmeans()), "missing column")
})

test_that("fq_run reuses a fitted analyzer without refitting", {
  dir <- withr::local_tempdir()
  paths <- file.path(dir, c("s1.png", "s2.png"))
  for (p in paths) {
    write_coloured_two_section_png(p)
  }
  manifest <- fq_manifest(dir)

  fit <- fq_run(manifest, fq_kmeans(k = 3), n_ref = 2)$fit
  out <- fq_run(manifest, fit)

  expect_identical(out$fit, fit)
  expect_equal(nrow(out$scores), 4L)
})

test_that("fq_run skips a slide that fails to score, with a warning", {
  dir <- withr::local_tempdir()
  write_coloured_two_section_png(file.path(dir, "good.png"))
  writeLines("not an image", file.path(dir, "bad.png")) # unreadable as an image

  manifest <- fq_manifest(dir)
  good_only <- manifest[manifest$slide_id == "good", ]
  fit <- fq_run(good_only, fq_kmeans(k = 3), n_ref = 1)$fit

  expect_warning(
    out <- fq_run(manifest, fit),
    "Skipped 1 slide"
  )
  expect_setequal(unique(out$scores$slide_id), "good")
  expect_false(".error" %in% names(out$scores))
})

test_that("fq_run leaves the caller's RNG untouched", {
  dir <- withr::local_tempdir()
  for (id in c("s1", "s2")) {
    write_coloured_two_section_png(file.path(dir, paste0(id, ".png")))
  }
  manifest <- fq_manifest(dir)

  set.seed(42)
  before <- .Random.seed
  fq_run(manifest, fq_kmeans(k = 3), n_ref = 2)
  expect_identical(.Random.seed, before)
})
