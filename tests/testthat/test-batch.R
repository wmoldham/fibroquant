# test-batch.R

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
