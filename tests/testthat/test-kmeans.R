# test-kmeans.R

test_that("fq_kmeans builds a spec with sensible defaults", {
  spec <- fq_kmeans()
  expect_true(S7::S7_inherits(spec, fq_spec))
  expect_true(S7::S7_inherits(spec, fq_kmeans))
  expect_equal(spec@k, 3)
  expect_equal(spec@channels, c("a", "b"))
  expect_equal(spec@smooth_sigma, 2)
  expect_equal(spec@nstart, 3)
})

test_that("fq_kmeans accepts overrides", {
  spec <- fq_kmeans(k = 4, channels = c("L", "a", "b"), nstart = 10)
  expect_equal(spec@k, 4)
  expect_equal(spec@channels, c("L", "a", "b"))
  expect_equal(spec@nstart, 10)
})

test_that("fq_kmeans rejects invalid arguments", {
  expect_error(fq_kmeans(k = 1))
  expect_error(fq_kmeans(k = 2.5))
  expect_error(fq_kmeans(k = c(2, 3)))
  expect_error(fq_kmeans(channels = "x"))
  expect_error(fq_kmeans(channels = character(0)))
  expect_error(fq_kmeans(smooth_sigma = -1))
  expect_error(fq_kmeans(nstart = 0))
})
