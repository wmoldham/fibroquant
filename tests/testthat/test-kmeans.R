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

test_that(".smooth blurs each channel and preserves shape and range", {
  rgb <- array(0, dim = c(20, 20, 3))
  rgb[, 1:10, ] <- 1 # a sharp vertical edge

  sm <- .smooth(rgb, sigma = 2)
  expect_equal(dim(sm), dim(rgb))
  expect_true(sm[10, 10, 1] < 1 && sm[10, 11, 1] > 0) # edge softened
  expect_true(all(sm >= 0 & sm <= 1))
})

test_that(".smooth reduces pixel-to-pixel variance", {
  set.seed(1)
  rgb <- array(runif(30 * 30 * 3), dim = c(30, 30, 3))
  sm <- .smooth(rgb, sigma = 2)
  expect_lt(var(as.vector(sm)), var(as.vector(rgb)))
})

test_that(".smooth with sigma 0 returns the image unchanged", {
  rgb <- array(runif(10 * 10 * 3), dim = c(10, 10, 3))
  expect_identical(.smooth(rgb, sigma = 0), rgb)
})
