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

test_that(".lab converts to CIELAB with L*a*b* in channel order", {
  rgb <- array(0, dim = c(2, 2, 3))
  rgb[1, 1, ] <- c(1, 1, 1) # white
  rgb[1, 2, ] <- c(0, 0, 0) # black
  rgb[2, 1, ] <- c(1, 0, 0) # red
  rgb[2, 2, ] <- c(0, 0, 1) # blue

  lab <- .lab(rgb)
  expect_equal(dim(lab), c(2L, 2L, 3L))
  expect_gt(lab[1, 1, 1], 99) # white L* ~ 100
  expect_lt(lab[1, 2, 1], 1)  # black L* ~ 0
  expect_gt(lab[2, 1, 2], 0)  # red has positive a*
  expect_lt(lab[2, 2, 3], 0)  # blue has negative b*
})

test_that(".lab preserves non-square image shape", {
  rgb <- array(runif(3 * 5 * 3), dim = c(3, 5, 3))
  expect_equal(dim(.lab(rgb)), c(3L, 5L, 3L))
})

test_that(".features selects channels at masked pixels", {
  lab <- array(0, dim = c(2, 3, 3))
  lab[, , 1] <- 1 # L*
  lab[, , 2] <- 2 # a*
  lab[, , 3] <- 3 # b*
  mask <- matrix(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE), nrow = 2)

  feat <- .features(lab, mask, channels = c("a", "b"))
  expect_equal(dim(feat), c(sum(mask), 2L))
  expect_equal(colnames(feat), c("a", "b"))
  expect_true(all(feat[, "a"] == 2))
  expect_true(all(feat[, "b"] == 3))
})

test_that(".features pulls the correct masked pixels in order", {
  lab <- array(0, dim = c(2, 2, 3))
  lab[, , 2] <- matrix(c(11, 21, 12, 22), nrow = 2) # a* coded by position

  mask <- matrix(c(TRUE, FALSE, FALSE, TRUE), nrow = 2)
  feat <- .features(lab, mask, channels = "a")
  expect_equal(dim(feat), c(2L, 1L))
  expect_equal(as.vector(feat), c(11, 22)) # [1,1] then [2,2], column-major
})

test_that(".features respects channel order and names", {
  lab <- array(0, dim = c(2, 2, 3))
  lab[, , 1] <- 1 # L*
  lab[, , 2] <- 2 # a*
  lab[, , 3] <- 3 # b*
  mask <- matrix(TRUE, 2, 2)

  feat <- .features(lab, mask, channels = c("b", "L"))
  expect_equal(colnames(feat), c("b", "L"))
  expect_true(all(feat[, "b"] == 3))
  expect_true(all(feat[, "L"] == 1))
})
