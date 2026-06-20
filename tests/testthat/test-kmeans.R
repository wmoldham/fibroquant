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

test_that(".smooth keeps airspace from bleeding into masked tissue", {
  rgb <-
    array(
      1,
      dim = c(20, 20, 3)
    )
  rgb[5:15, 5:15, ] <- 0.2   # tissue block on bright airspace
  mask <- matrix(FALSE, 20, 20)
  mask[5:15, 5:15] <- TRUE

  naive <-
    .smooth(
      rgb,
      sigma = 2
    )
  masked <-
    .smooth(
      rgb,
      sigma = 2,
      mask = mask
    )

  # Naive smoothing lightens the tissue rim toward the airspace.
  expect_gt(max(naive[, , 1][mask]), 0.2)

  # Mask-aware smoothing leaves uniform tissue uniform: airspace adds no weight.
  expect_equal(
    masked[, , 1][mask],
    rep(0.2, sum(mask)),
    tolerance = 1e-6
  )
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

# A minimal stand-in for fq_section: the kmeans methods only read @rgb and @mask.
section_stub <- S7::new_class(
  "section_stub",
  properties = list(
    rgb = S7::class_numeric,
    mask = S7::class_logical
  )
)

# Light cyan over dark red: two clusters separable in a*b*, cyan the brighter.
two_colour <- function(mask) {
  rgb <- array(0, dim = c(20, 20, 3))
  rgb[1:10, , 1] <- 0.6
  rgb[1:10, , 2] <- 0.9
  rgb[1:10, , 3] <- 0.9
  rgb[11:20, , 1] <- 0.4
  rgb[11:20, , 2] <- 0.1
  rgb[11:20, , 3] <- 0.1
  section_stub(rgb = rgb, mask = mask)
}

test_that("fq_fit on fq_kmeans returns a severity-ordered analyzer", {
  set.seed(1)
  fit <- fq_fit(fq_kmeans(k = 2, smooth_sigma = 0), list(two_colour(matrix(TRUE, 20, 20))))
  expect_true(S7::S7_inherits(fit, fq_kmeans_analyzer))
  expect_true(S7::S7_inherits(fit, fq_analyzer))
  expect_equal(dim(fit@centers), c(2L, 2L))
  expect_equal(length(fit@luminance), 2L)
  expect_true(all(diff(fit@luminance) <= 0)) # grade 1 brightest, k darkest
})

test_that("fq_fit on fq_kmeans pools across sections", {
  set.seed(1)
  make_section <- function() {
    rgb <- array(runif(15 * 15 * 3), dim = c(15, 15, 3))
    section_stub(rgb = rgb, mask = matrix(TRUE, 15, 15))
  }

  fit <- fq_fit(fq_kmeans(k = 3), list(make_section(), make_section()))
  expect_equal(dim(fit@centers), c(3L, 2L))
  expect_true(all(diff(fit@luminance) <= 0))
})

test_that("fq_fit on fq_kmeans errors with fewer pixels than clusters", {
  rgb <- array(0.5, dim = c(2, 2, 3))
  mask <- matrix(c(TRUE, FALSE, FALSE, FALSE), 2, 2) # one tissue pixel
  sec <- section_stub(rgb = rgb, mask = mask)
  expect_error(fq_fit(fq_kmeans(k = 3, smooth_sigma = 0), list(sec)))
})

test_that("fq_score returns area fractions and a severity index", {
  set.seed(1)
  sec <- two_colour(matrix(TRUE, 20, 20))
  fit <- fq_fit(fq_kmeans(k = 2, smooth_sigma = 0), list(sec))
  sc <- fq_score(fit, sec)

  expect_s3_class(sc, "tbl_df")
  expect_equal(nrow(sc), 1L)
  expect_named(sc, c("severity_index", "frac_sev_1", "frac_sev_2"))
  expect_equal(sc$frac_sev_1, 0.5) # half cyan -> grade 1
  expect_equal(sc$frac_sev_2, 0.5) # half red -> grade 2
  expect_equal(sc$severity_index, 5)
})

test_that("fq_score names fraction columns frac_sev_1..k", {
  set.seed(1)
  rgb <- array(runif(20 * 20 * 3), dim = c(20, 20, 3))
  sec <- section_stub(rgb = rgb, mask = matrix(TRUE, 20, 20))
  fit <- fq_fit(fq_kmeans(k = 3), list(sec))
  sc <- fq_score(fit, sec)

  expect_named(sc, c("severity_index", "frac_sev_1", "frac_sev_2", "frac_sev_3"))
  expect_equal(sc$frac_sev_1 + sc$frac_sev_2 + sc$frac_sev_3, 1)
})

test_that("fq_render returns a grade field with NA off tissue", {
  set.seed(1)
  mask <- matrix(TRUE, 20, 20)
  mask[, 1] <- FALSE # first column is background
  sec <- two_colour(mask)
  fit <- fq_fit(fq_kmeans(k = 2, smooth_sigma = 0), list(sec))
  fld <- fq_render(fit, sec)

  expect_true(S7::S7_inherits(fld, fq_field))
  expect_equal(dim(fld@values), dim(sec@mask))
  expect_equal(fld@k, 2L)
  expect_true(all(is.na(fld@values[, 1])))     # masked column
  expect_true(all(fld@values[1:10, 2:20] == 1)) # cyan -> grade 1
  expect_true(all(fld@values[11:20, 2:20] == 2)) # red -> grade 2
})

test_that("fq_kmeans carries and validates max_px", {
  expect_equal(fq_kmeans()@max_px, 1e5)
  expect_equal(fq_kmeans(max_px = 5e4)@max_px, 5e4)
  expect_error(fq_kmeans(max_px = 0), "max_px")
  expect_error(fq_kmeans(max_px = c(1, 2)), "max_px")
})
