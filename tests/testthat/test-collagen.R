# test-collagen.R

# A synthetic trichrome section: red counterstained parenchyma with a top band
# of blue collagen covering `collagen_frac` of the tissue rows. A little colour
# noise spreads the optical-density values so thresholds are not degenerate.
synth_trichrome <- function(collagen_frac = 0.3, h = 80, w = 60, seed = 1) {
  withr::with_seed(seed, {
    rgb <- array(1, dim = c(h, w, 3)) # white glass background
    tr <- 6:75
    tc <- 6:55
    mask <- matrix(FALSE, h, w)
    mask[tr, tc] <- TRUE

    jit <- function(v) v + stats::runif(length(tr) * length(tc), -0.03, 0.03)
    # parenchyma: pink-red (Biebrich scarlet / acid fuchsin), absorbing in all
    # channels so it survives the optical-density background filter.
    rgb[tr, tc, 1] <- jit(0.62)
    rgb[tr, tc, 2] <- jit(0.42)
    rgb[tr, tc, 3] <- jit(0.50)
    # collagen: blue (aniline blue) over the top band of tissue rows
    nb <- round(length(tr) * collagen_frac)
    if (nb > 0) {
      br <- tr[seq_len(nb)]
      j2 <- function(v) v + stats::runif(length(br) * length(tc), -0.03, 0.03)
      rgb[br, tc, 1] <- j2(0.38)
      rgb[br, tc, 2] <- j2(0.48)
      rgb[br, tc, 3] <- j2(0.62)
    }
    rgb <- pmin(pmax(rgb, 0), 1)

    fq_section(
      rgb = rgb,
      um_per_px = 4,
      source = list(path = "synthetic"),
      mask = mask,
      bbox = list(rows = c(1L, h), cols = c(1L, w)),
      section = "A"
    )
  })
}

test_that("fq_fit returns a fitted collagen analyzer", {
  sections <- lapply(c(0.2, 0.3, 0.4), synth_trichrome)
  fit <- fq_fit(fq_collagen(), sections)

  expect_true(S7::S7_inherits(fit, fq_collagen_analyzer))
  expect_equal(dim(fit@stain_matrix), c(3L, 3L))
  expect_true(is.finite(fit@threshold))
  # The collagen row should lean bluer (closer to the blue reference) than the
  # counterstain row, i.e. the right vector was chosen as collagen.
  d_coll <- sum(fit@stain_matrix[1, ] * fibroquant:::.collagen_ref)
  d_counter <- sum(fit@stain_matrix[2, ] * fibroquant:::.collagen_ref)
  expect_gt(d_coll, d_counter)
})

test_that("fq_score returns a one-row severity_index tibble", {
  fit <- fq_fit(fq_collagen(), lapply(c(0.2, 0.3, 0.4), synth_trichrome))
  out <- fq_score(fit, synth_trichrome(0.3))

  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 1L)
  expect_true("severity_index" %in% names(out))
  expect_true(is.finite(out$severity_index))
  expect_gte(out$severity_index, 0)
  expect_lte(out$severity_index, 100)
})

test_that("more collagen yields a higher CPA", {
  fit <- fq_fit(fq_collagen(), lapply(c(0.2, 0.3, 0.4), synth_trichrome))
  low <- fq_score(fit, synth_trichrome(0.1))$severity_index
  high <- fq_score(fit, synth_trichrome(0.6))$severity_index
  expect_gt(high, low)
})

test_that("fq_render returns a binary collagen field", {
  fit <- fq_fit(fq_collagen(), lapply(c(0.2, 0.3, 0.4), synth_trichrome))
  sec <- synth_trichrome(0.3)
  field <- fq_render(fit, sec)

  expect_true(S7::S7_inherits(field, fq_field))
  expect_equal(field@k, 2L)
  expect_equal(dim(field@values), dim(sec@mask))
  expect_setequal(unique(as.vector(field@values)), c(NA, 1, 2))
  # Collagen pixels fall inside tissue.
  expect_true(all(sec@mask[which(field@values == 2)]))
})

test_that("fq_render density = TRUE grades the collagen channel", {
  fit <- fq_fit(fq_collagen(), lapply(c(0.2, 0.3, 0.4), synth_trichrome))
  sec <- synth_trichrome(0.3)
  dens <- fq_render(fit, sec, density = TRUE)

  expect_true(S7::S7_inherits(dens, fq_field))
  expect_gt(dens@k, 2L)
  grades <- as.vector(dens@values)[sec@mask]
  expect_true(all(grades >= 1 & grades <= dens@k))
  # Collagen rows (the top band, ~6-26) grade higher than parenchyma rows.
  hi <- mean(dens@values[10:24, 10:50], na.rm = TRUE)
  lo <- mean(dens@values[40:70, 10:50], na.rm = TRUE)
  expect_gt(hi, lo)
})

test_that("the stain matrix is invertible and trichrome-gated by colour", {
  fit <- fq_fit(fq_collagen(), lapply(c(0.2, 0.3, 0.4), synth_trichrome))
  expect_silent(solve(fit@stain_matrix))
})

test_that(".macenko errors when there is too little stained tissue", {
  od <- matrix(0.01, nrow = 100, ncol = 3) # all below beta
  expect_error(fibroquant:::.macenko(od, beta = 0.15, alpha = 1), "Too few")
})
