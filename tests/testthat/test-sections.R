# test-sections.R

# A synthetic wide scan: two dark blobs on white with a thin streak in the gap.
fake_two_section_slide <- function() {
  rgb <-
    array(
      1,
      dim = c(100, 40, 3)
    )
  rgb[10:35, 10:30, ] <- 0.2 # left section
  rgb[65:90, 10:30, ] <- 0.2 # right section
  rgb[49:51, 8:32, ] <- 0.2  # thin streak in the gap

  fq_slide(
    rgb = rgb,
    um_per_px = 4,
    source = list(
      path = "fake.vsi",
      format = "bioformats",
      series = 3L,
      resolution = 5L,
      native_um_px = 0.5,
      dims = c(100L, 40L)
    )
  )
}

test_that(".rank_components keeps the n largest and orders left to right", {
  labels <- matrix(0L, nrow = 20, ncol = 10)
  labels[2:6, 2:8] <- 1L   # left, area 35
  labels[14:18, 2:8] <- 2L # right, area 35
  labels[10, 5] <- 3L      # speck, area 1
  expect_equal(.rank_components(labels, 2, 0.1), c(1L, 2L))

  flipped <- matrix(0L, nrow = 20, ncol = 10)
  flipped[14:18, 2:8] <- 1L # label 1 on the right
  flipped[2:6, 2:8] <- 2L   # label 2 on the left
  expect_equal(.rank_components(flipped, 2, 0.1), c(2L, 1L)) # left first
})

test_that("fq_split returns an ordered fq_sections collection, dropping the streak", {
  sections <- fq_split(fake_two_section_slide(), close_um = 4)
  expect_true(S7::S7_inherits(sections, fq_sections))
  expect_length(sections, 2)
  expect_equal(
    vapply(sections, function(s) s@section, character(1)),
    c("A", "B")
  )
  expect_true(
    all(vapply(sections, function(s) S7::S7_inherits(s, fq_section), logical(1)))
  )
  expect_output(print(sections), "<fq_sections> 2 section")

  a <- sections[[1]]
  expect_equal(a@bbox$rows, c(10L, 35L))
  expect_equal(dim(a@rgb)[1:2], dim(a@mask))
  expect_equal(a@um_per_px, 4)
})

test_that("fq_split splits a real .vsi into sections (integration)", {
  skip_if_no_vsi()

  sections <- fq_read(vsi_path()) |> fq_split()
  expect_true(S7::S7_inherits(sections, fq_sections))
  expect_true(
    all(vapply(sections, function(s) S7::S7_inherits(s, fq_section), logical(1)))
  )
  expect_gte(length(sections), 1)
})

# A dark tissue blob with a bright airway lumen punched through it. The split
# fills the lumen into the footprint; the analysis mask must drop it again.
fake_airway_slide <- function() {
  rgb <-
    array(
      1,
      dim = c(60, 40, 3)
    )
  rgb[15:45, 10:30, ] <- 0.2   # tissue blob
  rgb[27:33, 17:23, ] <- 0.95  # airway lumen inside it
  fq_slide(
    rgb = rgb,
    um_per_px = 4,
    source = list(
      path = "fake.vsi",
      format = "bioformats",
      dims = c(60L, 40L)
    )
  )
}

test_that("fq_split excludes airway lumen from the analysis mask", {
  sec <- fq_split(
    fake_airway_slide(),
    n = 1,
    close_um = 4
  )[[1]]

  # The footprint keeps the filled lumen; the mask drops it.
  expect_gt(sum(sec@footprint), sum(sec@mask))

  # The lumen centre is inside the footprint, outside the mask.
  bbox <- sec@bbox
  lum_r <- 30 - bbox$rows[1] + 1
  lum_c <- 20 - bbox$cols[1] + 1
  expect_true(sec@footprint[lum_r, lum_c])
  expect_false(sec@mask[lum_r, lum_c])

  # The dropped pixels are bright lumen, not tissue.
  dropped <- sec@footprint & !sec@mask
  lum <- (sec@rgb[, , 1] + sec@rgb[, , 2] + sec@rgb[, , 3]) / 3
  expect_gt(mean(lum[dropped]), mean(lum[sec@mask]))
})

test_that("fq_split sections carry a footprint matching rgb", {
  sec <- fq_split(
    fake_airway_slide(),
    n = 1,
    close_um = 4
  )[[1]]
  expect_identical(dim(sec@footprint), dim(sec@rgb)[1:2])
  expect_type(sec@footprint, "logical")
})
