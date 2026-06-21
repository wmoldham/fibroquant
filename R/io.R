# io.R

# fq_info ----------------------------------------------------------------------

#' Tabulate the image series and resolution levels in a slide file
#'
#' Whole-slide formats such as Olympus `.vsi` pack several images, including a
#' label, a tissue overview, the full-resolution scan, and a macro image. The
#' scan itself is stored as a pyramid of resolution levels. This returns one row
#' per series and level, with its pixel dimensions and effective microns per
#' pixel, so [fq_read()] can choose the scan series and a working resolution.
#'
#' @param path Path to a slide file readable by `RBioFormats` (e.g. `.vsi`).
#' @return A tibble with one row per resolution level and columns `series`,
#'   `res` (the level, where 1 is finest), `size_x` and `size_y` (pixel
#'   dimensions), and `um_px` (effective microns per pixel). At the finest level
#'   (`res == 1`), `um_px` is the native physical pixel size.
#' @export
fq_info <- function(path) {
  .info_table(RBioFormats::read.metadata(path))
}

# Build the (series, level, dimensions, um/px) table from an RBioFormats
# metadata list.
.info_table <- function(meta) {
  if (methods::is(meta, "ImageMetadata")) {
    meta <- list(meta)
  }
  core <- lapply(meta, function(m) m$coreMetadata)

  field <- function(name) {
    vapply(
      core,
      function(m) as.integer(m[[name]]),
      integer(1)
    )
  }

  # Normalise indices to start at 1 in case they arrive 0 based. Series is
  # global. Level resets within each series.
  series <- field("series")
  series <- series - min(series) + 1L

  res_raw <- field("resolutionLevel")
  res <- res_raw
  for (s in unique(series)) {
    in_series <- series == s
    res[in_series] <- res_raw[in_series] - min(res_raw[in_series]) + 1L
  }

  tab <-
    tibble::tibble(
      series = series,
      res = res,
      size_x = field("sizeX"),
      size_y = field("sizeY")
    )

  # A coarser level's um/px is the native scale times its downsampling,
  # full_x / size_x, within each series.
  full_um <- .physical_um_per_px(meta)
  tab$um_px <- NA_real_
  for (s in unique(tab$series)) {
    in_series <- tab$series == s
    full_x <- tab$size_x[in_series][which.min(tab$res[in_series])]
    tab$um_px[in_series] <- full_um * (full_x / tab$size_x[in_series])
  }

  tab[order(tab$series, tab$res), ]
}

# Parse the physical pixel size (x, in microns) from the global metadata.
# Skip series that lack it.
.physical_um_per_px <- function(meta) {
  for (m in meta) {
    gm <- m$globalMetadata
    if (is.null(gm)) {
      next
    }
    hit <- grep(
      "physical pixel size",
      names(gm),
      ignore.case = TRUE
    )
    if (length(hit) > 0) {
      clean <- gsub("[()]", "", gm[[hit[1]]])
      parts <- strsplit(clean, ",")[[1]]
      return(as.numeric(parts)[1])
    }
  }
  NA_real_
}

# fq_read ----------------------------------------------------------------------

# The series with the largest pixel area is the scan. as.numeric guards against
# integer overflow on the pixel counts.
.pick_scan_series <- function(info) {
  (as.numeric(info$size_x) * as.numeric(info$size_y)) |>
    tapply(info$series, max) |>
    which.max() |>
    names() |>
    as.integer()
}

# Resolution level whose effective um/px is closest to the target.
.pick_resolution <- function(info, series, target_um_px) {
  sub <- info[info$series == series, , drop = FALSE]
  sub$res[which.min(abs(sub$um_px - target_um_px))]
}

# Coerce an RBioFormats or EBImage image to an H x W x 3 array in [0, 1].
# Replicate a single channel to RGB. Drop any alpha channel.
.as_rgb_array <- function(img) {
  arr <- as.array(EBImage::imageData(img))
  if (length(dim(arr)) == 2L) {
    arr <-
      array(
        arr,
        dim = c(dim(arr), 3L)
      )
  }
  arr[, , seq_len(3L), drop = FALSE]
}

#' Read a slide at a working resolution
#'
#' Reads a whole-slide image into an [fq_slide]. For formats `RBioFormats`
#' handles, such as `.vsi`, it selects the scan series and the resolution level
#' nearest `target_um_px`. TIFF, PNG, and JPEG are read directly by EBImage as a
#' single, uncalibrated level.
#'
#' @param path Path to a slide file.
#' @param series Series index. `NULL` selects the scan series automatically, the
#'   one largest by area.
#' @param target_um_px Desired microns per pixel. The level nearest this value
#'   is chosen.
#' @param resolution An explicit resolution level. Overrides `target_um_px` when
#'   given.
#' @return An [fq_slide].
#' @export
fq_read <- function(
    path,
    series = NULL,
    target_um_px = 4,
    resolution = NULL
) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tif", "tiff", "png", "jpg", "jpeg")) {
    return(.read_ebimage(path))
  }

  info <- fq_info(path)
  if (is.null(series)) {
    series <- .pick_scan_series(info)
  }
  if (is.null(resolution)) {
    resolution <- .pick_resolution(info, series, target_um_px)
  }

  img <-
    RBioFormats::read.image(
      path,
      series = series,
      resolution = resolution,
      normalize = TRUE
    )
  rgb <- .as_rgb_array(img)
  in_series <- info$series == series
  eff_um_px <- info$um_px[in_series & info$res == resolution][1]
  native_um_px <- info$um_px[in_series & info$res == 1L][1]

  fq_slide(
    rgb = rgb,
    um_per_px = eff_um_px,
    source = list(
      path = path,
      format = "bioformats",
      series = series,
      resolution = resolution,
      native_um_px = native_um_px,
      dims = dim(rgb)[1:2]
    )
  )
}

.read_ebimage <- function(path, um_per_px = NA_real_) {
  rgb <- .as_rgb_array(EBImage::readImage(path))

  fq_slide(
    rgb = rgb,
    um_per_px = um_per_px,
    source = list(
      path = path,
      format = "ebimage",
      series = NA_integer_,
      resolution = NA_integer_,
      native_um_px = um_per_px,
      dims = dim(rgb)[1:2]
    )
  )
}
