#' Tabulate the image series and resolution levels in a slide file
#'
#' Whole-slide formats such as Olympus `.vsi` pack several images (a label, a
#' tissue overview, the high-resolution scan, a macro), and the scan itself is
#' stored as a pyramid of resolution levels. This returns one row per
#' (series, level) with its pixel dimensions and effective microns-per-pixel, so
#' [fq_read_slide()] can choose the scan series and a working resolution.
#'
#' @param path Path to a slide file readable by `RBioFormats` (e.g. `.vsi`).
#' @return A tibble with columns `series`, `res`, `size_x`, `size_y`, `um_px`,
#'   one row per resolution level. The finest level (`res == 1`) carries the
#'   native physical pixel size in `um_px`.
#' @export
fq_slide_info <- function(path) {
    .slide_info_table(RBioFormats::read.metadata(path))
  }

# Build the (series, level, dimensions, um/px) table from an RBioFormats
# metadata list.
.slide_info_table <- function(meta) {
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

    # Normalise to 1-based in case these fields arrive 0-based: series
    # globally, level within each series.
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

    # A coarser level's um/px is the native size scaled by its downsampling,
    # full_x / size_x, taken within each series.
    full_um <- .physical_um_per_px(meta)
    tab$um_px <- NA_real_
    for (s in unique(tab$series)) {
      in_series <- tab$series == s
      full_x <- tab$size_x[in_series][which.min(tab$res[in_series])]
      tab$um_px[in_series] <- full_um * (full_x / tab$size_x[in_series])
    }

    tab[order(tab$series, tab$res), ]
  }

# Parse the full-resolution physical pixel size (x, microns) from the global
# metadata, skipping series that lack it.
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
