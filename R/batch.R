# batch.R

#' Build a manifest of slide files for batch processing
#'
#' Scans a directory (or takes an explicit list of paths) for slide files
#' [fq_read()] can open and returns one row per file. The manifest fixes the set
#' of slides for a batch and gives each a stable `slide_id` for joining scores
#' back to their file and naming rendered maps.
#'
#' @param path A directory to scan, or a character vector of slide file paths.
#' @param recursive Recurse into sub-directories when `path` is a directory.
#' @param extensions Slide file extensions to match, lower-case, no leading dot.
#' @return A tibble with one row per slide: `slide_id` (file name without
#'   extension) and `path` (normalised absolute path), ordered by `path`.
#' @export
fq_manifest <- function(
    path,
    recursive = FALSE,
    extensions = c("vsi", "tif", "tiff", "png", "jpg", "jpeg")
) {
  if (length(path) == 1L && dir.exists(path)) {
    pattern <- sprintf("\\.(%s)$", paste(extensions, collapse = "|"))
    files <- list.files(
      path,
      pattern = pattern,
      recursive = recursive,
      ignore.case = TRUE,
      full.names = TRUE
    )
  } else {
    files <- path
  }

  if (length(files) == 0L) {
    stop("No slide files found.", call. = FALSE)
  }

  files <- sort(normalizePath(files, mustWork = FALSE))
  missing <- !file.exists(files)
  if (any(missing)) {
    stop(
      "No such file(s): ",
      paste(files[missing], collapse = ", "),
      call. = FALSE
    )
  }

  slide_id <- tools::file_path_sans_ext(basename(files))
  if (anyDuplicated(slide_id)) {
    warning(
      "Duplicate slide_id(s) disambiguated with a suffix; ",
      "rename files for stable ids.",
      call. = FALSE
    )
    slide_id <- make.unique(slide_id)
  }

  tibble::tibble(
    slide_id = slide_id,
    path = files
  )
}
