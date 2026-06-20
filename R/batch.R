# batch.R

#' Build a manifest of slide files for batch processing
#'
#' Scans a directory (or takes an explicit list of paths) for slide files
#' [fq_read()] can open and returns one row per file. The manifest fixes the set
#' of slides for a batch and gives each a stable `slide_id` for joining scores
#' back to their file. Add experimental covariates (treatment, sex, age, ...) as
#' further columns; [fq_run()] propagates them to the results.
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

#' Score a folder of slides against one analyzer
#'
#' Runs the full pipeline over a manifest: fits the analyzer once on a
#' representative subsample of slides (or reuses a fit you pass in), then splits
#' and scores every slide, returning one row per section with the manifest's
#' covariates joined on. Maps are not produced here -- regenerate any you want
#' with [fq_render()] from the returned fit and a section.
#'
#' Parallelism is opt-in and entirely the caller's: set `mirai::daemons(n)` in
#' your session and the per-slide work runs across those daemons; with no
#' daemons set it runs sequentially. A real parallel run needs the package
#' installed (`R CMD INSTALL`), as the daemons load `fibroquant` to score.
#'
#' @param manifest A manifest tibble from [fq_manifest()] (columns `slide_id`,
#'   `path`, plus any covariates).
#' @param analyzer An analyzer spec (e.g. [fq_kmeans()]) to fit on the
#'   subsample, or an already-fitted [fq_analyzer] to apply as-is.
#' @param n Number of tissue sections to keep per slide, passed to [fq_split()].
#' @param close_um Morphological closing radius in microns for [fq_split()];
#'   bridges alveolar airspace within a section without bridging the gap between
#'   sections.
#' @param min_area_frac Minimum connected-component area for [fq_split()], as a
#'   fraction of the largest; smaller components (debris, scan streaks) are
#'   dropped.
#' @param target_um_px Working resolution in microns/pixel for [fq_read()]; the
#'   nearest level is chosen. Keep it matched between the fit and the scored
#'   slides so the basis and the data it scores share a resolution.
#' @param n_ref Number of slides to subsample for the fit, or `NULL` (default) to
#'   use every slide. Ignored when `analyzer` is already fitted. The first
#'   section of each reference slide is pooled, so near-duplicate sections from
#'   one slide are not double-counted.
#' @param stratify Optional manifest column to balance the subsample across,
#'   e.g. `"treatment"` so the fit spans the injury range.
#' @param seed Seed for the fit's slide subsample and k-means restarts, for a
#'   reproducible basis. Scoring is deterministic and needs no seed.
#' @return A list with `scores` (one row per section: `slide_id`, `section`,
#'   `n_sections`, the analyzer's metrics, and the joined covariates) and `fit`
#'   (the fitted [fq_analyzer], for reuse or rendering).
#' @export
fq_run <- function(
    manifest,
    analyzer,
    n = 2,
    close_um = 50,
    min_area_frac = 0.05,
    target_um_px = 4,
    n_ref = NULL,
    stratify = NULL,
    seed = 1
) {
  if (!is.data.frame(manifest)) {
    stop("`manifest` must be a data frame from fq_manifest().", call. = FALSE)
  }
  needed <- setdiff(c("slide_id", "path"), names(manifest))
  if (length(needed) > 0) {
    stop(
      "manifest is missing column(s): ",
      paste(needed, collapse = ", "),
      call. = FALSE
    )
  }

  fit <-
    if (S7::S7_inherits(analyzer, fq_analyzer)) {
      analyzer
    } else if (S7::S7_inherits(analyzer, fq_spec)) {
      .fit_reference(
        manifest = manifest,
        spec = analyzer,
        n = n,
        close_um = close_um,
        min_area_frac = min_area_frac,
        target_um_px = target_um_px,
        n_ref = n_ref,
        stratify = stratify,
        seed = seed
      )
    } else {
      stop(
        "`analyzer` must be an fq_spec or a fitted fq_analyzer.",
        call. = FALSE
      )
    }

  # One self-contained worker per slide. in_parallel() crates it (via carrier)
  # so it serialises to mirai daemons; it references only installed-package
  # functions plus the captured fit and split parameters. With no daemons set,
  # purrr runs it sequentially.
  worker <-
    purrr::in_parallel(
      function(row) {
        slide <- fibroquant::fq_read(row$path, target_um_px = target_um_px)
        sections <-
          fibroquant::fq_split(
            slide,
            n = n,
            close_um = close_um,
            min_area_frac = min_area_frac
          )
        n_sections <- length(sections)
        graded <-
          lapply(
            seq_len(n_sections),
            function(i) {
              section <- sections[[i]]
              keys <-
                tibble::tibble(
                  slide_id = row$slide_id,
                  section = section@section,
                  n_sections = n_sections
                )
              dplyr::bind_cols(keys, fibroquant::fq_score(fit, section))
            }
          )
        dplyr::bind_rows(graded)
      },
      fit = fit,
      n = n,
      close_um = close_um,
      min_area_frac = min_area_frac,
      target_um_px = target_um_px
    )

  rows <-
    lapply(
      seq_len(nrow(manifest)),
      function(i) {
        list(
          slide_id = manifest$slide_id[i],
          path = manifest$path[i]
        )
      }
    )

  # Parallelise only if the caller has set daemons; otherwise run in-process.
  parts <- purrr::map(rows, worker)
  scores <- dplyr::bind_rows(parts)

  covariates <- manifest[, setdiff(names(manifest), "path"), drop = FALSE]
  result <-
    dplyr::left_join(
      scores,
      covariates,
      by = "slide_id"
    )

  list(
    scores = result,
    fit = fit
  )
}

# Fit the analyzer on a representative subsample: pick the reference slides,
# read and split each, pool their first sections, and fit. Serial and seeded so
# the basis is reproducible (the seed also governs any per-section pixel cap).
.fit_reference <- function(
    manifest,
    spec,
    n,
    close_um,
    min_area_frac,
    target_um_px,
    n_ref,
    stratify,
    seed
) {
  set.seed(seed)
  ref <- .sample_reference(manifest, n_ref, stratify)

  ref_sections <-
    lapply(
      seq_len(nrow(ref)),
      function(i) {
        slide <- fq_read(ref$path[i], target_um_px = target_um_px)
        sections <-
          fq_split(
            slide,
            n = n,
            close_um = close_um,
            min_area_frac = min_area_frac
          )
        sections[[1]]
      }
    )

  fq_fit(spec, ref_sections)
}

# Choose the reference slides for the fit: all of them when n_ref covers the
# set, else a random draw, or an even draw across the strata of one column.
.sample_reference <- function(manifest, n_ref, stratify) {
  n_slides <- nrow(manifest)
  if (is.null(n_ref) || n_ref >= n_slides) {
    return(manifest)
  }
  if (is.null(stratify)) {
    rows <- sort(sample.int(n_slides, n_ref))
    return(manifest[rows, , drop = FALSE])
  }
  if (!stratify %in% names(manifest)) {
    stop("`stratify` column not found in manifest: ", stratify, call. = FALSE)
  }

  groups <- split(seq_len(n_slides), manifest[[stratify]])
  per <- max(1L, floor(n_ref / length(groups)))
  picked <-
    unlist(
      lapply(
        groups,
        function(g) if (length(g) <= per) g else sample(g, per)
      ),
      use.names = FALSE
    )
  manifest[sort(picked), , drop = FALSE]
}
