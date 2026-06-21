# batch.R

# fq_manifest ------------------------------------------------------------------

#' Build a manifest of slide files for batch processing
#'
#' Scans a directory, or takes an explicit vector of paths, for slide files that
#' [fq_read()] can open, and returns one row per file. The manifest fixes the
#' set of slides for a batch and gives each a stable `slide_id` for joining
#' scores back to their file. Add experimental covariates such as treatment,
#' sex, or age as further columns, and [fq_run()] carries them into the results.
#'
#' @param path A directory to scan, or a character vector of slide file paths.
#' @param recursive Recurse into sub-directories when `path` is a directory.
#' @param extensions File extensions to match, in lower case, with no leading
#'   dot.
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

# fq_run -----------------------------------------------------------------------

#' Score a folder of slides against one analyzer
#'
#' Runs the full pipeline over a manifest. It fits the analyzer once on a
#' representative subsample of slides, or reuses a fit you pass in, then splits
#' and scores every slide. It returns one row per section with the manifest's
#' covariates joined on. A slide that fails to read or score is skipped with a
#' warning rather than stopping the run. It does not produce severity maps.
#' Generate any you need with [fq_render()] from the returned fit and a section.
#'
#' Parallelism is optional and controlled by the caller. Set `mirai::daemons(n)`
#' in your session and the work for each slide runs across those daemons. With
#' no daemons set it runs sequentially. A parallel run needs the package
#' installed with `R CMD INSTALL`, because each daemon loads `fibroquant` to
#' score.
#'
#' @param manifest A manifest tibble from [fq_manifest()] (columns `slide_id`,
#'   `path`, plus any covariates).
#' @param analyzer An analyzer spec such as [fq_kmeans()] to fit on the
#'   subsample, or an already fitted [fq_analyzer] to apply directly.
#' @param n Number of tissue sections to keep per slide, passed to [fq_split()].
#' @param close_um Morphological closing radius in microns, passed to
#'   [fq_split()]. It bridges alveolar airspace within a section without
#'   bridging the gap between sections.
#' @param min_area_frac Minimum area of a connected component for [fq_split()],
#'   as a fraction of the largest. Smaller components such as debris and scan
#'   streaks are dropped.
#' @param target_um_px Working resolution in microns per pixel for [fq_read()].
#'   The nearest level is chosen. Keep this matched between the fit and the
#'   scored slides so the basis and the data it scores share a resolution.
#' @param n_ref Total number of slides to subsample for the fit (default 25), or
#'   `Inf`/`NULL` to use every slide. This budget bounds the read-and-split cost
#'   of fitting. Ignored when `analyzer` is already fitted. Only the first
#'   section of each reference slide is pooled, so similar sections from one
#'   slide are not counted twice.
#' @param stratify Optional manifest column to balance the subsample across,
#'   such as `"treatment"`, so the fit spans the range of injury. The `n_ref`
#'   budget is spread as evenly as possible across the column's groups. When
#'   there are more groups than `n_ref`, a random `n_ref` of them are used, one
#'   slide each.
#' @param seed Seed for the fit's slide subsample and k-means restarts, for a
#'   reproducible basis. Scoring is deterministic and needs no seed.
#' @param progress Show a progress bar that advances as each slide finishes, in
#'   both sequential and parallel runs. `FALSE` suppresses it.
#' @return A list with two elements. `scores` has one row per section, with
#'   columns `slide_id`, `section` (the section label), the analyzer's metrics
#'   from [fq_score()], and the manifest's covariates joined on. `fit` is the
#'   fitted [fq_analyzer], for reuse or rendering.
#' @export
fq_run <- function(
    manifest,
    analyzer,
    n = 2,
    close_um = 50,
    min_area_frac = 0.05,
    target_um_px = 4,
    n_ref = 25,
    stratify = NULL,
    seed = 1,
    progress = TRUE
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

  # One worker per slide. in_parallel() crates it via carrier so it serialises
  # to mirai daemons. It references only installed package functions plus the
  # captured fit and split parameters. With no daemons set, purrr runs it
  # sequentially. A slide that fails to read or score returns a row carrying its
  # error instead of stopping the batch.
  worker <-
    purrr::in_parallel(
      function(row) {
        tryCatch(
          {
            slide <- fibroquant::fq_read(row$path, target_um_px = target_um_px)
            sections <-
              fibroquant::fq_split(
                slide,
                n = n,
                close_um = close_um,
                min_area_frac = min_area_frac
              )
            graded <-
              lapply(
                seq_along(sections),
                function(i) {
                  section <- sections[[i]]
                  keys <-
                    tibble::tibble(
                      slide_id = row$slide_id,
                      section = section@section
                    )
                  dplyr::bind_cols(keys, fibroquant::fq_score(fit, section))
                }
              )
            dplyr::bind_rows(graded)
          },
          error = function(e) {
            tibble::tibble(
              slide_id = row$slide_id,
              .error = conditionMessage(e)
            )
          }
        )
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

  parts <- purrr::map(rows, worker, .progress = progress)
  scores <- dplyr::bind_rows(parts)

  # Failed slides carry a .error column and no metrics. Warn and drop them so
  # one bad slide does not sink the whole batch.
  if (".error" %in% names(scores)) {
    failed <- scores[!is.na(scores$.error), , drop = FALSE]
    if (nrow(failed) > 0) {
      warning(
        "Skipped ", nrow(failed), " slide(s) that failed to score: ",
        paste0(failed$slide_id, " (", failed$.error, ")", collapse = "; "),
        call. = FALSE
      )
    }
    scores <- scores[
      is.na(scores$.error),
      setdiff(names(scores), ".error"),
      drop = FALSE
    ]
  }

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

# Internal helpers -------------------------------------------------------------

# Fit the analyzer on a representative subsample. Pick the reference slides,
# read and split each, pool their first sections, and fit. The seed governs the
# subsample, the k-means restarts, and any pixel cap. with_seed scopes it so the
# caller's RNG is left untouched.
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
  withr::with_seed(seed, {
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
  })
}

# Choose the reference slides for the fit. Use all of them when the budget
# covers the set. Otherwise draw n_ref at random, or spread n_ref as evenly as
# possible across the strata of one column, capping the total either way.
.sample_reference <- function(manifest, n_ref, stratify) {
  n_slides <- nrow(manifest)
  if (is.null(n_ref) || n_ref >= n_slides) {
    return(manifest)
  }
  if (is.null(stratify)) {
    return(manifest[sort(sample.int(n_slides, n_ref)), , drop = FALSE])
  }
  if (!stratify %in% names(manifest)) {
    stop("`stratify` column not found in manifest: ", stratify, call. = FALSE)
  }

  groups <- split(seq_len(n_slides), manifest[[stratify]])
  if (n_ref < length(groups)) {
    message(
      "n_ref (", n_ref, ") is fewer than the ", length(groups),
      " strata; fitting on a random ", n_ref, " of them."
    )
  }
  alloc <- .allocate_budget(lengths(groups), n_ref)
  picked <-
    unlist(
      Map(
        function(rows, k) if (k >= length(rows)) rows else sample(rows, k),
        groups,
        alloc
      ),
      use.names = FALSE
    )
  manifest[sort(picked), , drop = FALSE]
}

# Spread a total budget across groups as evenly as their sizes allow. Fill the
# least-filled groups first, breaking ties at random, until the budget is spent.
# When the budget is smaller than the group count, a random subset of groups
# each get one. Assumes budget is less than the total of sizes.
.allocate_budget <- function(sizes, budget) {
  alloc <- integer(length(sizes))
  while (budget > 0) {
    open <- which(alloc < sizes)
    least <- open[alloc[open] == min(alloc[open])]
    if (length(least) > 1) {
      least <- sample(least)
    }
    take <- min(length(least), budget)
    winners <- least[seq_len(take)]
    alloc[winners] <- alloc[winners] + 1L
    budget <- budget - take
  }
  alloc
}
