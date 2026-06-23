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

# fq_results -------------------------------------------------------------------

#' The result of a scoring run
#'
#' The object [fq_run()] returns: the per-section scores together with
#' everything needed to redraw any section's severity map after the fact. It is
#' a typed container over plain tibbles, so `@scores` stays a `dplyr`-able table
#' for the downstream stats join while the run carries its fit and the recipe
#' that produced it.
#'
#' Because the read-and-split path is deterministic, the run does not cache any
#' pixels. The `fit` plus the split `params` and the slide `paths` are the
#' closed set needed to replay a section and render its map. [fq_map()] does
#' that replay in one call.
#'
#' @param scores A tibble with one row per section: `slide_id`, `section`, the
#'   analyzer's metrics from [fq_score()], and the manifest's covariates.
#' @param fit The fitted [fq_analyzer] the run scored with, for reuse or to
#'   render maps. Typed to the [fq_analyzer] superclass so any analyzer fits.
#' @param params The split recipe that produced the scored sections, a named
#'   list with `n`, `close_um`, `min_area_frac`, and `target_um_px`.
#' @param paths A two-column tibble mapping each `slide_id` to its `path`, so a
#'   scored section can be replayed from its source file.
#' @return An `fq_results` object.
#' @seealso [fq_run()], [fq_map()]
#' @export
fq_results <-
  S7::new_class(
    "fq_results",
    properties = list(
      scores = S7::class_data.frame,
      fit = fq_analyzer,
      params = S7::class_list,
      paths = S7::class_data.frame
    ),
    validator = function(self) {
      recipe <- c("n", "close_um", "min_area_frac", "target_um_px")
      missing_recipe <- setdiff(recipe, names(self@params))
      if (length(missing_recipe) > 0L) {
        return(
          paste0(
            "@params must carry the split recipe; missing: ",
            paste(missing_recipe, collapse = ", ")
          )
        )
      }
      if (!all(c("slide_id", "path") %in% names(self@paths))) {
        return("@paths must have columns `slide_id` and `path`")
      }
      NULL
    }
  )

# A run's analyzer label for printing: the spec's class name (e.g. "fq_kmeans"),
# falling back to the fit class with its "_analyzer" suffix stripped when a fit
# carries no spec.
# TODO: the no-spec fallback branch is untested while fq_kmeans_analyzer (which
# always carries a spec) is the only analyzer. Add a test once an analyzer that
# omits @spec lands (cpa/density/texture).
.analyzer_label <- function(fit) {
  spec <- tryCatch(fit@spec, error = function(e) NULL)
  if (!is.null(spec)) {
    S7::S7_class(spec)@name
  } else {
    sub("_analyzer$", "", S7::S7_class(fit)@name)
  }
}

S7::method(print, fq_results) <- function(x, ...) {
  n_sections <- nrow(x@scores)
  n_slides <- nrow(x@paths)
  cat(
    sprintf(
      "%d %s from %d %s \u00b7 %s \u00b7 %g \u00b5m/px\n",
      n_sections, if (n_sections == 1L) "section" else "sections",
      n_slides, if (n_slides == 1L) "slide" else "slides",
      .analyzer_label(x@fit),
      x@params$target_um_px
    )
  )
  invisible(x)
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
#' @return An [fq_results] object. Its `@scores` slot has one row per section,
#'   with columns `slide_id`, `section` (the section label), the analyzer's
#'   metrics from [fq_score()], and the manifest's covariates joined on. `@fit`
#'   is the fitted [fq_analyzer], for reuse or rendering. `@params` is the split
#'   recipe (`n`, `close_um`, `min_area_frac`, `target_um_px`) and `@paths` maps
#'   each `slide_id` to its file. Together `@fit`, `@params`, and `@paths` are
#'   everything [fq_map()] needs to redraw any section's severity map after the
#'   run.
#' @seealso [fq_map()]
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

  fq_results(
    scores = result,
    fit = fit,
    params = list(
      n = n,
      close_um = close_um,
      min_area_frac = min_area_frac,
      target_um_px = target_um_px
    ),
    paths = tibble::tibble(
      slide_id = manifest$slide_id,
      path = manifest$path
    )
  )
}

# fq_map -----------------------------------------------------------------------

#' Redraw a scored section's severity map after a run
#'
#' [fq_run()] keeps the scores and the fit, not the sections themselves, so a
#' row in the scores table is not directly tied to its pixels. `fq_map()` closes
#' that gap for post-hoc review: it replays [fq_read()] and [fq_split()] with
#' the recipe the run recorded, then renders the section through the run's fit,
#' returning a plottable [fq_field] in one call. Because the read-and-split path
#' is deterministic, nothing is cached between runs; the map is recomputed on
#' demand from the run's `@params` and `@paths`.
#'
#' `run[[slide_id, section]]` is shorthand for `fq_map(run, slide_id, section)`,
#' and `run[[slide_id]]` for `fq_map(run, slide_id)`.
#'
#' @param run An [fq_results] from [fq_run()].
#' @param slide_id A `slide_id` from `run@scores` (or `run@paths`).
#' @param section A section label such as `"A"` or `"B"`, as it appears in the
#'   `section` column of `run@scores`. `NULL` returns a named list of maps, one
#'   per section of the slide.
#' @return An [fq_field], or a named list of `fq_field`s when `section` is
#'   `NULL`.
#' @seealso [fq_run()], [fq_render()]
#' @export
fq_map <- function(run, slide_id, section = NULL) {
  if (!S7::S7_inherits(run, fq_results)) {
    stop("`run` must be an fq_results object from fq_run().", call. = FALSE)
  }
  if (is.null(section)) {
    sections <- .fq_replay(run, slide_id)
    maps <- lapply(sections, function(s) fq_render(run@fit, s))
    names(maps) <- vapply(sections, function(s) s@section, character(1))
    return(maps)
  }
  # Force the replay before dispatch. Passing it as a lazy argument into the
  # fq_render generic lets S7 force the section promise mid-dispatch, turning a
  # clean "not found" stop() into a re-entrant promise error.
  sec <- .fq_replay(run, slide_id, section)
  fq_render(run@fit, sec)
}

# run[["Image_3470", "A"]] and run[["Image_3470"]] delegate to fq_map(); the
# terse idiom for interactive review lands on the same code path as the verb.
S7::method(`[[`, fq_results) <- function(x, i, j, ...) {
  if (missing(j)) fq_map(x, i) else fq_map(x, i, j)
}

# Internal helpers -------------------------------------------------------------

# Replay read -> split -> pick for one slide of a run. The split is
# deterministic, so this rebuilds the exact section(s) the run scored from the
# recipe in @params. Returns the fq_sections collection when section is NULL,
# otherwise the single matching fq_section.
.fq_replay <- function(run, slide_id, section = NULL) {
  i <- match(slide_id, run@paths$slide_id)
  if (is.na(i)) {
    stop("slide_id not found in run: ", slide_id, call. = FALSE)
  }

  p <- run@params
  slide <- fq_read(run@paths$path[i], target_um_px = p$target_um_px)
  sections <-
    fq_split(
      slide,
      n = p$n,
      close_um = p$close_um,
      min_area_frac = p$min_area_frac
    )
  if (is.null(section)) {
    return(sections)
  }

  labels <- vapply(sections, function(s) s@section, character(1))
  j <- match(section, labels)
  if (is.na(j)) {
    stop(
      "section '", section, "' not found for slide ", slide_id,
      "; available: ", paste(labels, collapse = ", "),
      call. = FALSE
    )
  }
  sections[[j]]
}

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
