# analyzer.R

#' Analyzer specification
#'
#' Abstract parent of the analyzer specs. A spec is the recipe -- which analyzer
#' to run and its hyperparameters -- and holds no fitted state. Construct a
#' concrete spec with a family constructor such as `fq_kmeans()`.
#'
#' @export
fq_spec <-
  S7::new_class(
    "fq_spec",
    abstract = TRUE
  )

#' Fitted analyzer
#'
#' Abstract parent of the fitted-analyzer objects returned by [fq_fit()]. A fit
#' pairs a spec with the basis learned from a batch of sections -- cluster
#' centres, severity ordering, thresholds -- and is the input to [fq_score()]
#' and [fq_render()].
#'
#' @export
fq_analyzer <-
  S7::new_class(
    "fq_analyzer",
    abstract = TRUE
  )

#' Fit an analyzer on a batch of sections
#'
#' Fits a [spec][fq_spec] once on the pooled sections and learns the shared
#' basis reused across the batch (fit-once / apply-many). Dispatches on `spec`.
#'
#' @param spec An analyzer spec, e.g. `fq_kmeans()`.
#' @param sections A list of `fq_section`s (or an `fq_sections`).
#' @param ... Passed on to methods.
#' @return An [fq_analyzer].
#' @export
fq_fit <-
  S7::new_generic("fq_fit", "spec", function(spec, sections, ...) {
    S7::S7_dispatch()
  })

#' Score a section against a fitted analyzer
#'
#' Applies a [fit][fq_analyzer] to one section and returns its per-section
#' metrics. Dispatches on `fit`.
#'
#' @param fit An [fq_analyzer] from [fq_fit()].
#' @param section An `fq_section`.
#' @param ... Passed on to methods.
#' @return A one-row tibble: `severity_index` plus analyzer-specific columns.
#' @export
fq_score <-
  S7::new_generic("fq_score", "fit", function(fit, section, ...) {
    S7::S7_dispatch()
  })

#' Render a section's severity field
#'
#' Applies a [fit][fq_analyzer] to one section and returns its per-pixel severity
#' field for pseudocolouring. Dispatches on `fit`.
#'
#' @param fit An [fq_analyzer] from [fq_fit()].
#' @param section An `fq_section`.
#' @param ... Passed on to methods.
#' @return An H x W severity field, `NA` off tissue.
#' @export
fq_render <-
  S7::new_generic("fq_render", "fit", function(fit, section, ...) {
    S7::S7_dispatch()
  })
