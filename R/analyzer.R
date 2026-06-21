# analyzer.R

# Abstract classes -------------------------------------------------------------

#' Analyzer specification
#'
#' Abstract parent of the analyzer specs. A spec records which analyzer to run
#' and its hyperparameters. It holds no fitted state. Construct a concrete spec
#' with a constructor such as `fq_kmeans()`.
#'
#' @export
fq_spec <-
  S7::new_class(
    "fq_spec",
    abstract = TRUE
  )

#' Fitted analyzer
#'
#' Abstract parent of the fitted analyzer objects returned by [fq_fit()]. A fit
#' pairs a spec with the basis learned from a batch of sections, such as cluster
#' centres, a severity ordering, and thresholds. It is the input to [fq_score()]
#' and [fq_render()].
#'
#' @export
fq_analyzer <-
  S7::new_class(
    "fq_analyzer",
    abstract = TRUE
  )

# Generics ---------------------------------------------------------------------

#' Fit an analyzer on a batch of sections
#'
#' Fits a [spec][fq_spec] once on the pooled sections and learns a basis that is
#' reused to score every section in the batch. Dispatches on `spec`.
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
#' Applies a [fit][fq_analyzer] to one section and returns its metrics as a
#' single row. Dispatches on `fit`.
#'
#' Every method must return a one-row tibble with a numeric `severity_index`
#' column. That shared column lets scores from any analyzer be stacked into one
#' table and is what [fq_run()] joins covariates onto. Methods may add any
#' number of further columns.
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
#' Applies a [fit][fq_analyzer] to one section and returns a severity grade for
#' every pixel, ready for pseudocolouring. Dispatches on `fit`.
#'
#' @param fit An [fq_analyzer] from [fq_fit()].
#' @param section An `fq_section`.
#' @param ... Passed on to methods.
#' @return A severity field with one grade per pixel, `NA` off tissue.
#' @export
fq_render <-
  S7::new_generic("fq_render", "fit", function(fit, section, ...) {
    S7::S7_dispatch()
  })
