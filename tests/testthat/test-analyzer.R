# test-analyzer.R

# A throwaway analyzer that exercises the fit / score / render contract without
# a real backend.
demo_spec <- S7::new_class("demo_spec", parent = fq_spec)
demo_fit <- S7::new_class(
  "demo_fit",
  parent = fq_result,
  properties = list(value = S7::class_numeric)
)

S7::method(fq_fit, demo_spec) <- function(spec, sections, ...) {
  demo_fit(value = length(sections))
}
S7::method(fq_score, demo_fit) <- function(fit, section, ...) {
  tibble::tibble(severity_index = fit@value)
}
S7::method(fq_render, demo_fit) <- function(fit, section, ...) {
  matrix(fit@value, nrow = 1, ncol = 1)
}

test_that("abstract specs and fits cannot be constructed", {
  expect_error(fq_spec())
  expect_error(fq_result())
})

test_that("fq_fit dispatches on the spec and returns a fit", {
  fit <- fq_fit(demo_spec(), list(1, 2, 3))
  expect_true(S7::S7_inherits(fit, fq_result))
  expect_true(S7::S7_inherits(fit, demo_fit))
  expect_equal(fit@value, 3)
})

test_that("fq_score and fq_render dispatch on the fit", {
  fit <- fq_fit(demo_spec(), list(1, 2))
  expect_equal(fq_score(fit, NULL)$severity_index, 2)
  expect_equal(fq_render(fit, NULL)[1, 1], 2)
})

test_that("a generic errors when no method matches", {
  expect_error(fq_fit("not a spec", list()))
})
