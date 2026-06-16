# helper-vsi.R

# Path to a local .vsi for integration tests, or "" when unset.
vsi_path <- function() {
  Sys.getenv("FIBROQUANT_TEST_VSI")
}

skip_if_no_vsi <- function() {
  if (vsi_path() == "") {
    testthat::skip("Set FIBROQUANT_TEST_VSI to a .vsi path to run integration tests.")
  }
}
