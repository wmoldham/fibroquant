# .Rprofile

suppressPackageStartupMessages({
  library(usethis)
  library(devtools)
})

Sys.setenv(
  FIBROQUANT_TEST_VSI =
    switch(
      Sys.info()[["sysname"]],
      Darwin  = "/Volumes/Will/Mouse lung 6.10.26/Image_3470.vsi",
      Linux   = "/run/media/will/Will/Mouse lung 6.10.26/Image_3470.vsi",
      ""
    )
)
