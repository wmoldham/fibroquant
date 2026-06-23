# fibroquant 0.1.0

First release: an end-to-end pipeline for quantifying lung fibrosis from
whole-slide histology images, built around a pluggable analyzer interface.

## Pipeline

* `fq_read()` loads a whole-slide image at a chosen working resolution, and
  `fq_split()` separates a multi-section slide into per-section `fq_section`s by
  connected components.
* `fq_manifest()` builds a stable slide index for batch work, carrying
  experimental covariates through to the results.

## Analyzers

* A pluggable analyzer interface — `fq_fit()`, `fq_score()`, and `fq_render()`
  over the `fq_spec`/`fq_analyzer` classes — so new analyzers drop in beside the
  first one.
* `fq_kmeans()`: the first analyzer, unsupervised CIELAB colour clustering of
  tissue ranked by lightness into severity grades.

## Running and reviewing

* `fq_run()` scores a folder of slides against one analyzer and returns an
  `fq_results` object (`@scores`, `@fit`, `@params`, `@paths`) rather than a bare
  list.
* `fq_map(run, slide_id, section)` — and the `run[[slide_id, section]]`
  shorthand — redraws any scored section's severity map after a run, replaying
  read-and-split from the recipe the run recorded.
* `plot()` methods for slides, section collections, and severity fields.
