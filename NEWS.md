# UPDhmm 0.99.0

First release of UPDhmm package

# UPDhmm 0.99.1

First pre-review corrections

# UPDhmm 0.99.2

Correction of reviewer's suggestions

# UPDhmm 0.99.3

2nd correction of reviewer's suggestions

# UPDhmm 0.99.4

Final corrections and final editing of vignette

# UPDhmm 1.5.0

- Added parallelization option `BPPARAM` in `calculateEvents()` and some new columns in the output.
- Added new vignette: "Preprocessing Guide".
- Update User Guide vignette
- Add 3 new functions `markRecurentRegions()`, `identifyRecurrentRegions()` and `collapseEvents()`

# UPDhmm 1.5.1
- Update vignette

# UPDhmm 1.5.2
- Fix vignette for LatEX error.

# UPDhmm 1.7.1

NEW FEATURES

- Added a long-read processing example in a new vignette.

IMPROVEMENTS

- Refined recurrent annotation in markRecurrentRegions().
Recurrent status is now determined using a minimum overlap fraction between events and recurrent regions, preventing marginal overlaps from being annotated as recurrent.

- More precise recurrent region definition in identifyRecurrentRegions().
Replaced direct interval merging with an overlap-aware, chromosome-wise agglomerative hierarchical clustering approach. Recurrent regions are now defined based on overlap similarity between events, with configurable clustering stringency.

BUG FIXES

- Corrected collapsed event size calculation in collapseEvents().
The total size of collapsed events is now computed as the true genomic span of the region rather than the sum of individual event sizes.
