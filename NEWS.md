# metalcor 0.0.0.9000 (2026-05-05)

- First commit, includes all functions well documented with examples, thoroughly tested, and exported as needed; also description, readme, and this news; passes check.

# metalcor 0.0.1.9000 (2026-05-05)

- Function `pprodcor` corrected low accuracy of calculations for `abs(rho)` values between 0.95 and 1, which gets particularly worse as 1 is approached, and caused large errors at 0.999, for example.  Correction solves problems upstream, in `qprodcor` and `rho_from_median`.
- Added vignette

# metalcor 1.0.0 (2026-05-06)

- CRAN submission
  - Minor spellcheck edits
- Added logo, displayed on readme

# metalcor 1.0.1 (2026-08-12)

- Added accepted paper citation to README, updated DOI in DESCRIPTION.

# metalcor 1.0.2 (2026-09-14)

- Functions `*prodcor` fixed a bug revealed by a new warning in R-devel only: "object length is not a multiple of subscript length".  The bug was there previously, though it occurred rarely, but it was silent in previous R versions.
  - From R-devel news, this is the change that triggers the new warnings addressed by this fix: "x[logi], with logi a logical index, now warns if length(x) is not a multiple of length(logi), fixing PR#19155 thanks to several participants in two R Dev Days and the R Sprint in 2026."
