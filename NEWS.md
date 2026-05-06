# metalcor 0.0.0.9000 (2026-05-05)

- First commit, includes all functions well documented with examples, thoroughly tested, and exported as needed; also description, readme, and this news; passes check.

# metalcor 0.0.1.9000 (2026-05-05)

- Function `pprodcor` corrected low accuracy of calculations for `abs(rho)` values between 0.95 and 1, which gets particularly worse as 1 is approached, and caused large errors at 0.999, for example.  Correction solves problems upstream, in `qprodcor` and `rho_from_median`.
- Added vignette
