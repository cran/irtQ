
<!-- README.md is generated from README.Rmd. Please edit that file -->

# irtQ 0.1.1

- Resolved the misalignment issue of standard errors in the output of
  the `est_irt()` function when `fix.a.1pl = TRUE` is specified and the
  items are calibrated using the 1PLM.

- Added a new function, `grdif()`, to perform differential item
  functioning (DIF) analysis across multiple groups. This function
  calculates three generalized IRT residual DIF (GRDIF) statistics. For
  more information about the function and its usage, please refer to the
  accompanying documentation.

- Fixed several typos in the manual documentation

# irtQ 0.1.0

- Initial release on CRAN

- The `irtQ` package is a successor of the `irtplay` package which was
  retracted from R CRAN due to the intellectual property (IP) violation.
  All issues of the IP violation have been clearly resolved in the
  `irtQ` package.

- Most of the functions the `irtQ` package are identical in appearance
  and functionality to those of `irtplay` package except a few functions
  (e.g., `shape_df()`, `est_score()`). However, the computing speed of
  several functions (e.g., `est_irt()`, `est_score()`, `lwrc()`) in the
  `irtQ` package are faster than the previous ones in the `irtplay`
  package. Read the documentation carefully prior to using the
  functions.
