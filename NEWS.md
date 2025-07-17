
<!-- README.md is generated from README.Rmd. Please edit that file -->

# irtQ 1.0.0

- The documentation for the `irtQ` package has been revised to reflect
  updates to function behavior, fix typos, and provide more relevant and
  detailed information for existing functions.

- A new function, `crdif()`, has been added. This function computes
  three statistics from the residual-based DIF detection framework using
  categorical residuals (RDIF-CR). It allows for the detection of global
  DIF, particularly in polytomously scored items.

- A new function, `shape_df_fipc()`, has been introduced. This function
  merges fixed-item metadata with automatically generated metadata for
  new items and produces a single data frame ordered by test position.
  It is designed to support fixed item parameter calibration (FIPC) via
  the `est_irt()` function.

- The `plot()` method has been enhanced to support the display of all
  item characteristic curves for a given item in a single panel.

- The `pcd2()` function has been updated to include a purification
  procedure.

- The `rdif()` and `catsib()` functions now include an `item.skip`
  argument. This allows users to specify a numeric vector of item
  indices to exclude from the DIF analysis.

# irtQ 0.2.1

- Enhanced functionality of the `bind.fill()` function by adding a new
  argument `fill`. The value in the argument is used to fill in missing
  data when aligning datasets.

- Fixed a bug within the `est_irt()` function that was previously unable
  to implement the fixed item parameter calibration (FIPC) when only
  freely estimating a single item given that all other items are fixed.

- Added a new function, `reval_mst()`, which evaluates the measurement
  precision and bias in Multistage-adaptive Test (MST) panels using a
  recursion-based evaluation method introduced by Lim et al. (2020).

- Added a new function, `pcd2()`, which computes the Pseudo-count
  $`D^{2}`$ statistics (Cappaert et al., 2018; Stone, 2000) to detect
  item parameter drift.

# irtQ 0.2.0

- Introduced Warm’s (1989) Weighted Likelihood (WL) estimation method to
  the `est_score()` function. This WL scoring method can now be utilized
  by setting `method = "WL"`.

- Enhanced the speed of ability parameter estimation in the
  `est_score()` function when using the ML, MLF, or MAP methods for the
  `method` argument. The updated version performs approximately three
  times faster than its predecessor.

- Addressed a bug within the `est_score()` function that was previously
  unable to accurately compute scores when only a single item data was
  provided. This issue was occurring with the EAP.SUM and INV.TCC
  estimation methods.

- Added two new functions for computing classification accuracy and
  consistency: `cac_rud()` and `cac_lee()`.

  - `cac_rud`: This function implements Rudner’s (2001, 2005) method for
    computing classification accuracy and consistency. It takes cut
    scores, ability estimates, standard errors, and optional weights as
    inputs and returns a list containing a confusion matrix, marginal
    and conditional classification accuracy and consistency indices, the
    probability of being assigned to each level category, and the cut
    scores used in the analysis.
  - `cac_lee`: This function implements Lee’s (2010) method for
    computing classification accuracy and consistency. It takes a data
    frame containing item metadata, cut scores, optional ability
    estimates, optional weights, a scaling factor, and a logical value
    indicating the cut score metric as inputs. It returns a list similar
    to `cac_rud`.

- Added a new function, `llike_score()`, which computes the
  loglikelihood of ability parameters given the item parameters and
  response data.

- Enhanced functionality of the `rdif()` and `grdif()` functions: Both
  now support the graded response model (GRM) and generalized partial
  credit model (GPCM).

- Fixed an issue in the `grdif()` function that inaccurately calculated
  the GRDIF statistics when group membership was specified in a
  non-standard way. Specifically, the problem arose when 0 wasn’t used
  as the reference group and consecutive numbers (e.g., 1, 2, 3) weren’t
  used to represent focal groups in the `group` argument.

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
