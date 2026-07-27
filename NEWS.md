
<!-- README.md is generated from README.Rmd. Please edit that file -->

# irtQ 1.2.0

## New Features

- Added a new function, `find_cut()`, which identifies TIF-crossing
  routing cut scores for MST panels. For each adjacent module pair
  within a stage, the function locates the theta where the two modules'
  test information functions (TIFs) intersect. A warning is issued when
  the mean difficulty order of modules within a stage differs from their
  input index order.

- Added a new S3 method, `plot.find_cut()`, which visualizes TIF curves
  and routing cut scores stage by stage using **ggplot2** facets.
  Proper, anomalous, and unselected cut scores are distinguished by line
  type. A `layout` argument (`"vertical"` / `"horizontal"`) controls the
  facet orientation.

- Added a new function, `run_mst()`, which simulates MST administrations
  for a given panel structure and returns response data along with
  ability and routing information for each simulated examinee.

- Added a new exported dataset, `simIPD`, which contains simulated CAT
  response data for illustrating IPD detection with `ripd()` and
  `pcd2()`. The dataset represents one replication of a CAT simulation
  (N = 3,000; test length = 30; 360-item 3PLM pool) in which 5% of items
  (18 items) had both $a$ and $b$ parameters decreased by 0.5.

## Minor Improvements

- `est_irt()`, `est_item()`, and `est_mg()` now accept a partial
  `control` list. Users can specify only the arguments they wish to
  override (e.g., `control = list(iter.max = 500)`); unspecified
  arguments fall back to their defaults via `modifyList()`.

- Reorganized the **pkgdown** reference page: `ripd()` and `pcd2()` are
  now grouped under a new *Item Parameter Drift (IPD)* section, and
  `reval_mst()`, `panel_info()`, `find_cut()`, `plot.find_cut()`, and
  `run_mst()` are grouped under a new *Multistage-Adaptive Test (MST)*
  section.

## Documentation

- Expanded the documentation for `ripd()` by adding a `@details` section
  that covers the theoretical background of the RIPD framework,
  asymptotic distributions of the three RIPD statistics, a drift-type
  diagnostic guide, the CAT-specific three-step workflow, and the
  purification procedure. Also added a `\donttest{}` example
  demonstrating a complete CAT-based IPD detection workflow using the
  `simIPD` dataset. Updated `@references` to include Lim & Choe (2023)
  and replaced the previous conference paper citation with the in-press
  journal reference (Lim & Han, in press).

- Added a `\donttest{}` example to `pcd2()` demonstrating CAT-based IPD
  detection using the `simIPD` dataset, including the bootstrap critical
  value procedure described in Lim & Han (in press).

- Updated the *MST Panel Evaluation and Simulation* article
  (`vignettes/articles/mst-panel-evaluation.Rmd`) to introduce
  `run_mst()` and extend the existing `reval_mst()` content with routing
  and scoring examples, a `find_cut()`-based principled cut score
  derivation, a side-by-side routing method comparison, and a Monte
  Carlo-vs-analytical validation example (Example 6).

## Bug Fixes

- Rebuilt the `simMST` dataset: the previous version had 9 items
  duplicated across non-adjacent modules because the original assembly
  only enforced no-overlap within a single routing pathway. The new
  version enforces a global no-overlap constraint across all 7 modules
  (56 unique items total) and adds a mean(*b*) per-module band
  constraint; cut scores were regenerated via `find_cut()` on the
  rebuilt modules.

# irtQ 1.1.0

## Major Improvements

- Improved the speed and reduced memory usage of item parameter
  estimation and standard error computation in `est_irt()`,
  `est_item()`, and `est_mg()`.
- Improved the computational speed of `est_score()` by up to 52% for
  dichotomous items (N = 10,000) and up to 27% for mixed-format tests,
  through a series of optimizations.
- Improved the computational speed of `sx2_fit()` substantially by
  replacing the O(J^2) Lord-Wingersky recursion with a forward-backward
  pass (up to 11x faster for mixed-format tests with J = 55 items) and
  vectorizing internal helper functions `expFreq()`, `obsFreq()`, and
  the PRM category-collapsing routine.

## New Features

- Added a unit test suite using the **testthat** 3rd edition
  (`testthat >= 3.0.0`). Tests cover core functions including `drm()`,
  `prm()`, `est_irt()`, `est_score()`, `est_mg()`, `rdif()`, `crdif()`,
  and `catsib()`, with the relevant tests across dichotomous,
  polytomous, and mixed-format item scenarios (355 tests total).
- Added a new function, `ripd()`, which implements the Residual-based
  Item Parameter Drift (RIPD) detection framework. The function computes
  three RIPD statistics: $RIPD_R$, $RIPD_S$, and $RIPD_{RS}$ (one for
  each item). $RIPD_R$ captures uniform item parameter drift (IPD) via
  differences in mean raw residuals between groups, $RIPD_S$ captures
  nonuniform IPD via differences in mean squared residuals, and
  $RIPD_{RS}$ is a combined chi-square-based statistic sensitive to both
  types of drift. An optional purification procedure is also supported.

## New Articles

- Launched the irtQ documentation website at
  <https://hwangQ.github.io/irtQ/>, built with **pkgdown**. The site
  includes a full function reference index and the following vignettes
  covering the complete irtQ workflow:
  - *Getting Started with irtQ*: an end-to-end overview of the package
    workflow.
  - *Item Parameter Estimation*: detailed guidance on `est_irt()`,
    `est_item()`, and `est_mg()`.
  - *Ability Estimation*: scoring methods available in `est_score()`.
  - *Model-Data Fit Evaluation*: using `irtfit()` and `sx2_fit()` to
    assess model fit.
  - *DIF Detection*: applying `rdif()`, `grdif()`, and `catsib()` to
    detect item bias.
  - *Classification Accuracy and Consistency*: computing indices via
    `cac_lee()` and `cac_rud()`.
  - *Utility Functions*: usage of `info()`, `traceline()`, `lwrc()`,
    `simdat()`, and related helpers.
  - *Evaluating MST Panels with `reval_mst()`*: measurement precision
    and bias evaluation for multistage adaptive tests.

## Bug Fixes

- Fixed a minor bug in `sx2_fit()` that caused incorrect cell collapsing
  between two adjacent score categories for polytomous items when
  computing the S-$X^2$ item fit statistic.
- Fixed minor bugs in `est_score()` and `info()`.
- Resolved an issue in `catsib()` where the function failed when all
  responses were missing (NA) in either the reference or focal group.
- Updated `cac_rud()` to include the `x` argument, allowing users to
  pass item metadata data frames directly.
- Revised default `control` parameters in `est_irt()`, `est_item()`, and
  `est_mg()`, and updated the documentation accordingly.
- Fixed a minor bug in `est_score()` function in terms of Newton-Raphson
  method.
- Fixed multiple stability issues in `catsib()`:
  - The final bin exclusion step in `catsib_item()` used a hardcoded
    threshold of 3 instead of the user-supplied `min.binsize` argument,
    causing inconsistent bin filtering behavior.
  - The reliability estimate `rho2` in `catsib_one()` was not clamped to
    $[0, 1]$, so when `errvar > sigma2` (e.g., very few items or
    purification cascade), a negative `rho2` reversed the regression
    correction direction, inflating the Type I error rate.
  - When `errvar >= sigma2` during purification, `rho2` collapsed to 0,
    causing all corrected scores to converge to the group mean. With
    group mean differences (impact), this produced empty bin data frames
    and an invalid purification result. A minimum floor of 0.05 is now
    enforced for `rho2` to preserve score spread.
- Fixed a critical bug in `covirt()` where the guessing parameter
  (`par[,3]`) was incorrectly passed as the difficulty parameter
  (`par[,2]`) to the `integrand()` function for DRM items. This caused
  the gradient computation to receive `c = b`, which zeroed out the
  $\partial P/\partial a$ and $\partial P/\partial b$ gradient
  components and produced a singular Fisher information matrix.
  Additionally, `NA` values in `par.3` for 1PLM and 2PLM items are now
  substituted with 0 prior to gradient evaluation to prevent `NA`
  propagation.

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
  precision and bias in Multistage-Adaptive Test (MST) panels using a
  recursion-based evaluation method introduced by Lim et al. (2020).

- Added a new function, `pcd2()`, which computes the Pseudo-count
  $D^{2}$ statistics (Cappaert et al., 2018; Stone, 2000) to detect item
  parameter drift.

# irtQ 0.2.0

- Introduced Warm's (1989) Weighted Likelihood (WL) estimation method to
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

  - `cac_rud`: This function implements Rudner's (2001, 2005) method for
    computing classification accuracy and consistency. It takes cut
    scores, ability estimates, standard errors, and optional weights as
    inputs and returns a list containing a confusion matrix, marginal
    and conditional classification accuracy and consistency indices, the
    probability of being assigned to each level category, and the cut
    scores used in the analysis.
  - `cac_lee`: This function implements Lee's (2010) method for
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
  non-standard way. Specifically, the problem arose when 0 wasn't used
  as the reference group and consecutive numbers (e.g., 1, 2, 3) weren't
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
