#' Item parameter estimation using MMLE-EM algorithm
#'
#' This function fits unidimensional item response theory (IRT) models to
#' mixed-format data comprising both dichotomous and polytomous items, using
#' marginal maximum likelihood estimation via the expectation–maximization
#' (MMLE-EM) algorithm (Bock & Aitkin, 1981). It also supports fixed item
#' parameter calibration (FIPC; Kim, 2006), a practical method for pretest (or
#' newly developed) item calibration in computerized adaptive testing (CAT).
#' FIPC enables the parameter estimates of pretest items to be placed on the
#' same scale as those of operational items (Ban et al., 2001). For dichotomous
#' items, the function supports the one-, two-, and three-parameter logistic
#' models. For polytomous items, it supports the graded response model (GRM) and
#' the (generalized) partial credit model (GPCM).
#'
#' @param x A data frame containing item metadata. This metadata is required to
#'   retrieve essential information for each item (e.g., number of score
#'   categories, IRT model type, etc.) necessary for calibration. You can create
#'   an empty item metadata frame using the function [irtQ::shape_df()].
#'
#'   When `use.startval = TRUE`, the item parameters specified in the metadata
#'   will be used as starting values for parameter estimation. If `x = NULL`,
#'   both `model` and `cats` arguments must be specified. Note that when `fipc =
#'   TRUE` to implement FIPC, item metadata for the test form must be supplied
#'   via the `x` argument. See **below** for more details. Default is `NULL`.
#' @param data A matrix of examinees' item responses corresponding to the items
#'   specified in the `x` argument. Rows represent examinees and columns
#'   represent items.
#' @param D A scaling constant used in IRT models to make the logistic function
#'   closely approximate the normal ogive function. A value of 1.7 is commonly
#'   used for this purpose. Default is 1.
#' @param model A character vector specifying the IRT model to fit each item.
#'   Available values are:
#'   - `"1PLM"`, `"2PLM"`, `"3PLM"`, `"DRM"` for dichotomous items
#'   - `"GRM"`, `"GPCM"` for polytomous items
#'
#'   Here, `"GRM"` denotes the graded response model and `"GPCM"` the
#'   (generalized) partial credit model. Note that `"DRM"` serves as a general
#'   label covering all three dichotomous IRT models. If a single model name is
#'   provided, it is recycled for all items. This argument is only used when `x
#'   = NULL` and `fipc = FALSE`. Default is `NULL`.
#' @param cats Numeric vector specifying the number of score categories per
#'   item. For dichotomous items, this should be 2. If a single value is
#'   supplied, it will be recycled across all items. When `cats = NULL` and all
#'   models specified in the `model` argument are dichotomous (`"1PLM"`,
#'   `"2PLM"`, `"3PLM"`, or `"DRM"`), the function defaults to 2 categories per
#'   item. This argument is used only when `x = NULL` and `fipc = FALSE`.
#'   Default is `NULL`.
#' @param item.id Character vector of item identifiers. If `NULL`, IDs are
#'   generated automatically. When `fipc = TRUE`, a provided `item.id` will
#'   override any IDs present in `x`. Default is `NULL`.
#' @param fix.a.1pl Logical. If `TRUE`, the slope parameters of all 1PLM items
#'   are fixed to `a.val.1pl`; otherwise, they are constrained to be equal and
#'   estimated. Default is `FALSE`.
#' @param fix.a.gpcm Logical. If `TRUE`, GPCM items are calibrated as PCM with
#'   slopes fixed to `a.val.gpcm`; otherwise, each item's slope is estimated.
#'   Default is `FALSE`.
#' @param fix.g Logical. If `TRUE`, all 3PLM guessing parameters are fixed to
#'   `g.val`; otherwise, each guessing parameter is estimated. Default is
#'   `FALSE`.
#' @param a.val.1pl Numeric. Value to which the slope parameters of 1PLM items
#'   are fixed when `fix.a.1pl = TRUE`. Default is 1.
#' @param a.val.gpcm Numeric. Value to which the slope parameters of GPCM items
#'   are fixed when `fix.a.gpcm = TRUE`. Default is 1.
#' @param g.val Numeric. Value to which the guessing parameters of 3PLM items
#'   are fixed when `fix.g = TRUE`. Default is 0.2.
#' @param use.aprior Logical. If `TRUE`, applies a prior distribution to all
#'   item discrimination (slope) parameters during calibration. Default is
#'   `FALSE`.
#' @param use.bprior Logical. If `TRUE`, applies a prior distribution to all
#'   item difficulty (or threshold) parameters during calibration. Default is
#'   `FALSE`.
#' @param use.gprior Logical. If `TRUE`, applies a prior distribution to all
#'   3PLM guessing parameters during calibration. Default is `TRUE`.
#' @param aprior,bprior,gprior A list specifying the prior distribution for all
#'   item discrimination (slope), difficulty (or threshold), guessing
#'   parameters. Three distributions are supported: Beta, Log-normal, and
#'   Normal. The list must have two elements:
#'   - `dist`: A character string, one of `"beta"`, `"lnorm"`, or `"norm"`.
#'   - `params`: A numeric vector of length two giving the distribution’s
#'   parameters. For details on each parameterization, see [stats::dbeta()],
#'   [stats::dlnorm()], and [stats::dnorm()].
#'
#'   Defaults are:
#'   - `aprior = list(dist = "lnorm", params = c(0.0, 0.5))`
#'   - `bprior = list(dist = "norm", params = c(0.0, 1.0))`
#'   - `gprior = list(dist = "beta", params = c(5, 16))`
#'
#'   for discrimination, difficulty, and guessing parameters, respectively.
#' @param missing  A value indicating missing responses in the data set. Default
#'   is `NA`.
#' @param Quadrature  A numeric vector of length two:
#'   - first element: number of quadrature points
#'   - second element: symmetric bound (absolute value) for those points
#'   For example, `c(49, 6)` specifies 49 evenly spaced points from –6 to 6.
#'   These points are used in the E-step of the EM algorithm. Default is `c(49,
#'   6)`.
#' @param weights A two-column matrix or data frame containing the quadrature
#'   points (in the first column) and their corresponding weights (in the second
#'   column) for the latent variable prior distribution. If not `NULL`, the
#'   scale of the latent ability distribution is fixed to match the scale of the
#'   provided quadrature points and weights. The weights and points can be
#'   conveniently generated using the function [irtQ::gen.weight()].
#'
#'   If `NULL`, a normal prior density is used instead, based on the
#'   information provided in the `Quadrature`, `group.mean`, and `group.var`
#'   arguments. Default is `NULL`.
#' @param group.mean A numeric value specifying the mean of the latent variable
#'   prior distribution when `weights = NULL`. Default is 0. This value is fixed
#'   to resolve the indeterminacy of the item parameter scale during
#'   calibration. However, the scale of the prior distribution is updated when
#'   FIPC is implemented.
#' @param group.var A positive numeric value specifying the variance of the
#'   latent variable prior distribution when `weights = NULL`. Default is 1.
#'   This value is fixed to resolve the indeterminacy of the item parameter
#'   scale during calibration. However, the scale of the prior distribution is
#'   updated when FIPC is implemented.
#' @param EmpHist Logical. If `TRUE`, the empirical histogram of the latent
#'   variable prior distribution is estimated simultaneously with the item
#'   parameters using the approach proposed by Woods (2007). Item calibration is
#'   conducted relative to the estimated empirical prior. See below for details.
#' @param use.startval Logical. If `TRUE`, the item parameters provided in the
#'   item metadata (i.e., the `x` argument) are used as starting values for item
#'   parameter estimation. Otherwise, internally generated starting values are
#'   used. Default is `FALSE`.
#' @param Etol A positive numeric value specifying the convergence criterion for
#'   the E-step of the EM algorithm. Default is `1e-4`.
#' @param MaxE A positive integer specifying the maximum number of iterations
#'   for the E-step in the EM algorithm. Default is `500`.
#' @param control A list of control parameters to be passed to the optimization
#'   function [stats::nlminb()]. These parameters control the M-step of the EM
#'   algorithm. For example, the maximum number of iterations in each M-step can
#'   be specified using `control = list(iter.max = 200)`. The default maximum
#'   number of iterations per M-step is 200. See [stats::nlminb()] for
#'   additional control options.
#' @param fipc Logical. If `TRUE`, fixed item parameter calibration (FIPC) is
#'   applied during item parameter estimation. When `fipc = TRUE`, the
#'   information on which items are fixed must be provided via either `fix.loc`
#'   or `fix.id`. See below for details.
#' @param fipc.method A character string specifying the FIPC method. Available
#'   options are:
#'   - `"OEM"`: No Prior Weights Updating and One EM Cycle (NWU-OEM; Wainer &
#'   Mislevy, 1990)
#'   - `"MEM"`: Multiple Prior Weights Updating and Multiple EM Cycles
#'   (MWU-MEM; Kim, 2006) When `fipc.method = "OEM"`, the maximum number of
#'   E-steps is automatically set to 1, regardless of the value specified in
#'   `MaxE`.
#' @param fix.loc A vector of positive integers specifying the row positions of
#'   the items to be fixed in the item metadata (i.e., `x`) when FIPC is
#'   implemented (i.e., `fipc = TRUE`). For example, suppose that five items
#'   located in the 1st, 2nd, 4th, 7th, and 9th rows of `x` should be fixed.
#'   Then use `fix.loc = c(1, 2, 4, 7, 9)`. Note that if `fix.id` is not `NULL`,
#'   the information provided in `fix.loc` is ignored. See below for details.
#' @param fix.id A character vector specifying the IDs of the items to be fixed
#'   when FIPC is implemented (i.e., `fipc = TRUE`). For example, suppose five
#'   items with IDs "CMC1", "CMC2", "CMC3", "CMC4", and "CMC5" are to be fixed,
#'   and that all item IDs are supplied via `item.id` column in the `x`
#'   argument. Then use `fix.id = c("CMC1", "CMC2", "CMC3", "CMC4", "CMC5")`.
#'   Note that if `fix.id` is not `NULL`, the information in `fix.loc` is
#'   ignored. See below for details.
#' @param se Logical. If `FALSE`, standard errors of the item parameter
#'   estimates are not computed. Default is `TRUE`.
#' @param verbose Logical. If `FALSE`, all progress messages, including
#'   information about the EM algorithm process, are suppressed. Default is
#'   `TRUE`.
#' @details A specific format of data frame should be used for the argument `x`.
#' The first column should contain item IDs, the second column should contain
#' the number of unique score categories for each item, and the third column
#' should specify the IRT model to be fitted to each item. Available IRT models
#' are:
#'   - `"1PLM"`, `"2PLM"`, `"3PLM"`, and `"DRM"` for dichotomous item data
#'   - `"GRM"` and `"GPCM"` for polytomous item data
#'
#' Note that `"DRM"` serves as a general label covering all dichotomous IRT
#' models (i.e., `"1PLM"`, `"2PLM"`, and `"3PLM"`), while `"GRM"` and `"GPCM"`
#' represent the graded response model and (generalized) partial credit model,
#' respectively.
#'
#' The subsequent columns should contain the item parameters for the specified
#' models. For dichotomous items, the fourth, fifth, and sixth columns represent
#' item discrimination (slope), item difficulty, and item guessing parameters,
#' respectively. When `"1PLM"` or `"2PLM"` is specified in the third column,
#' `NA`s must be entered in the sixth column for the guessing parameters.
#'
#' For polytomous items, the item discrimination (slope) parameter should appear
#' in the fourth column, and the item difficulty (or threshold) parameters for
#' category boundaries should occupy the fifth through the last columns. When
#' the number of unique score categories differs across items, unused parameter
#' cells should be filled with `NA`s.
#'
#' In the \pkg{irtQ} package, the threshold parameters for GPCM items are
#' expressed as the item location (or overall difficulty) minus the threshold
#' values for each score category. Note that when a GPCM item has *K* unique
#' score categories, *K - 1* threshold parameters are required, since the
#' threshold for the first category boundary is always fixed at 0. For example,
#' if a GPCM item has five score categories, four threshold parameters must be
#' provided.
#'
#' An example of a data frame for a single-format test is shown below:
#' \tabular{lrlrrrrr}{
#'   ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \cr
#'   ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \cr
#'   ITEM3  \tab 2 \tab 3PLM \tab 1.736 \tab  1.501 \tab  0.203 \cr
#'   ITEM4  \tab 2 \tab 3PLM \tab 0.835 \tab -1.049 \tab  0.182 \cr
#'   ITEM5  \tab 2 \tab DRM \tab 0.926 \tab  0.394 \tab  0.099
#' }
#'
#' An example of a data frame for a mixed-format test is shown below:
#' \tabular{lrlrrrrr}{
#'   ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \tab         NA \tab         NA\cr
#'   ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \tab         NA \tab         NA\cr
#'   ITEM3  \tab 2 \tab 3PLM \tab 0.926 \tab  0.394 \tab  0.099 \tab         NA \tab         NA\cr
#'   ITEM4  \tab 2 \tab DRM \tab 1.052 \tab -0.407 \tab  0.201 \tab         NA \tab         NA\cr
#'   ITEM5  \tab 4 \tab GRM  \tab 1.913 \tab -1.869 \tab -1.238 \tab -0.714 \tab         NA \cr
#'   ITEM6  \tab 5 \tab GRM  \tab 1.278 \tab -0.724 \tab -0.068 \tab  0.568 \tab  1.072\cr
#'   ITEM7  \tab 4 \tab GPCM  \tab 1.137 \tab -0.374 \tab  0.215 \tab  0.848 \tab         NA \cr
#'   ITEM8  \tab 5 \tab GPCM  \tab 1.233 \tab -2.078 \tab -1.347 \tab -0.705 \tab -0.116
#' }
#'
#' See the *IRT Models* section in the [irtQ-package] documentation for more
#' details about the IRT models used in the \pkg{irtQ} package. A convenient way
#' to create a data frame for the argument `x` is by using the function
#' [irtQ::shape_df()].
#'
#' To fit IRT models to data, the item response data must be accompanied by
#' information on the IRT model and the number of score categories for each
#' item. There are two ways to provide this information:
#'    1. Supply item metadata to the argument `x`. As explained above, such
#'    metadata can be easily created using [irtQ::shape_df()].
#'    2. Specify the IRT models and score category information directly
#'    through the arguments `model` and `cats`.
#'
#' If `x = NULL`, the function uses the information specified in `model` and
#' `cats`.
#'
#' To implement FIPC, the item metadata must be provided via the `x` argument.
#' This is because the item parameters of the fixed items in the metadata are
#' used to estimate the characteristics of the underlying latent variable prior
#' distribution when calibrating the remaining (freely estimated) items. More
#' specifically, the latent prior distribution is estimated based on the
#' fixed items, and then used to calibrate the new (pretest) items so that their
#' parameters are placed on the same scale as those of the fixed items (Kim,
#' 2006).The full item metadata, including both fixed and non-fixed items, can be
#' conveniently created using the [irtQ::shape_df_fipc()] function.
#'
#' In terms of approaches for FIPC, Kim (2006) described five different methods.
#' Among them, two methods are available in the [irtQ::est_irt()] function. The
#' first method is `"NWU-OEM"`, which uses a single E-step in the EM
#' algorithm (involving only the fixed items) followed by a single M-step
#' (involving only the non-fixed items). This method was proposed by Wainer and
#' Mislevy (1990) in the context of online calibration and can be implemented by
#' setting `fipc.method = "OEM"`.
#'
#' The second method is `"MWU-MEM"`, which iteratively updates the latent
#' variable prior distribution and estimates the parameters of the non-fixed
#' items. In this method, the same procedure as the NWU-OEM approach is applied
#' during the first EM cycle. From the second cycle onward, both the parameters
#' of the non-fixed items and the weights of the prior distribution are
#' concurrently updated. This method can be implemented by setting `fipc.method
#' = "MEM"`. See Kim (2006) for more details.
#'
#' When `fipc = TRUE`, information about which items are to be fixed must be
#' provided via either the `fix.loc` or `fix.id` argument. For example, suppose
#' that five items with IDs "CMC1", "CMC2", "CMC3", "CMC4", and "CMC5" should
#' be fixed, and all item IDs are provided via the `x` or `item.id` argument.
#' Also, assume these five items are located in the 1st through 5th rows of the
#' item metadata (i.e., `x`). In this case, the fixed items can be specified
#' using either `fix.loc = c(1, 2, 3, 4, 5)` or
#' `fix.id = c("CMC1", "CMC2", "CMC3", "CMC4", "CMC5")`. Note that if both
#' `fix.loc` and `fix.id` are not `NULL`, the information in `fix.loc` is
#' ignored.
#'
#' When `EmpHist = TRUE`, the empirical histogram of the latent variable prior
#' distribution (i.e., the densities at the quadrature points) is estimated
#' simultaneously with the item parameters. If `EmpHist = TRUE` and
#' `fipc = TRUE`, the scale parameters of the empirical
#' prior distribution (e.g., mean and variance) are also estimated.
#' If `EmpHist = TRUE` and `fipc = FALSE`, the scale parameters are fixed to
#' the values specified in `group.mean` and `group.var`. When `EmpHist = FALSE`,
#' a normal  prior distribution is used instead. If `fipc = TRUE`, the scale
#' parameters of this normal prior are estimated along with the item parameters.
#' If `fipc = FALSE`, they are fixed to the values specified in `group.mean` and
#' `group.var`.
#'
#' @return This function returns an object of class `est_irt`. The returned
#'   object contains the following components:
#'
#'   \item{estimates}{A data frame containing both the item parameter estimates
#'   and their corresponding standard errors.}
#'
#'   \item{par.est}{A data frame of item parameter estimates, structured
#'   according to the item metadata format.}
#'
#'   \item{se.est}{A data frame of standard errors for the item parameter
#'   estimates, computed using the cross-product approximation method
#'   (Meilijson, 1989).}
#'
#'   \item{pos.par}{A data frame indicating the position index of each estimated
#'   item parameter. The position information is useful for interpreting
#'   the variance-covariance matrix of item parameter estimates}
#'
#'   \item{covariance}{A variance-covariance matrix of the item parameter
#'   estimates.}
#'
#'   \item{loglikelihood}{The marginal log-likelihood, calculated as the sum of
#'   the log-likelihoods across all items.}
#'
#'   \item{aic}{Akaike Information Criterion (AIC) based on the log-likelihood.}
#'
#'   \item{bic}{Bayesian Information Criterion (BIC) based on the
#'   log-likelihood.}
#'
#'   \item{group.par}{A data frame containing the mean, variance, and standard
#'   deviation of the latent variable prior distribution.}
#'
#'   \item{weights}{A two-column data frame of quadrature points (column 1) and
#'   corresponding weights (column 2) of the (updated) latent prior distribution.}
#'
#'   \item{posterior.dist}{A matrix of normalized posterior densities for all
#'   response patterns at each quadrature point. Rows and columns represent
#'   response patterns and quadrature points, respectively.}
#'
#'   \item{data}{A data frame of examinees' response data.}
#'
#'   \item{scale.D}{The scaling factor used in the IRT model.}
#'
#'   \item{ncase}{The total number of response patterns.}
#'
#'   \item{nitem}{The total number of items in the response data.}
#'
#'   \item{Etol}{The convergence criterion for the E-step of the EM algorithm.}
#'
#'   \item{MaxE}{The maximum number of E-steps allowed in the EM algorithm.}
#'
#'   \item{aprior}{A list describing the prior distribution used for discrimination
#'    parameters.}
#'
#'   \item{bprior}{A list describing the prior distribution used for difficulty
#'    parameters.}
#'
#'   \item{gprior}{A list describing the prior distribution used for guessing
#'    parameters.}
#'
#'   \item{npar.est}{The total number of parameters estimated.}
#'
#'   \item{niter}{The number of completed EM cycles.}
#'
#'   \item{maxpar.diff}{The maximum absolute change in parameter estimates at
#'   convergence.}
#'
#'   \item{EMtime}{Time (in seconds) spent on EM cycles.}
#'
#'   \item{SEtime}{Time (in seconds) spent computing standard errors.}
#'
#'   \item{TotalTime}{Total computation time (in seconds).}
#'
#'   \item{test.1}{First-order test result indicating whether the gradient
#'   sufficiently vanished for solution stability.}
#'
#'   \item{test.2}{Second-order test result indicating whether the information matrix
#'   is positive definite, a necessary condition for identifying a local maximum.}
#'
#'   \item{var.note}{A note indicating whether the variance-covariance matrix
#'   was successfully obtained from the information matrix.}
#'
#'   \item{fipc}{Logical. Indicates whether FIPC was used.}
#'
#'   \item{fipc.method}{The method used for FIPC.}
#'
#'   \item{fix.loc}{A vector of integers specifying the row locations of fixed items
#'   when FIPC was applied.}
#'
#'   Note that you can easily extract components from the output using the
#'   [irtQ::getirt()] function.
#'
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::shape_df()], [irtQ::shape_df_fipc()], [irtQ::getirt()]
#'
#' @references Ban, J. C., Hanson, B. A., Wang, T., Yi, Q., & Harris, D., J.
#'   (2001) A comparative study of on-line pretest item calibration/scaling
#'   methods in computerized adaptive testing. *Journal of Educational
#'   Measurement, 38*(3), 191-212.
#'
#'   Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of
#'   item parameters: Application of an EM algorithm. *Psychometrika, 46*,
#'   443-459.
#'
#'   Kim, S. (2006). A comparative study of IRT fixed parameter calibration
#'   methods. *Journal of Educational Measurement, 43*(4), 355-381.
#'
#'   Meilijson, I. (1989). A fast improvement to the EM algorithm on its own
#'   terms. *Journal of the Royal Statistical Society: Series B
#'   (Methodological), 51*, 127-138.
#'
#'   Stocking, M. L. (1988). *Scale drift in on-line calibration* (Research Rep.
#'   88-28). Princeton, NJ: ETS.
#'
#'   Wainer, H., & Mislevy, R. J. (1990). Item response theory, item
#'   calibration, and proficiency estimation. In H. Wainer (Ed.), *Computer
#'   adaptive testing: A primer* (Chap. 4, pp.65-102). Hillsdale, NJ: Lawrence
#'   Erlbaum.
#'
#'   Woods, C. M. (2007). Empirical histograms in item response theory with
#'   ordinal data. *Educational and Psychological Measurement, 67*(1), 73-87.
#'
#' @examples
#' \donttest{
#'
#' ## --------------------------------------------------------------
#' ## 1. Item parameter estimation for dichotomous item data (LSAT6)
#' ## --------------------------------------------------------------
#' # Fit the 1PL model to LSAT6 data and estimate a common slope parameter
#' # (i.e., constrain slope parameters to be equal)
#' (mod.1pl.c <- est_irt(data = LSAT6, D = 1, model = "1PLM", cats = 2,
#'                       fix.a.1pl = FALSE))
#'
#' # Display a summary of the estimation results
#' summary(mod.1pl.c)
#'
#' # Extract the item parameter estimates
#' getirt(mod.1pl.c, what = "par.est")
#'
#' # Extract the standard error estimates
#' getirt(mod.1pl.c, what = "se.est")
#'
#' # Fit the 1PL model to LSAT6 data and fix slope parameters to 1.0
#' (mod.1pl.f <- est_irt(data = LSAT6, D = 1, model = "1PLM", cats = 2,
#'                       fix.a.1pl = TRUE, a.val.1pl = 1))
#'
#' # Display a summary of the estimation results
#' summary(mod.1pl.f)
#'
#' # Fit the 2PL model to LSAT6 data
#' (mod.2pl <- est_irt(data = LSAT6, D = 1, model = "2PLM", cats = 2))
#'
#' # Display a summary of the estimation results
#' summary(mod.2pl)
#'
#' # Assess the model fit for the 2PL model using the S-X2 fit statistic
#' (sx2fit.2pl <- sx2_fit(x = mod.2pl))
#'
#' # Compute item and test information functions at a range of theta values
#' theta <- seq(-4, 4, 0.1)
#' (info.2pl <- info(x = mod.2pl, theta = theta))
#'
#' # Plot the test characteristic curve (TCC)
#' (trace.2pl <- traceline(x = mod.2pl, theta = theta))
#' plot(trace.2pl)
#'
#' # Plot the item characteristic curve (ICC) for the first item
#' plot(trace.2pl, item.loc = 1)
#'
#' # Fit the 2PL model and simultaneously estimate an empirical histogram
#' # of the latent variable prior distribution
#' # Also apply a looser convergence threshold for the E-step
#' (mod.2pl.hist <- est_irt(data = LSAT6, D = 1, model = "2PLM", cats = 2,
#'                          EmpHist = TRUE, Etol = 0.001))
#' (emphist <- getirt(mod.2pl.hist, what = "weights"))
#' plot(emphist$weight ~ emphist$theta, type = "h")
#'
#' # Fit the 3PL model and apply a Beta prior to the guessing parameters
#' (mod.3pl <- est_irt(
#'   data = LSAT6, D = 1, model = "3PLM", cats = 2, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 16))
#' ))
#'
#' # Display a summary of the estimation results
#' summary(mod.3pl)
#'
#' # Fit the 3PL model and fix the guessing parameters at 0.2
#' (mod.3pl.f <- est_irt(data = LSAT6, D = 1, model = "3PLM", cats = 2,
#'                       fix.g = TRUE, g.val = 0.2))
#'
#' # Display a summary of the estimation results
#' summary(mod.3pl.f)
#'
#' # Fit different dichotomous models to each item in the LSAT6 data:
#' # Fit the constrained 1PL model to items 1–3, the 2PL model to item 4,
#' # and the 3PL model with a Beta prior on guessing to item 5
#' (mod.drm.mix <- est_irt(
#'   data = LSAT6, D = 1, model = c("1PLM", "1PLM", "1PLM", "2PLM", "3PLM"),
#'   cats = 2, fix.a.1pl = FALSE, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 16))
#' ))
#'
#' # Display a summary of the estimation results
#' summary(mod.drm.mix)
#'
#' ## -------------------------------------------------------------------
#' ## 2. Item parameter estimation for mixed-format data (simulated data)
#' ## -------------------------------------------------------------------
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Extract item metadata
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df
#'
#' # Modify the item metadata so that the 39th and 40th items use the GPCM
#' x[39:40, 3] <- "GPCM"
#'
#' # Generate 1,000 examinees' latent abilities from N(0, 1)
#' set.seed(37)
#' score1 <- rnorm(1000, mean = 0, sd = 1)
#'
#' # Simulate item response data
#' sim.dat1 <- simdat(x = x, theta = score1, D = 1)
#'
#' # Fit the 3PL model to all dichotomous items, the GPCM to items 39 and 40,
#' # and the GRM to items 53, 54, and 55.
#' # Use a Beta prior for guessing parameters, a log-normal prior for slope
#' # parameters, and a normal prior for difficulty (threshold) parameters.
#' # Also, specify the argument `x` to provide IRT model and score category information.
#' item.meta <- shape_df(item.id = x$id, cats = x$cats, model = x$model,
#'   default.par = TRUE)
#' (mod.mix1 <- est_irt(
#'   x = item.meta, data = sim.dat1, D = 1, use.aprior = TRUE, use.bprior = TRUE,
#'   use.gprior = TRUE,
#'   aprior = list(dist = "lnorm", params = c(0.0, 0.5)),
#'   bprior = list(dist = "norm", params = c(0.0, 2.0)),
#'   gprior = list(dist = "beta", params = c(5, 16))
#' ))
#'
#' # Display a summary of the estimation results
#' summary(mod.mix1)
#'
#' # Estimate examinees' latent scores using MLE and the estimated item parameters
#' (score.mle <- est_score(x = mod.mix1, method = "ML", range = c(-4, 4), ncore = 2))
#'
#' # Compute traditional model-fit statistics
#' (fit.mix1 <- irtfit(
#'   x = mod.mix1, score = score.mle$est.theta, group.method = "equal.width",
#'   n.width = 10, loc.theta = "middle"
#' ))
#'
#' # Residual plot for the first item (dichotomous)
#' plot(
#'   x = fit.mix1, item.loc = 1, type = "both", ci.method = "wald",
#'   show.table = TRUE, ylim.sr.adjust = TRUE
#' )
#'
#' # Residual plot for the last item (polytomous)
#' plot(
#'   x = fit.mix1, item.loc = 55, type = "both", ci.method = "wald",
#'   show.table = FALSE, ylim.sr.adjust = TRUE
#' )
#'
#' # Fit the 2PL model to all dichotomous items, the GPCM to items 39 and 40,
#' # and the GRM to items 53, 54, and 55.
#' # Provide IRT model and score category information via `model` and `cats`
#' # arguments.
#' (mod.mix2 <- est_irt(
#'   data = sim.dat1, D = 1,
#'   model = c(rep("2PLM", 38), rep("GPCM", 2), rep("2PLM", 12), rep("GRM", 3)),
#'   cats = c(rep(2, 38), rep(5, 2), rep(2, 12), rep(5, 3))
#' ))
#'
#' # Display a summary of the estimation results
#' summary(mod.mix2)
#'
#' # Fit the 2PL model to all dichotomous items, the GPCM to items 39 and 40,
#' # and the GRM to items 53, 54, and 55.
#' # Also estimate the empirical histogram of the latent prior distribution.
#' # Provide IRT model and score category information via `model` and `cats` arguments.
#' (mod.mix3 <- est_irt(
#'   data = sim.dat1, D = 1,
#'   model = c(rep("2PLM", 38), rep("GPCM", 2), rep("2PLM", 12), rep("GRM", 3)),
#'   cats = c(rep(2, 38), rep(5, 2), rep(2, 12), rep(5, 3)), EmpHist = TRUE
#' ))
#' (emphist <- getirt(mod.mix3, what = "weights"))
#' plot(emphist$weight ~ emphist$theta, type = "h")
#'
#' # Fit the 2PL model to all dichotomous items, the PCM to items 39 and 40 by
#' # fixing slope parameters to 1, and the GRM to items 53, 54, and 55.
#' # Provide IRT model and score category information via `model` and `cats` arguments.
#' (mod.mix4 <- est_irt(
#'   data = sim.dat1, D = 1,
#'   model = c(rep("2PLM", 38), rep("GPCM", 2), rep("2PLM", 12), rep("GRM", 3)),
#'   cats = c(rep(2, 38), rep(5, 2), rep(2, 12), rep(5, 3)),
#'   fix.a.gpcm = TRUE, a.val.gpcm = 1
#' ))
#'
#' # Display a summary of the estimation results
#' summary(mod.mix4)
#'
#' ## ----------------------------------------------------------------
#' ## 3. Fixed item parameter calibration (FIPC) for mixed-format data
#' ##    (simulated)
#' ## ----------------------------------------------------------------
#' ## Import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
#'
#' # Select item metadata
#' x <- bring.flexmirt(file = flex_sam, "par")$Group1$full_df
#'
#' # Generate 1,000 examinees' latent abilities from N(0.4, 1.3)
#' set.seed(20)
#' score2 <- rnorm(1000, mean = 0.4, sd = 1.3)
#'
#' # Simulate response data
#' sim.dat2 <- simdat(x = x, theta = score2, D = 1)
#'
#' # Fit the 3PL model to all dichotomous items and the GRM to all polytomous items
#' # Fix five 3PL items (1st–5th) and three GRM items (53rd–55th)
#' # Also estimate the empirical histogram of the latent variable distribution
#' # Use the MEM method
#' fix.loc <- c(1:5, 53:55)
#' (mod.fix1 <- est_irt(
#'   x = x, data = sim.dat2, D = 1, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 16)), EmpHist = TRUE,
#'   Etol = 1e-3, fipc = TRUE, fipc.method = "MEM", fix.loc = fix.loc
#' ))
#'
#' # Extract group-level parameter estimates
#' (prior.par <- mod.fix1$group.par)
#'
#' # Visualize the empirical prior distribution
#' (emphist <- getirt(mod.fix1, what = "weights"))
#' plot(emphist$weight ~ emphist$theta, type = "h")
#'
#' # Display a summary of the estimation results
#' summary(mod.fix1)
#'
#' # Alternatively, fix the same items by providing their item IDs
#' # using the `fix.id` argument. In this case, set `fix.loc = NULL`
#' fix.id <- c(x$id[1:5], x$id[53:55])
#' (mod.fix1 <- est_irt(
#'   x = x, data = sim.dat2, D = 1, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 16)), EmpHist = TRUE,
#'   Etol = 1e-3, fipc = TRUE, fipc.method = "MEM", fix.loc = NULL,
#'   fix.id = fix.id
#' ))
#'
#' # Display a summary of the estimation results
#' summary(mod.fix1)
#'
#' # Fit the 3PL model to all dichotomous items and the GRM to all polytomous items
#' # Fix the same items as before (1st–5th and 53rd–55th)
#' # This time, do not estimate the empirical histogram of the latent prior
#' # Instead, estimate the scale of the normal prior distribution
#' # Use the MEM method
#' fix.loc <- c(1:5, 53:55)
#' (mod.fix2 <- est_irt(
#'   x = x, data = sim.dat2, D = 1, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 16)), EmpHist = FALSE,
#'   Etol = 1e-3, fipc = TRUE, fipc.method = "MEM", fix.loc = fix.loc
#' ))
#'
#' # Extract group-level parameter estimates
#' (prior.par <- mod.fix2$group.par)
#'
#' # Visualize the prior distribution
#' (emphist <- getirt(mod.fix2, what = "weights"))
#' plot(emphist$weight ~ emphist$theta, type = "h")
#'
#' # Fit the 3PL model to all dichotomous items and the GRM to all polytomous items
#' # Fix only the five 3PL items (1st–5th) and estimate the empirical histogram
#' # Use the OEM method (i.e., only one EM cycle is used)
#' fix.loc <- c(1:5)
#' (mod.fix3 <- est_irt(
#'   x = x, data = sim.dat2, D = 1, use.gprior = TRUE,
#'   gprior = list(dist = "beta", params = c(5, 16)), EmpHist = TRUE,
#'   Etol = 1e-3, fipc = TRUE, fipc.method = "OEM", fix.loc = fix.loc
#' ))
#'
#' # Extract group-level parameter estimates
#' (prior.par <- mod.fix3$group.par)
#'
#' # Visualize the prior distribution
#' (emphist <- getirt(mod.fix3, what = "weights"))
#' plot(emphist$weight ~ emphist$theta, type = "h")
#'
#' # Display a summary of the estimation results
#' summary(mod.fix3)
#'
#' # Fit the 3PL model to all dichotomous items and the GRM to all polytomous items
#' # Fix all 55 items and estimate only the latent ability distribution
#' # Use the MEM method
#' fix.loc <- c(1:55)
#' (mod.fix4 <- est_irt(
#'   x = x, data = sim.dat2, D = 1, EmpHist = TRUE,
#'   Etol = 1e-3, fipc = TRUE, fipc.method = "MEM", fix.loc = fix.loc
#' ))
#'
#' # Extract group-level parameter estimates
#' (prior.par <- mod.fix4$group.par)
#'
#' # Visualize the prior distribution
#' (emphist <- getirt(mod.fix4, what = "weights"))
#' plot(emphist$weight ~ emphist$theta, type = "h")
#'
#' # Display a summary of the estimation results
#' summary(mod.fix4)
#'
#' # Alternatively, fix all 55 items by providing their item IDs
#' # using the `fix.id` argument. In this case, set `fix.loc = NULL`
#' fix.id <- x$id
#' (mod.fix4 <- est_irt(
#'   x = x, data = sim.dat2, D = 1, EmpHist = TRUE,
#'   Etol = 1e-3, fipc = TRUE, fipc.method = "MEM", fix.loc = NULL,
#'   fix.id = fix.id
#' ))
#'
#' # Display a summary of the estimation results
#' summary(mod.fix4)
#'
#' }
#'
#' @export
#'
est_irt <- function(x = NULL,
                    data,
                    D = 1,
                    model = NULL,
                    cats = NULL,
                    item.id = NULL,
                    fix.a.1pl = FALSE,
                    fix.a.gpcm = FALSE,
                    fix.g = FALSE,
                    a.val.1pl = 1,
                    a.val.gpcm = 1,
                    g.val = .2,
                    use.aprior = FALSE,
                    use.bprior = FALSE,
                    use.gprior = TRUE,
                    aprior = list(dist = "lnorm", params = c(0.0, 0.5)),
                    bprior = list(dist = "norm", params = c(0.0, 1.0)),
                    gprior = list(dist = "beta", params = c(5, 16)),
                    missing = NA,
                    Quadrature = c(49, 6.0),
                    weights = NULL,
                    group.mean = 0.0,
                    group.var = 1.0,
                    EmpHist = FALSE,
                    use.startval = FALSE,
                    Etol = 1e-04,
                    MaxE = 500,
                    control = list(iter.max = 200),
                    fipc = FALSE,
                    fipc.method = "MEM",
                    fix.loc = NULL,
                    fix.id = NULL,
                    se = TRUE,
                    verbose = TRUE) {

  # match.call
  cl <- match.call()

  # item parameter estimation
  if (!fipc) {
    # item parameter estimation using MMLE-EM algorithm
    est_par <- est_irt_em(
      x = x, data = data, D = D, model = model, cats = cats, item.id = item.id,
      fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g,
      a.val.1pl = a.val.1pl, a.val.gpcm = a.val.gpcm, g.val = g.val,
      use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
      aprior = aprior, bprior = bprior, gprior = gprior, missing = missing,
      Quadrature = Quadrature, weights = weights, group.mean = group.mean,
      group.var = group.var, EmpHist = EmpHist, use.startval = use.startval,
      Etol = Etol, MaxE = MaxE, control = control, se = se, verbose = verbose
    )
  } else {
    # implement FIPC method
    est_par <- est_irt_fipc(
      x = x, data = data, D = D, item.id = item.id, fix.a.1pl = fix.a.1pl,
      fix.a.gpcm = fix.a.gpcm, fix.g = fix.g, a.val.1pl = a.val.1pl,
      a.val.gpcm = a.val.gpcm, g.val = g.val, use.aprior = use.aprior,
      use.bprior = use.bprior, use.gprior = use.gprior, aprior = aprior,
      bprior = bprior, gprior = gprior, missing = missing, Quadrature = Quadrature,
      weights = weights, group.mean = group.mean, group.var = group.var,
      EmpHist = EmpHist, use.startval = use.startval, Etol = Etol, MaxE = MaxE,
      control = control, fipc = TRUE, fipc.method = fipc.method, fix.loc = fix.loc,
      fix.id = fix.id, se = se, verbose = verbose
    )
  }

  # return the estimation results
  class(est_par) <- "est_irt"
  est_par$call <- cl
  est_par
}



# This function estimates item parameters using MMLE-EM algorithm
#' @importFrom Rfast colsums
#' @import dplyr
est_irt_em <- function(x = NULL,
                       data,
                       D = 1,
                       model = NULL,
                       cats = NULL,
                       item.id = NULL,
                       fix.a.1pl = FALSE,
                       fix.a.gpcm = FALSE,
                       fix.g = FALSE,
                       a.val.1pl = 1,
                       a.val.gpcm = 1,
                       g.val = .2,
                       use.aprior = FALSE,
                       use.bprior = FALSE,
                       use.gprior = TRUE,
                       aprior = list(dist = "lnorm", params = c(0.0, 0.5)),
                       bprior = list(dist = "norm", params = c(0.0, 1.0)),
                       gprior = list(dist = "beta", params = c(5, 16)),
                       missing = NA,
                       Quadrature = c(49, 6.0),
                       weights = NULL,
                       group.mean = 0,
                       group.var = 1,
                       EmpHist = FALSE,
                       use.startval = FALSE,
                       Etol = 1e-04,
                       MaxE = 500,
                       control = list(eval.max = 200, iter.max = 200),
                       se = TRUE,
                       verbose = TRUE) {

  # check start time
  start.time <- Sys.time()

  ## ---------------------------------------------------------------------
  # prepare the item parameter estimation
  ## ---------------------------------------------------------------------
  # check if the starting values are available
  if (use.startval & is.null(x)) {
    stop(paste0(
      "To use starting values for item parameter estimation, \n",
      "the item metadata must be specified in the argument 'x'."
    ), call. = FALSE)
  }

  # transform a data set to matrix
  data <- data.matrix(data)

  # extract information about the number of score categories and models
  if (verbose) {
    cat("Parsing input...", "\n")
  }
  if (!is.null(x)) {
    # confirm and correct all item metadata information
    x <- confirm_df(x)

    # extract item ids for all items
    id <- x$id

    # extract score categories for all items
    cats <- x$cats

    # extract model names for all items
    model <- x$model

    # if use.startval = FALSE
    # create a new item medadata object with starting values
    if (!use.startval) {
      x <- startval_df(cats = cats, model = model, item.id = id)
    }
  } else {
    # make the model names as upper cases
    model <- toupper(model)

    # when a character string scalar is provided in the model argument
    # create a vector of the same model names
    if (length(model) == 1) {
      model <- rep(model, ncol(data))
    }

    # check if the score category information is provided in the cats argument
    if (is.null(cats)) {
      if (all(model %in% c("1PLM", "2PLM", "3PLM", "DRM"))) {
        cats <- rep(2, ncol(data))
      } else {
        stop("The number of score categories for the items should be specified in the argument 'cats'.", call. = FALSE)
      }
    }

    # when an integer scalar is provided in the cats argument
    # create a vector of the same category numbers
    if (length(cats) == 1) {
      cats <- rep(cats, ncol(data))
    }

    # create an empty item metadata data frame
    x <- startval_df(cats = cats, model = model, item.id = item.id)

    # confirm and correct all item metadata information
    x <- confirm_df(x)

    # extract item ids for all items
    id <- x$id

    # extract score categories for all items
    cats <- x$cats

    # extract model names for all items
    model <- x$model
  }

  # check the total number of item in the response data set
  nitem <- ncol(data)

  # check the total number of examinees
  nstd <- nrow(data)

  # check whether included data are correct
  if (nrow(x) != nitem) stop("The number of items included in 'x' and 'data' must be the same.", call. = FALSE)

  # recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # check the number of item responses across all items
  n.resp <- Rfast::colsums(!is.na(data))

  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)
  if (length(loc_allmiss) > 0L) {
    memo2 <- paste0(paste0("item ", loc_allmiss, collapse = ", "), " has/have no item response data. \n")
    stop(memo2, call. = FALSE)
  }

  # find the location of 1PLM items in which slope parameters should be constrained to be equal
  # also, find the location of items with other models
  if ("1PLM" %in% model & !fix.a.1pl) {
    loc_1p_const <- which(model == "1PLM")
    loc_else <- which(model != "1PLM")

    # count the number of 1PLM items to be constrained
    n.1PLM <- length(loc_1p_const)
  } else {
    loc_1p_const <- NULL
    loc_else <- 1:nrow(x)
    n.1PLM <- NULL
  }

  # record the original location of item parameters to be estimated, and
  # the relocated position of item parameters when computing
  # the variance-covariance matrix of item parameter estimates
  param_loc <- parloc(
    x = x, loc_1p_const = loc_1p_const, loc_else = loc_else,
    fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g
  )

  ## ---------------------------------------------------------------------
  # conduct item parameter estimation using MMLE-EM algorithm
  ## ---------------------------------------------------------------------
  # create initial weights of prior ability distribution when it is not specified
  if (is.null(weights)) {
    # create a vector of quad-points
    quadpt <- seq(-Quadrature[2], Quadrature[2], length.out = Quadrature[1])

    # create a two column data frame to contain the quad-points and weights
    weights <- gen.weight(dist = "norm", mu = group.mean, sigma = sqrt(group.var), theta = quadpt)
    n.quad <- length(quadpt)
  } else {
    quadpt <- weights[, 1]
    n.quad <- length(quadpt)
    moments.tmp <- cal_moment(node = quadpt, weight = weights[, 2])
    group.mean <- moments.tmp[1]
    group.var <- moments.tmp[2]
  }

  # factorize the response values
  resp <- purrr::map2(
    .x = data.frame(data, stringsAsFactors = FALSE), .y = cats,
    .f = function(k, m) factor(k, levels = (seq_len(m) - 1))
  )

  # create a contingency table of score categories for each item
  # and then, transform the table to a matrix format
  std.id <- 1:nstd
  freq.cat <-
    purrr::map(
      .x = resp,
      .f = function(k) {
        stats::xtabs(~ std.id + k,
          na.action = stats::na.pass, addNA = FALSE
        ) %>%
          # as.numeric() %>%
          matrix(nrow = length(k))
      }
    )

  # delete 'resp' object
  rm(resp, envir = environment(), inherits = FALSE)

  # break down the item metadata into several elements
  elm_item <- breakdown(x)

  # classify the items into DRM and PRM item groups
  idx.item <- idxfinder(elm_item)
  idx.drm <- idx.item$idx.drm
  idx.prm <- idx.item$idx.prm

  # then, divide items into three groups
  # (1) DRM 1PL items with the constrained slop
  # (2) other DRM items
  # (3) PRM items
  if (sum(loc_else) == 0) {
    drm.else <- NULL
  } else {
    drm.else <- loc_else
  }
  idx4est <- list(
    drm.slc = loc_1p_const,
    drm.else = drm.else,
    prm = idx.prm
  )

  # divide the data set for the mixed-item format
  datlist <- divide_data(data = data, idx.item = idx.item, freq.cat = freq.cat)
  data_drm <- cbind(datlist$data_drm_q, datlist$data_drm_p)
  data_prm <- datlist$data_prm
  data_all <- datlist$data_all

  # delete 'datlist' object
  rm(datlist, envir = environment(), inherits = FALSE)

  # create the lower and upper bounds of the item parameters
  parbd <- lubound(model, cats, n.1PLM, idx4est, fix.a.1pl, fix.g, fix.a.gpcm)

  # find the columns of the frequency matrix corresponding to all items
  cols.item <- cols4item(nitem, cats, loc_1p_const)

  # estimation
  if (verbose) {
    cat("Estimating item parameters...", "\n")
  }

  # implement EM algorithm
  time1 <- Sys.time()
  # par.history <- list()
  for (r in 1:MaxE) {
    # implement E-step
    estep <- Estep(
      elm_item = elm_item, idx.drm = idx.drm, idx.prm = idx.prm,
      data_drm = data_drm, data_prm = data_prm, data_all = data_all,
      weights = weights, D = D
    )

    # implement M-step
    mstep <- Mstep(
      estep = estep, id = id, cats = cats, model = model, quadpt = quadpt,
      n.quad = n.quad, D = D, cols.item = cols.item, loc_1p_const = loc_1p_const,
      loc_else = loc_else, idx4est = idx4est, n.1PLM = n.1PLM, EmpHist = EmpHist,
      weights = weights, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g,
      a.val.1pl = a.val.1pl, a.val.gpcm = a.val.gpcm, g.val = g.val,
      use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
      aprior = aprior, bprior = bprior, gprior = gprior, group.mean = group.mean,
      group.var = group.var, nstd = nstd, Quadrature = Quadrature, control = control,
      iter = r, fipc = FALSE, reloc.par = param_loc$reloc.par, parbd = parbd
    )

    # compute the difference between previous and updated item parameter estimates
    diff_par <- mstep$elm_item$pars - elm_item$pars
    max.diff <- abs(max(diff_par, na.rm = TRUE))

    # loglikelihood value
    llike <- mstep$loglike

    # print
    if (verbose) {
      cat("\r", paste0(
        "EM iteration: ", r, ", Loglike: ", format(round(mstep$loglike, 4), nsmall = 4),
        ", Max-Change: ", format(round(max.diff, 6), nsmall = 5)
      ))
    }

    # check the convergence of EM algorithm
    converge <- max.diff <= Etol

    # extract the updated item parameter estimates
    elm_item$pars <- mstep$elm_item$pars
    # par.history[[r]] <- mstep$par_vec

    # extract the updated quadrature points and
    # the corresponding weights of the prior population density
    weights <- mstep$weights

    # delete 'estep' and 'mstep' object
    # rm(estep, mstep, envir=environment(), inherits = FALSE)

    # terminate the EM step if the convergence criterion is satisfied
    if (converge | r == MaxE) {
      break
    }
  }

  if (verbose) {
    cat("", "\n")
  }
  time2 <- Sys.time()

  # record the item parameter estimation time
  est_time1 <- round(as.numeric(difftime(time2, time1, units = "secs")), 2)

  # the first order test: check convergence-criteria test
  test_1st <- all(c(all(mstep$convergence == 0L), r < MaxE))
  if (test_1st) {
    memo3 <- "Convergence criteria are satisfied."
  } else {
    memo3 <- "Convergence criteria are not satisfied."
    warning(paste0(memo3, " \n"), call. = FALSE)
  }

  # conduct one more E-step to update the posterior distribution
  # using the final item parameter estimates
  estep <- Estep(
    elm_item = elm_item, idx.drm = idx.drm, idx.prm = idx.prm,
    data_drm = data_drm, data_prm = data_prm, data_all = data_all,
    weights = weights, D = D
  )

  # compute the final log of marginal likelihood
  llike <- sum(log(estep$likehd %*% matrix(weights[, 2])))

  ## ---------------------------------------------------------------------
  # estimates the information matrix and standard errors of item parameter estimates
  ## ---------------------------------------------------------------------
  # extract the finalized posterior density
  post_dist <- estep$post_dist

  # delete 'estep' and 'mstep' object
  rm(estep, mstep, data_drm, data_prm, data_all, envir = environment(), inherits = FALSE)

  # compute the information matrix of item parameter estimates using the cross-product method
  if (se) {
    if (verbose) {
      cat("Computing item parameter var-covariance matrix...", "\n")
    }
    time1 <- Sys.time()

    # create a vector of the quadrature points (length = nstd by n.quad)
    quadpt.vec <- rep(quadpt, each = nstd)

    # compute the information matrix of item parameters
    info.data <- info_xpd(
      elm_item = elm_item, freq.cat = freq.cat, post_dist = post_dist,
      quadpt.vec = quadpt.vec, n.quadpt.vec = length(quadpt.vec), nstd = nstd,
      D = D, loc_1p_const = loc_1p_const, loc_else = loc_else, n.1PLM = n.1PLM,
      fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g, a.val.1pl = a.val.1pl,
      a.val.gpcm = a.val.gpcm, g.val = g.val, reloc.par = param_loc$reloc.par
    )

    # compute the information matrix of item parameter priors
    info.prior <- info_prior(
      elm_item = elm_item, D = D, loc_1p_const = loc_1p_const,
      loc_else = loc_else, n.1PLM = n.1PLM, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g,
      a.val.1pl = a.val.1pl, a.val.gpcm = a.val.gpcm, g.val = g.val, aprior = aprior, bprior = bprior,
      gprior = gprior, use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
      reloc.par = param_loc$reloc.par
    )

    # sum of two information matrices
    info.mat <- info.data + info.prior

    # delete 'info.data', 'info.prior', and 'quadpt.vec' objects
    rm(info.data, info.prior, quadpt.vec, envir = environment(), inherits = FALSE)

    # the second-order test: check if the information matrix is positive definite
    test_2nd <- all(eigen(info.mat, only.values = TRUE)$values > 1e-20)
    if (test_2nd) {
      if (test_1st) {
        memo4 <- "Solution is a possible local maximum."
      } else {
        memo4 <- "Information matrix of item parameter estimates is positive definite."
      }
    } else {
      memo4 <- "Information matrix of item parameter estimates is not positive definite; unstable solution."
      warning(paste0(memo4, " \n"), call. = FALSE)
    }

    # compute the variance-covariance matrix of the item parameter estimates, and
    # check if the hessian matrix can be inversed
    cov_mat <- suppressWarnings(tryCatch(
      {
        solve(info.mat, tol = 1e-200)
      },
      error = function(e) {
        NULL
      }
    ))

    # compute the standard errors of item parameter estimates
    if (is.null(cov_mat)) {
      se_par <- rep(99999, length(diag(info.mat)))
      memo5 <- "Variance-covariance matrix of item parameter estimates is not obtainable; unstable solution."
      warning(paste0(memo5, " \n"), call. = FALSE)
    } else {
      se_par <- suppressWarnings(sqrt(diag(cov_mat)))
      memo5 <- "Variance-covariance matrix of item parameter estimates is obtainable."
    }

    # prevent showing NaN values of standard errors
    if (any(is.nan(se_par))) {
      se_par[is.nan(se_par)] <- 99999
    }

    # set an upper bound of standard error
    se_par <- ifelse(se_par > 99999, 99999, se_par)
    time2 <- Sys.time()

    # record the standard error computation time
    est_time2 <- round(as.numeric(difftime(time2, time1, units = "secs")), 2)
  } else {
    memo4 <- "Information matrix of item parameter estimates is not computed."
    memo5 <- "Variance-covariance matrix of item parameter estimates is not computed."
    cov_mat <- NULL
    se_par <- NULL
    est_time2 <- NULL
  }

  # create an item metadata including the final item parameter estimates
  par_df <-
    data.frame(x[, 1:3], elm_item$pars) %>%
    confirm_df(g2na = TRUE)

  # deploy the standard errors on the location of matrix as the item parameter estimates
  se_df <- loc.par <- param_loc$loc.par
  for (i in 1:nrow(loc.par)) {
    num.loc <- which(!is.na(loc.par[i, ]))
    se.loc <- loc.par[i, ][num.loc]
    if (se) {
      se_df[i, num.loc] <- se_par[se.loc]
    } else {
      se_df[i, num.loc] <- NA_real_
    }
  }

  # create a full data.frame for the standard error estimates
  se_df <-
    data.frame(x[, 1:3], se_df) %>%
    confirm_df(g2na = TRUE)

  # create a full data.frame containing the position of item parameter estimates
  # this data.frame is useful when interpreting the variance-covariance matrix of item parameter estimates
  loc_df <-
    data.frame(x[, 1:3], loc.par) %>%
    confirm_df(g2na = TRUE)

  ## ---------------------------------------------------------------------
  # summarize the estimation results
  ## ---------------------------------------------------------------------
  # create a full data.frame including both the item parameter estimates and standard error estimates
  all_df <- data.frame(matrix(NA, nrow = nrow(loc.par), ncol = 2 * ncol(loc.par)))
  all_df[, seq(1, 2 * ncol(loc.par), 2)] <- par_df[, -c(1:3)]
  all_df[, seq(2, 2 * ncol(loc.par), 2)] <- se_df[, -c(1:3)]
  col.names <- rep(NA, 2 * ncol(loc.par))
  col.names[seq(1, 2 * ncol(loc.par), 2)] <- paste0("par.", 1:ncol(loc.par))
  col.names[seq(2, 2 * ncol(loc.par), 2)] <- paste0("se.", 1:ncol(loc.par))
  colnames(all_df) <- col.names
  full_all_df <- data.frame(x[, 1:3], all_df)

  # population density parameters
  moments <- c(mu = group.mean, sigma2 = group.var, sigma = sqrt(group.var))
  moments.se <- rep(NA, 3)

  # data.frame for the population density parameter estimates
  group.par <- data.frame(rbind(moments, moments.se))
  colnames(group.par) <- c("mu", "sigma2", "sigma")
  rownames(group.par) <- c("estimates", "se")

  # prior information
  if (use.aprior) aprior.dist <- aprior else aprior.dist <- NULL
  if (use.bprior) bprior.dist <- bprior else bprior.dist <- NULL
  if (use.gprior) gprior.dist <- gprior else gprior.dist <- NULL

  # statistics based on the loglikelihood of the fitted model:
  npar.est <- length(param_loc$reloc.par)
  neg2llke <- -2 * llike
  aic <- 2 * npar.est + neg2llke
  bic <- npar.est * log(nstd) + neg2llke

  ## ---------------------------------------------------------------
  # check end time
  end.time <- Sys.time()

  # record total computation time
  est_time3 <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 2)

  # return results
  rst <- list(
    estimates = full_all_df, par.est = par_df, se.est = se_df, pos.par = loc_df, covariance = cov_mat, loglikelihood = llike, aic = aic, bic = bic,
    group.par = group.par, weights = weights, posterior.dist = post_dist, data = data, scale.D = D, ncase = nstd, nitem = nrow(par_df),
    Etol = Etol, MaxE = MaxE, aprior = aprior.dist, bprior = bprior.dist, gprior = gprior.dist, npar.est = npar.est, niter = r, maxpar.diff = max.diff,
    EMtime = est_time1, SEtime = est_time2, TotalTime = est_time3, test.1 = memo3, test.2 = memo4, var.note = memo5, fipc = FALSE,
    fipc.method = NULL, fix.loc = NULL
  )

  if (verbose) {
    cat("Estimation is finished in", est_time3, "seconds.", "\n")
  }
  return(rst)
}



# This function implement the fixed item parameter calibration (FIPC) using MMLE-EM algorithm
#' @import dplyr
est_irt_fipc <- function(x = NULL,
                         data,
                         D = 1,
                         item.id = NULL,
                         fix.a.1pl = FALSE,
                         fix.a.gpcm = FALSE,
                         fix.g = FALSE,
                         a.val.1pl = 1,
                         a.val.gpcm = 1,
                         g.val = .2,
                         use.aprior = FALSE,
                         use.bprior = FALSE,
                         use.gprior = TRUE,
                         aprior = list(dist = "lnorm", params = c(0.0, 0.5)),
                         bprior = list(dist = "norm", params = c(0.0, 1.0)),
                         gprior = list(dist = "beta", params = c(5, 16)),
                         missing = NA,
                         Quadrature = c(49, 6.0),
                         weights = NULL,
                         group.mean = 0.0,
                         group.var = 1.0,
                         EmpHist = FALSE,
                         use.startval = FALSE,
                         Etol = 1e-04,
                         MaxE = 500,
                         control = list(eval.max = 200, iter.max = 200),
                         fipc = TRUE,
                         fipc.method = "MEM",
                         fix.loc = NULL,
                         fix.id = NULL,
                         se = TRUE,
                         verbose = TRUE) {

  # check start time
  start.time <- Sys.time()

  ## ---------------------------------------------------------------------
  # prepare the item parameter estimation
  ## ---------------------------------------------------------------------
  # check if item metadata argument of 'x' is provided
  if (is.null(x)) {
    stop(paste0(
      "To implement the fixed item parameter calibration, \n",
      "the item metadata must be specified in the argument 'x'."
    ), call. = FALSE)
  }

  # transform a data set to matrix
  data <- data.matrix(data)

  # extract information about the number of score categories and models
  if (verbose) {
    cat("Parsing input...", "\n")
  }

  # confirm and correct all item metadata information
  x <- confirm_df(x)

  # override the item id in x if item.id is not NULL
  if (!is.null(item.id)) {
    x$id <- item.id
  }

  # check the number of items
  nitem.all <- nrow(x)

  # find the items that should be fixed
  if (is.null(fix.loc) & is.null(fix.id)) {
    stop("When 'FIPC = TRUE', the information of which items are fixed \n",
      "must be provided via the 'fix.loc' or 'fix.id' argument.",
      call. = FALSE
    )
  }
  if (!is.null(fix.loc) & !is.null(fix.id)) {
    warning(paste0(
      "The information given to the 'fix.id' argument was used to fix the item parameters and \n",
      "the information given to the 'fix.loc' argument was ignored."
    ), call. = FALSE)
  }
  if (!is.null(fix.id)) {
    fix.loc <- which(x$id %in% unique(fix.id))
  }

  # check the location of items whose item parameters are estimated
  nofix.loc <- c(1:nitem.all)[!c(1:nitem.all) %in% fix.loc]
  if (length(nofix.loc) == 0L) {
    nofix.loc <- NULL
  }

  # divide the item metadata into two groups: fixed (x_fix) and non-fixed (x_new)
  x_fix <- x[fix.loc, ]
  nitem.fix <- nrow(x_fix)
  if (!is.null(nofix.loc)) {
    x_new <- x[nofix.loc, ]
    nitem.new <- nrow(x_new)
  } else {
    x_new <- NULL
  }

  # clear the two item metadata sets
  x_fix <- confirm_df(x_fix)
  if (!is.null(nofix.loc)) {
    x_new <- confirm_df(x_new)
  }

  # record the score categories and model information of the new items to be estimated
  if (!is.null(x_new)) {
    id <- x_new$id
    cats <- x_new$cats
    model <- x_new$model
  } else {
    id <- x_fix$id
    cats <- x_fix$cats
    model <- x_fix$model
  }

  # generate the empty metadata with starting values
  if (!is.null(x_new) && !use.startval) {
    x_new <- startval_df(cats = cats, model = model, item.id = id)
  }

  # create the total item metadata to be used in the further estimation process
  x_all <- startval_df(cats = x$cats, model = x$model, item.id = x$id)
  x_all[fix.loc, 1:ncol(x_fix)] <- x_fix
  if (!is.null(nofix.loc)) {
    x_all[nofix.loc, 1:ncol(x_new)] <- x_new
  }

  # check whether included data are correct
  # if(nitem.all != ncol(data)) stop("The number of items included in 'x' and 'data' must be the same.", call.=FALSE)

  # recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # check the number of item responses across all items
  n.resp <- Rfast::colsums(!is.na(data))

  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)
  if (length(loc_allmiss) > 0L) {
    memo2 <- paste0(paste0("item ", loc_allmiss, collapse = ", "), " has/have no item response data. \n")
    stop(memo2, call. = FALSE)
  }

  # save the item response data into the different object
  data_all <- data
  rm(data, envir = environment(), inherits = FALSE)

  # divide the item response data into two groups: fixed (x_fix) and non-fixed (x)
  data_fix <- data_all[, fix.loc, drop = FALSE]
  if (!is.null(x_new)) {
    data_new <- data_all[, nofix.loc, drop = FALSE]
  } else {
    data_new <- NULL
  }

  # find the location of 1PLM items in which slope parameters should be constrained to be equal
  # also, find the location of other items
  if ("1PLM" %in% model & !fix.a.1pl) {
    loc_1p_const <- which(model == "1PLM")
    loc_else <- which(model != "1PLM")

    # count the number of 1PLM items to be constrained
    n.1PLM <- length(loc_1p_const)
  } else {
    loc_1p_const <- NULL
    n.1PLM <- NULL
    if (!is.null(x_new)) {
      loc_else <- 1:nrow(x_new)
    } else {
      loc_else <- 1:nrow(x_all)
    }
  }

  # record the original location of item parameters to be estimated, and
  # the relocated position of item parameters when computing
  # the variance-covariance matrix of item parameter estimates
  if (!is.null(x_new)) {
    param_loc <- parloc(
      x = x_new, loc_1p_const = loc_1p_const, loc_else = loc_else,
      fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g
    )
  } else {
    param_loc <- parloc(
      x = x_all, loc_1p_const = loc_1p_const, loc_else = loc_else,
      fix.a.1pl = FALSE, fix.a.gpcm = FALSE, fix.g = FALSE
    )
    param_loc$loc.par[!is.na(param_loc$loc.par)] <- NA_real_
    param_loc$reloc.par <- NULL
  }


  ## ---------------------------------------------------------------------
  # conduct item parameter estimation using MMLE-EM algorithm
  ## ---------------------------------------------------------------------
  # check the total number of examinees
  nstd <- nrow(data_all)

  # create initial weights of prior ability distribution when it is not specified
  if (is.null(weights)) {
    # create quadrature points
    quadpt <- seq(-Quadrature[2], Quadrature[2], length.out = Quadrature[1])
    n.quad <- length(quadpt)

    # create the data.frame containing the quadurature points and weights
    weights <- gen.weight(dist = "norm", mu = group.mean, sigma = sqrt(group.var), theta = quadpt)
  } else {
    quadpt <- weights[, 1]
    n.quad <- length(quadpt)
  }

  # factorize the response values
  if (!is.null(x_new)) {
    resp_new <-
      purrr::map2(
        .x = data.frame(data_new, stringsAsFactors = FALSE), .y = cats,
        .f = function(k, m) factor(k, levels = (seq_len(m) - 1))
      )
    resp_fix <-
      purrr::map2(
        .x = data.frame(data_fix, stringsAsFactors = FALSE), .y = x_fix$cats,
        .f = function(k, m) factor(k, levels = (seq_len(m) - 1))
      )
  } else {
    resp_new <- NULL
    resp_fix <- NULL
  }
  resp_all <-
    purrr::map2(
      .x = data.frame(data_all, stringsAsFactors = FALSE), .y = x_all$cats,
      .f = function(k, m) factor(k, levels = (seq_len(m) - 1))
    )

  # create a contingency table of score categories for each item
  # and then, transform the table to a matrix format
  std.id <- 1:nstd
  if (!is.null(x_new)) {
    freq_new.cat <- purrr::map(
      .x = resp_new,
      .f = function(k) {
        stats::xtabs(~ std.id + k,
          na.action = stats::na.pass, addNA = FALSE
        ) %>%
          matrix(nrow = length(k))
      }
    )
    freq_fix.cat <- purrr::map(
      .x = resp_fix,
      .f = function(k) {
        stats::xtabs(~ std.id + k,
          na.action = stats::na.pass, addNA = FALSE
        ) %>%
          matrix(nrow = length(k))
      }
    )
  } else {
    freq_new.cat <- NULL
    freq_fix.cat <- NULL
  }
  freq_all.cat <- purrr::map(
    .x = resp_all,
    .f = function(k) {
      stats::xtabs(~ std.id + k,
        na.action = stats::na.pass, addNA = FALSE
      ) %>%
        matrix(nrow = length(k))
    }
  )

  # delete 'resp' object
  if (!is.null(x_new)) {
    rm(resp_new, envir = environment(), inherits = FALSE)
    rm(resp_fix, envir = environment(), inherits = FALSE)
  }
  rm(resp_all, envir = environment(), inherits = FALSE)

  # break down the item metadata into several elements
  if (!is.null(x_new)) {
    elm_item_new <- breakdown(x_new)
    elm_item_fix <- breakdown(x_fix)
  }
  elm_item_all <- breakdown(x_all)

  # classify the items into DRM and PRM item groups
  if (!is.null(x_new)) {
    idx.item.new <- idxfinder(elm_item_new)
    idx.item.fix <- idxfinder(elm_item_fix)
    idx.drm.new <- idx.item.new$idx.drm
    idx.prm.new <- idx.item.new$idx.prm
    idx.drm.fix <- idx.item.fix$idx.drm
    idx.prm.fix <- idx.item.fix$idx.prm
  }
  idx.item.all <- idxfinder(elm_item_all)
  idx.drm.all <- idx.item.all$idx.drm
  idx.prm.all <- idx.item.all$idx.prm

  # then, divide items in the new test form into three groups
  # (1) DRM 1PL items with the constrained slop
  # (2) other DRM items
  # (3) PRM items
  if (!is.null(x_new)) {
    if (sum(loc_else) == 0) {
      drm.else <- NULL
    } else {
      drm.else <- loc_else
    }
    idx4est <- list(
      drm.slc = loc_1p_const,
      drm.else = drm.else,
      prm = idx.prm.new
    )
  } else {
    idx4est <- NULL
  }

  # divide the data set for the mixed-item format
  if (!is.null(x_new)) {
    datlist.new <- divide_data(data = data_new, idx.item = idx.item.new, freq.cat = freq_new.cat)
    datlist.fix <- divide_data(data = data_fix, idx.item = idx.item.fix, freq.cat = freq_fix.cat)
    data_drm.new <- cbind(datlist.new$data_drm_q, datlist.new$data_drm_p)
    data_drm.fix <- cbind(datlist.fix$data_drm_q, datlist.fix$data_drm_p)
    data_prm.new <- datlist.new$data_prm
    data_prm.fix <- datlist.fix$data_prm
    data.new <- datlist.new$data_all
    data.fix <- datlist.fix$data_all
  }
  datlist.all <- divide_data(data = data_all, idx.item = idx.item.all, freq.cat = freq_all.cat)
  data_drm.all <- cbind(datlist.all$data_drm_q, datlist.all$data_drm_p)
  data_prm.all <- datlist.all$data_prm
  data.all <- datlist.all$data_all

  # delete 'datlist' object
  if (!is.null(x_new)) {
    rm(datlist.new, datlist.fix, freq_fix.cat, envir = environment(), inherits = FALSE)
  }
  rm(datlist.all, freq_all.cat, envir = environment(), inherits = FALSE)

  # moments of population prior distribution
  if (is.null(x_new)) {
    mmt_dist_old <- c(group.mean, group.var)
  }

  # create the lower and upper bounds of the item parameters
  if (!is.null(x_new)) {
    parbd <- lubound(model, cats, n.1PLM, idx4est, fix.a.1pl, fix.g, fix.a.gpcm)
  } else {
    parbd <- NULL
  }

  # find the columns of the frequency matrix corresponding to all items
  if (!is.null(x_new)) {
    cols.item <- cols4item(nitem.new, cats, loc_1p_const)
  } else {
    cols.item <- NULL
  }

  # estimation
  if (verbose) {
    cat("Estimating item parameters...", "\n")
  }

  # set the number of EM iteration to one when OEM (Wainer & Mislevy, 1990) method is used
  if (fipc.method == "OEM") MaxE <- 1

  # implement EM algorithm
  time1 <- Sys.time()
  for (r in 1:MaxE) {
    # implement E-step
    if (!is.null(x_new)) {
      if (r == 1L) {
        estep <- Estep_fipc(
          elm_item1 = elm_item_new, elm_item2 = elm_item_fix, idx.drm2 = idx.drm.fix,
          idx.prm2 = idx.prm.fix, data_drm2 = data_drm.fix, data_prm2 = data_prm.fix,
          data_all1 = data.new, weights = weights, D = D
        )
      } else {
        estep <- Estep_fipc(
          elm_item1 = elm_item_new, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
          idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
          data_all1 = data.new, weights = weights, D = D
        )
      }
    } else {
      estep <- Estep_fipc(
        elm_item1 = elm_item_all, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
        idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
        data_all1 = data.all, weights = weights, D = D
      )
      estep$elm_item <- NULL
    }

    # implement M-step
    mstep <- Mstep(
      estep = estep, id = id, cats = cats, model = model, quadpt = quadpt, n.quad = n.quad,
      D = D, cols.item = cols.item, loc_1p_const = loc_1p_const, loc_else = loc_else, idx4est = idx4est,
      n.1PLM = n.1PLM, EmpHist = EmpHist, weights = weights, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g,
      a.val.1pl = a.val.1pl, a.val.gpcm = a.val.gpcm, g.val = g.val, use.aprior = use.aprior, use.bprior = use.bprior,
      use.gprior = use.gprior, aprior = aprior, bprior = bprior, gprior = gprior, group.mean = group.mean,
      group.var = group.var, nstd = nstd, Quadrature = Quadrature, control = control,
      iter = r, fipc = TRUE, reloc.par = param_loc$reloc.par, parbd = parbd
    )

    if (!is.null(x_new)) {
      # compute the difference between previous and updated item parameter estimates
      diff_par <- mstep$elm_item$pars - elm_item_new$pars
      max.diff <- abs(max(diff_par, na.rm = TRUE))
    } else {
      # compute the mean and sd of the updated prior distribution
      mmt_dist_new <- cal_moment(node = mstep$weights$theta, weight = mstep$weights$weight)
      diff_par <- mmt_dist_new - mmt_dist_old
      max.diff <- abs(max(diff_par, na.rm = TRUE))
    }

    # log-likelihood value
    llike <- mstep$loglike

    # print
    if (verbose) {
      cat("\r", paste0(
        "EM iteration: ", r, ", Loglike: ", format(round(mstep$loglike, 4), nsmall = 4),
        ", Max-Change: ", format(round(max.diff, 6), nsmall = 5)
      ))
    }

    # check the convergence of EM algorithm
    converge <- max.diff <= Etol

    # extract the updated item (or group) parameter estimates
    # and update the new and all item parameters
    if (!is.null(x_new)) {
      elm_item_new$pars <- mstep$elm_item$pars
      elm_item_all$pars[nofix.loc, 1:ncol(elm_item_new$pars)] <- elm_item_new$pars
    } else {
      mmt_dist_old <- mmt_dist_new
    }

    # extract the updated quadrature points and
    # the corresponding weights of the prior population density
    weights <- mstep$weights

    # terminate the EM step if the convergence criterion is satisfied
    if (converge | r == MaxE) {
      break
    }
  }

  if (verbose) {
    cat("", "\n")
  }
  time2 <- Sys.time()

  # record the item parameter estimation time
  est_time1 <- round(as.numeric(difftime(time2, time1, units = "secs")), 2)

  # the first order test: check convergence-criteria test
  test_1st <- all(c(all(mstep$convergence == 0L), r < MaxE))
  if (test_1st) {
    memo3 <- "Convergence criteria are satisfied."
  } else {
    memo3 <- "Convergence criteria are not satisfied."
    warning(paste0(memo3, " \n"), call. = FALSE)
  }

  # conduct one more E-step to update the posterior distribution
  # using the final item parameter estimates
  if (!is.null(x_new)) {
    estep <- Estep_fipc(
      elm_item1 = elm_item_new, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
      idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
      data_all1 = data.new, weights = weights, D = D
    )
  } else {
    estep <- Estep_fipc(
      elm_item1 = elm_item_all, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
      idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
      data_all1 = data.all, weights = weights, D = D
    )
  }

  # compute the final log of marginal likelihood
  llike <- sum(log(estep$likehd %*% matrix(weights[, 2])))

  # compute the mean and variance of the estimated density distribution
  pop_moments <- cal_moment(node = quadpt, weight = weights[, 2])

  ## ---------------------------------------------------------------------
  # estimates the standard errors of item parameter estimates
  ## ---------------------------------------------------------------------
  if (!is.null(x_new)) {
    # extract the finalized posterior density
    post_dist <- estep$post_dist

    # delete 'estep' and 'mstep' object
    rm(estep, mstep, data_drm.new, data_drm.fix,
      data_prm.new, data_prm.fix, data.new, data.fix,
      data_drm.all, data_prm.all, data.all,
      envir = environment(), inherits = FALSE
    )

    # compute the information matrix of item parameter estimates using the cross-product method
    if (se) {
      if (verbose) {
        cat("Computing item parameter var-covariance matrix...", "\n")
      }
      time1 <- Sys.time()

      # create a vector of the quadrature points (length = nstd by n.quad)
      quadpt.vec <- rep(quadpt, each = nstd)

      # compute the information matrix of item parameters
      info.data <- info_xpd(
        elm_item = elm_item_new, freq.cat = freq_new.cat, post_dist = post_dist,
        quadpt.vec = quadpt.vec, n.quadpt.vec = length(quadpt.vec), nstd = nstd,
        D = D, loc_1p_const = loc_1p_const, loc_else = loc_else, n.1PLM = n.1PLM,
        fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g, a.val.1pl = a.val.1pl,
        a.val.gpcm = a.val.gpcm, g.val = g.val, reloc.par = param_loc$reloc.par
      )

      # compute the information matrix of item parameter priors
      info.prior <- info_prior(
        elm_item = elm_item_new, D = D, loc_1p_const = loc_1p_const,
        loc_else = loc_else, n.1PLM = n.1PLM, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g,
        a.val.1pl = a.val.1pl, a.val.gpcm = a.val.gpcm, g.val = g.val, aprior = aprior, bprior = bprior,
        gprior = gprior, use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
        reloc.par = param_loc$reloc.par
      )

      # sum of two information matrices
      info.mat <- info.data + info.prior

      # the second-order test: check if the information matrix is positive definite
      test_2nd <- all(eigen(info.mat, only.values = TRUE)$values > 1e-20)
      if (test_2nd) {
        if (test_1st) {
          memo4 <- "Solution is a possible local maximum."
        } else {
          memo4 <- "Information matrix of item parameter estimates is positive definite."
        }
      } else {
        memo4 <- "Information matrix of item parameter estimates is not positive definite; unstable solution."
        warning(paste0(memo4, " \n"), call. = FALSE)
      }

      # compute the variance-covariance matrix of the item parameter estimates, and
      # check if the hessian matrix can be inversed
      cov_mat <- suppressWarnings(tryCatch(
        {
          solve(info.mat, tol = 1e-200)
        },
        error = function(e) {
          NULL
        }
      ))

      # compute the standard errors of item parameter estimates
      if (is.null(cov_mat)) {
        se_par <- rep(99999, length(diag(info.mat)))
        memo5 <- "Variance-covariance matrix of item parameter estimates is not obtainable; unstable solution."
        warning(paste0(memo5, " \n"), call. = FALSE)
      } else {
        se_par <- suppressWarnings(sqrt(diag(cov_mat)))
        memo5 <- "Variance-covariance matrix of item parameter estimates is obtainable."
      }

      # prevent showing NaN values of standard errors
      if (any(is.nan(se_par))) {
        se_par[is.nan(se_par)] <- 99999
      }

      # set an upper bound of standard error
      se_par <- ifelse(se_par > 99999, 99999, se_par)
      time2 <- Sys.time()

      # record the standard error computation time
      est_time2 <- round(as.numeric(difftime(time2, time1, units = "secs")), 2)
    } else {
      memo4 <- "Information matrix of item parameter estimates is not computed."
      memo5 <- "Variance-covariance matrix of item parameter estimates is not computed."
      cov_mat <- NULL
      se_par <- NULL
      est_time2 <- NULL
    }

    # create an item metadata including the final item parameter estimates for all items
    x_all <-
      data.frame(x_all[, 1:3], elm_item_all$pars) %>%
      confirm_df(g2na = TRUE)

    # deploy the standard errors on the location of matrix as the item parameter estimates
    # 1) for the only new items
    se_df <- loc.par <- param_loc$loc.par
    for (i in 1:nrow(loc.par)) {
      num.loc <- which(!is.na(loc.par[i, ]))
      se.loc <- loc.par[i, ][num.loc]
      if (se) {
        se_df[i, num.loc] <- se_par[se.loc]
      } else {
        se_df[i, num.loc] <- NA_real_
      }
    }

    # 2) for the a total test form
    all.col <- max(x_all$cats)
    all.col <- ifelse(all.col == 2, 3, all.col)
    se_all_df <- loc_all.par <- matrix(NA, nrow = nitem.all, ncol = all.col)
    if (ncol(se_df) < ncol(se_all_df)) {
      n2add <- ncol(se_all_df) - ncol(se_df)
      se_df <- cbind(se_df, matrix(NA, nrow = nrow(se_df), ncol = n2add))
    }
    se_all_df[nofix.loc, 1:ncol(se_df)] <- se_df
    loc_all.par[nofix.loc, 1:ncol(loc.par)] <- loc.par

    # create a full data.frame for the standard error estimates
    se_all_df <-
      data.frame(x_all[, 1:3], se_all_df) %>%
      confirm_df(g2na = TRUE)

    # create a full data.frame containing the position of item parameter estimates
    # this data.frame is useful when interpreting the variance-covariance matrix of item parameter estimates
    loc_all_df <-
      data.frame(x_all[, 1:3], loc_all.par) %>%
      confirm_df(g2na = TRUE)
  } else {
    # extract the finalized posterior density
    post_dist <- estep$post_dist

    # the second-order test
    if (test_1st) {
      memo4 <- "Solution is a possible local maximum."
    } else {
      memo4 <- "Solution is not a possible local maximum because convergence criteria are not satisfied."
    }

    # compute the standard errors of item parameter estimates
    cov_mat <- NULL
    memo5 <- "Variance-covariance matrix of item parameter estimates was not estimated."

    # record the standard error computation time (NULL)
    est_time2 <- NULL

    # deploy the standard errors on the location of matrix as the item parameter estimates
    se_df <- loc_all.par <- param_loc$loc.par

    # create a full data.frame for the standard error estimates
    # however, SEs are all NA when only group parameters are estimated
    se_all_df <-
      data.frame(x_all[, 1:3], se_df) %>%
      confirm_df(g2na = TRUE)

    # create a full data.frame containing the position of item parameter estimates
    # however, this should be NULL when only group parameters are estimated
    loc_all_df <- NULL
  }

  ## ---------------------------------------------------------------------
  # summarize the estimation results
  ## ---------------------------------------------------------------------
  # create a full data.frame including both the item parameter estimates and standard error estimates
  all_df <- data.frame(matrix(NA, nrow = nrow(loc_all.par), ncol = 2 * ncol(loc_all.par)))
  all_df[, seq(1, 2 * ncol(loc_all.par), 2)] <- x_all[, -c(1:3)]
  all_df[, seq(2, 2 * ncol(loc_all.par), 2)] <- se_all_df[, -c(1:3)]
  col.names <- rep(NA, 2 * ncol(loc_all.par))
  col.names[seq(1, 2 * ncol(loc_all.par), 2)] <- paste0("par.", 1:ncol(loc_all.par))
  col.names[seq(2, 2 * ncol(loc_all.par), 2)] <- paste0("se.", 1:ncol(loc_all.par))
  colnames(all_df) <- col.names
  full_all_df <- data.frame(x_all[, 1:3], all_df)

  # population density parameters
  mu <- pop_moments[1]
  sigma2 <- pop_moments[2]
  sigma <- sqrt(pop_moments[2])
  moments.est <- c(mu, sigma2, sigma)

  # compute the standard errors of population density parameters
  se.mu <- sigma / sqrt(nstd)
  se.sigma2 <- sigma2 * sqrt(2 / (nstd - 1))
  se.sigma <- (1 / (2 * sigma)) * se.sigma2 # using Delta method
  moments.se <- c(se.mu, se.sigma2, se.sigma)

  # data.frame for the population density parameter estimates
  group.par <- data.frame(rbind(moments.est, moments.se))
  colnames(group.par) <- c("mu", "sigma2", "sigma")
  rownames(group.par) <- c("estimates", "se")

  # prior information
  if (use.aprior) aprior.dist <- aprior else aprior.dist <- NULL
  if (use.bprior) bprior.dist <- bprior else bprior.dist <- NULL
  if (use.gprior) gprior.dist <- gprior else gprior.dist <- NULL

  # statistics based on the loglikelihood of the fitted model:
  if (!is.null(x_new)) {
    npar.est <- length(param_loc$reloc.par) + 2
  } else {
    npar.est <- 2
  }
  neg2llke <- -2 * llike
  aic <- 2 * npar.est + neg2llke
  bic <- npar.est * log(nstd) + neg2llke

  ## ---------------------------------------------------------------
  # check end time
  end.time <- Sys.time()

  # record total computation time
  est_time3 <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 2)

  # return results
  rst <- list(
    estimates = full_all_df, par.est = x_all, se.est = se_all_df, pos.par = loc_all_df, covariance = cov_mat, loglikelihood = llike,
    aic = aic, bic = bic, group.par = group.par, weights = weights, posterior.dist = post_dist, data = data_all, scale.D = D, ncase = nstd,
    nitem = nitem.all, Etol = Etol, MaxE = MaxE, aprior = aprior.dist, bprior = bprior.dist, gprior = gprior.dist, npar.est = npar.est,
    niter = r, maxpar.diff = max.diff, EMtime = est_time1, SEtime = est_time2, TotalTime = est_time3, test.1 = memo3, test.2 = memo4,
    var.note = memo5, fipc = TRUE, fipc.method = fipc.method, fix.loc = fix.loc
  )

  if (verbose) {
    cat("Estimation is finished in", est_time3, "seconds.", "\n")
  }
  return(rst)
}
