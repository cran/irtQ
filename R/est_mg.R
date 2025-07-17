#' Multiple-group item calibration using MMLE-EM algorithm
#'
#' This function performs multiple-group item calibration (Bock & Zimowski, 1997)
#' using marginal maximum likelihood estimation via the expectation-maximization
#' (MMLE-EM) algorithm (Bock & Aitkin, 1981). It also supports multiple-group
#' fixed item parameter calibration (MG-FIPC; e.g., Kim & Kolen, 2016), which
#' extends the single-group FIPC method (Kim, 2006) to multiple-group settings.
#' For dichotomous items, the function supports one-, two-, and three-parameter
#' logistic IRT models. For polytomous items, the graded response model (GRM)
#' and the (generalized) partial credit model (GPCM) are available.
#'
#' @inheritParams est_irt
#' @param x A list containing item metadata for all groups to be analyzed.
#'   For example, if five groups are analyzed, the list should contain five
#'   elements, each representing the item metadata for one group. The order of
#'   the elements in the list must match the order of group names specified in
#'   the `group.name` argument.
#'
#'   Each group's item metadata includes essential information for each item
#'   (e.g., number of score categories, IRT model type, etc.) required for
#'   calibration. See [irtQ::est_irt()] or [irtQ::simdat()] for more
#'   details about the item metadata.
#'
#'   When `use.startval = TRUE`, the item parameters specified in the metadata
#'   will be used as starting values for parameter estimation. If `x = NULL`,
#'   both `model` and `cats` arguments must be specified. Note that when
#'   `fipc = TRUE` to implement MG-FIPC, the `x` argument must be specified
#'   and cannot be NULL. Default is `NULL`.
#' @param data A list containing item response matrices for all groups to be
#'   analyzed. For example, if five groups are analyzed, the list should include
#'   five elements, each representing the response data matrix for one group.
#'   The elements in the list must be ordered to match the group names specified
#'   in the `group.name` argument. Each matrix contains examinees' item responses
#'   corresponding to the item metadata for that group. In each matrix, rows
#'   represent examinees and columns represent items.
#' @param group.name A character vector indicating the names of the groups.
#'   For example, if five groups are analyzed, use
#'   `group.name = c("G1", "G2", "G3", "G4", "G5")`. Group names can be any
#'   valid character strings.
#' @param model A list containing character vectors specifying the IRT models
#'   used to calibrate items across all groups. For example, if five groups are
#'   analyzed, the list should contain five elements, each being a character
#'   vector of IRT model names for one group. The elements in the list must be
#'   ordered to match the group names specified in the `group.name` argument.
#'
#'   Available IRT models include:
#'   - `"1PLM"`, `"2PLM"`, `"3PLM"`, `"DRM"` for dichotomous items
#'   - `"GRM"`, `"GPCM"` for polytomous items
#'
#'   Here, `"GRM"` denotes the graded response model and `"GPCM"` the
#'   (generalized) partial credit model. Note that `"DRM"` serves as a general
#'   label covering all three dichotomous IRT models.If a single model name is
#'   provided in any element of the list, it will be recycled across all items
#'   within that group.
#'
#'   This argument is used only when `x = NULL` and `fipc = FALSE`. Default
#'   is `NULL`.
#' @param cats A list containing numeric vectors specifying the number of score
#'   categories for items in each group. For example, if five groups are
#'   analyzed, the list should contain five numeric vectors corresponding to
#'   the five groups. The elements in the list must be ordered consistently with
#'   the group names specified in the `group.name` argument.
#'
#'   If a single numeric value is specified in any element of the list, it will
#'   be recycled across all items in the corresponding group. If `cats = NULL`
#'   and all models specified in the `model` argument are dichotomous
#'   (i.e., "1PLM", "2PLM", "3PLM", or "DRM"), the function assumes that all
#'   items have two score categories across all groups.
#'
#'   This argument is used only when `x = NULL` and `fipc = FALSE`. Default
#'   is `NULL`.
#' @param item.id A list containing character vectors of item IDs for each group
#'   to be analyzed. For example, if five groups are analyzed, the list should
#'   contain five character vectors corresponding to the five groups. The elements
#'   in the list must be ordered consistently with the group names specified in
#'   the `group.name` argument.
#'
#'   When `fipc = TRUE` and item IDs are provided via the `item.id` argument,
#'   the item IDs in the `x` argument will be overridden. Default is `NULL`.
#' @param free.group A numeric or character vector indicating the groups for which
#'   the scales (i.e., mean and standard deviation) of the latent ability distributions
#'   are freely estimated. The scales of the remaining groups (those not specified in this
#'   argument) are fixed using the values provided in the `group.mean` and `group.var`
#'   arguments, or from the `weights` argument.
#'
#'   For example, suppose that five groups are analyzed with group names "G1", "G2",
#'   "G3", "G4", and "G5". To freely estimate the scales for groups 2 through 5, set
#'   `free.group = c(2, 3, 4, 5)` or `free.group = c("G2", "G3", "G4", "G5")`.
#'   In this case, the first group ("G1") will have a fixed scale
#'   (e.g., a mean of 0 and variance of 1 when `group.mean = 0`, `group.var = 1`,
#'   and `weights = NULL`).
#' @param weights A two-column matrix or data frame containing the quadrature
#'   points (in the first column) and the corresponding weights (in the second
#'   column) for the latent ability prior distribution. If not `NULL`, the
#'   latent ability distributions for the groups not specified in the `free.group`
#'   argument are fixed to match the scale defined by the provided quadrature points
#'   and weights. The weights and points can be conveniently generated using
#'   the function [irtQ::gen.weight()].
#'
#'   If `NULL`, a normal prior density is used instead, based on the
#'   information provided in the `Quadrature`, `group.mean`, and `group.var`
#'   arguments. Default is `NULL`.
#' @param group.mean A numeric value specifying the mean of the latent variable
#'   prior distribution when `weights = NULL`. Default is 0. For groups not
#'   specified in the `free.group` argument, their distribution means are fixed
#'   to this value in order to resolve the indeterminacy of the item parameter
#'   scale.
#' @param group.var A positive numeric value specifying the variance of the
#'   latent variable prior distribution when `weights = NULL`. Default is 1.
#'   For groups not specified in the `free.group` argument, their distribution
#'   variances are fixed to this value in order to resolve the indeterminacy of
#'   the item parameter scale.
#' @param EmpHist Logical. If `TRUE`, the empirical histograms of the latent
#'   ability prior distributions across all groups are estimated simultaneously
#'   with the item parameters using the approach proposed by Woods (2007).
#'   Item calibration is then performed relative to the estimated empirical priors.
#' @param Etol A positive numeric value specifying the convergence criterion for
#'   the E-step of the EM algorithm. Default is `1e-4`.
#' @param fipc Logical. If `TRUE`, multiple-group fixed item parameter
#'   calibration (MG-FIPC) is applied during item parameter estimation.
#'   When `fipc = TRUE`, the information on which items are fixed
#'   must be provided via either `fix.loc` or `fix.id`. See below for details.
#' @param fix.loc A list of positive integer vectors. Each internal vector specifies
#'   the positions of the items to be fixed in the item metadata (i.e., `x`)
#'   for each group when MG-FIPC is implemented (i.e., `fipc = TRUE`).
#'   The internal objects in the list must follow the same order as the group names
#'   provided in the `group.name` argument.
#'
#'   For example, suppose three groups are analyzed. In the first group,
#'   the 1st, 3rd, and 5th items are fixed; in the second group, the 2nd, 3rd,
#'   4th, and 7th items are fixed; and in the third group, the 1st, 2nd, and
#'   6th items are fixed. Then
#'   `fix.loc = list(c(1, 3, 5), c(2, 3, 4, 7), c(1, 2, 6))`.
#'   Note that if the `fix.id` argument is not NULL, the information in `fix.loc`
#'   will be ignored. See below for details.
#' @param fix.id A vector of character strings specifying the IDs of items to be
#'   fixed when MG-FIPC is implemented (i.e., `fipc = TRUE`).
#'
#'   For example, suppose that three groups are analyzed. In the first group,
#'   three items with IDs G1I1, C1I1, and C1I2 are fixed. In the second group,
#'   four items with IDs C1I1, C1I2, C2I1, and C2I2 are fixed. In the third
#'   group, three items with IDs C2I1, C2I2, and G3I1 are fixed.
#'
#'   In this case, there are six unique items fixed across the groups—namely,
#'   G1I1, C1I1, C1I2, C2I1, C2I2, and G3I1, because C1I1 and C1I2 appear in
#'   both the first and second groups, while C2I1 and C2I2 appear in both the
#'   second and third groups. Thus, you should specify
#'   `fix.id = c("G1I1", "C1I1", "C1I2", "C2I1", "C2I2", "G3I1")`.
#'   Note that if the `fix.id` argument is not NULL, the information provided in
#'   `fix.loc` is ignored. See below for details.
#'
#' @details Multiple-group (MG) item calibration (Bock & Zimowski, 1996)
#' provides a unified framework for handling testing scenarios involving
#' multiple groups, such as nonequivalent groups equating, vertical scaling,
#' and the identification of differential item functioning (DIF). In such
#' applications, examinees from different groups typically respond to either
#' the same test form or to different forms that share common (anchor) items.
#'
#' The goal of MG item calibration is to estimate both item parameters and
#' latent ability distributions for all groups simultaneously (Bock & Zimowski, 1996).
#' The \pkg{irtQ} package implements MG calibration via the [irtQ::est_mg()] function,
#' which uses marginal maximum likelihood estimation through the
#' expectation-maximization (MMLE-EM) algorithm (Bock & Aitkin, 1981).
#' In addition, the function supports multiple-group fixed item parameter
#' calibration (MG-FIPC; e.g., Kim & Kolen, 2016), which allows the parameters
#' of specific items to be fixed across groups.
#'
#' In MG IRT analyses, it is common for multiple groups' test forms to share
#' some common (anchor) items. By default, the [irtQ::est_mg()] function
#' automatically constrains items with identical item IDs across groups
#' to share the same parameter estimates.
#'
#' Most of the features of the [irtQ::est_mg()] function are similar to those of the
#' [irtQ::est_irt()] function. The main difference is that several arguments in
#' [irtQ::est_mg()] accept list objects containing elements for each group to be
#' analyzed. These arguments include `x`, `data`, `model`, `cats`,
#' `item.id`, `fix.loc` and `fix.id`.
#'
#' Additionally, [irtQ::est_mg()] introduces two new arguments: `group.name`
#' and `free.group`. The `group.name` argument is required to assign a unique
#' identifier to each group. The order of the list elements provided in
#' `x`, `data`, `model`, `cats`, `item.id`, `fix.loc` and `fix.id`
#' must match the order of group names specified in the `group.name` argument.
#'
#' The `free.group` argument is required to indicate which groups have their
#' latent ability distribution scales (i.e., mean and standard deviation)
#' freely estimated. When no item parameters are fixed (i.e., `fipc = FALSE`),
#' at least one group must have a fixed latent ability scale (e.g., mean = 0 and
#' variance = 1) among the multiple groups sharing common items, in order to
#' resolve the scale indeterminacy inherent in IRT estimation. By specifying
#' the groups in the `free.group` argument, the scales for those groups will be
#' freely estimated, while the scales for all other groups not included in
#' `free.group` will be fixed using the values provided in the `group.mean` and
#' `group.var` arguments or from the `weights` argument.
#'
#' Situations requiring the implementation of MG-FIPC typically arise when
#' new latent ability scales from multiple-group (MG) test data need to be
#' linked to an established scale (e.g., that of an existing item bank).
#' In a single run of the MG-FIPC procedure, the parameters of non-fixed
#' (freed) items across multiple test forms, as well as the latent ability
#' distributions for multiple groups, can be estimated on the same scale
#' as the fixed items (Kim & Kolen, 2016).
#'
#' For example, suppose that three different test forms—Form 1, Form 2, and
#' Form 3—are administered to three nonequivalent groups: Group1, Group2,
#' and Group3. Form 1 and Form 2 share 12 common items (C1I1 to C1I12),
#' while Form 2 and Form 3 share 10 common items (C2I1 to C2I10). There are
#' no common items between Form 1 and Form 3. Also, assume that all unique
#' items in Form 1 are from an existing item bank and have already been
#' calibrated on the item bank's scale.
#'
#' In this case, the goal of MG-FIPC is to estimate the parameters of all
#' items across the three test forms—except the unique items in Form 1—
#' and the latent ability distributions of the three groups, all on the
#' same scale as the item bank. To achieve this, the unique items in Form 1
#' must be fixed during MG-FIPC to link the current MG test data to the
#' item bank scale.
#'
#' The [irtQ::est_mg()] function can implement MG-FIPC by setting
#' `fipc = TRUE`. In this case, the information on which items to fix must be
#' provided through either the `fix.loc` or `fix.id` argument. When using
#' `fix.loc`, you must supply a list of item positions (locations)
#' to be fixed in each group’s test form. For example, suppose that the test
#' data from the three groups above are analyzed. In the first group,
#' the 1st, 3rd, and 5th items are fixed; in the second group, the 2nd, 3rd,
#' 4th, and 7th items are fixed; and in the third group, the 1st, 2nd, and
#' 6th items are fixed. In this case, you would specify:
#' `fix.loc = list(c(1, 3, 5), c(2, 3, 4, 7), c(1, 2, 6))`.
#'
#' Alternatively, you can use `fix.id` to specify a character vector of item IDs
#' to be fixed across groups. In the first group, the items with IDs G1I1, C1I1,
#' and C1I2 are fixed; in the second group, the items with IDs C1I1, C1I2,
#' C2I1, and C2I2 are fixed; and in the third group, the items with IDs C2I1,
#' C2I2, and G3I1 are fixed. In this case, there are six unique items to be
#' fixed across all groups: G1I1, C1I1, C1I2, C2I1, C2I2, and G3I1. You would
#' then specify:
#' `fix.id = c("G1I1", "C1I1", "C1I2", "C2I1", "C2I2", "G3I1")`.
#'
#' Note that when both `fix.loc` and `fix.id` are provided, the information
#' in `fix.id` takes precedence and overrides `fix.loc`.
#'
#' @return This function returns an object of class `est_irt`. The returned
#'   object contains the following components:
#'
#' \item{estimates}{A list containing two internal elements: `overall` and `group`.
#'   The `overall` element is a data frame with item parameter estimates and their
#'   standard errors, computed from the combined data across all groups. This data
#'   frame includes only the unique items across all groups. The `group` element is
#'   a list of group-specific data frames, each containing item parameter estimates
#'   and standard errors for that particular group.}
#'
#' \item{par.est}{A list with the same structure as `estimates`, containing only
#'   the item parameter estimates (excluding standard errors), formatted according to
#'   the item metadata structure.}
#'
#' \item{se.est}{A list with the same structure as `estimates`, but containing only
#'   the standard errors of the item parameter estimates. Note that the standard
#'   errors are calculated using the cross-product approximation method
#'   (Meilijson, 1989).}
#'
#' \item{pos.par}{A data frame indicating the position index of each estimated
#'   item parameter. This index is based on the combined data set across all groups
#'   (i.e., the first internal object of `estimates`). The position information is
#'   useful for interpreting the variance-covariance matrix of item parameter
#'   estimates.}
#'
#' \item{covariance}{A variance-covariance matrix of the item parameter
#'   estimates, based on the combined data set across all groups
#'   (i.e., the first internal object of `estimates`).}
#'
#' \item{loglikelihood}{A list containing two internal objects (i.e., overall
#'   and group) of marginal log-likelihood values based on the observed data.
#'   The structure of the list matches that of `estimates`. Specifically,
#'   the `overall` component contains the total log-likelihood summed across
#'   all unique items from all groups, while the `group` component provides
#'   group-specific log-likelihood values.}
#'
#' \item{aic}{A model fit statistic based on the Akaike Information Criterion (AIC),
#'   calculated from the log-likelihood of all unique items.}
#'
#' \item{bic}{A model fit statistic based on the Bayesian Information Criterion (BIC),
#'   calculated from the log-likelihood of all unique items.}
#'
#' \item{group.par}{A list containing summary statistics (i.e., mean, variance,
#'   and standard deviation) of the latent variable prior distributions across
#'   all groups.}
#'
#' \item{weights}{A list of two-column data frames, where the first column
#'   contains quadrature points and the second column contains the corresponding
#'   weights of the (updated) latent variable prior distributions for each group.}
#'
#' \item{posterior.dist}{A matrix of normalized posterior densities for all
#'   response patterns at each quadrature point. Rows and columns represent
#'   response patterns and quadrature points, respectively.}
#'
#' \item{data}{A list containing two internal objects (i.e., overall and group)
#'   representing the examinees' response data sets. The structure of this
#'   list matches that of the `estimates` component.}
#'
#' \item{scale.D}{The scaling factor used in the IRT model.}
#'
#' \item{ncase}{A list containing two internal objects (i.e., overall and group)
#'   representing the total number of response patterns. The structure of this
#'   list matches that of the `estimates` component.}
#'
#' \item{nitem}{A list containing two internal objects (i.e., overall and group)
#'   representing the total number of items included in the response data.
#'   The structure of this list matches that of the `estimates` component.}
#'
#' \item{Etol}{The convergence criterion for the E-step of the EM algorithm.}
#'
#' \item{MaxE}{The maximum number of E-steps allowed in the EM algorithm.}
#'
#' \item{aprior}{A list describing the prior distribution used for discrimination
#'    parameters.}
#'
#' \item{bprior}{A list describing the prior distribution used for difficulty
#'    parameters.}
#'
#' \item{gprior}{A list describing the prior distribution used for guessing
#'    parameters.}
#'
#' \item{npar.est}{The total number of parameters estimated across all unique items.}
#'
#' \item{niter}{The number of completed EM cycles.}
#'
#' \item{maxpar.diff}{The maximum absolute change in parameter estimates at
#'   convergence.}
#'
#' \item{EMtime}{Time (in seconds) spent on EM cycles.}
#'
#' \item{SEtime}{Time (in seconds) spent computing standard errors.}
#'
#' \item{TotalTime}{Total computation time (in seconds).}
#'
#' \item{test.1}{First-order test result indicating whether the gradient
#'   sufficiently vanished for solution stability.}
#'
#' \item{test.2}{Second-order test result indicating whether the information matrix
#'   is positive definite, a necessary condition for identifying a local maximum.}
#'
#' \item{var.note}{A note indicating whether the variance-covariance matrix
#'   was successfully obtained from the information matrix.}
#'
#' \item{fipc}{Logical. Indicates whether FIPC was used.}
#'
#' \item{fipc.method}{The method used for FIPC.}
#'
#' \item{fix.loc}{A list containing two internal objects (i.e., overall and group)
#'   indicating the locations of fixed items when FIPC is applied. The structure
#'   of the list matches that of the 'estimates' component.}
#'
#'   Note that you can easily extract components from the output using the
#'   [irtQ::getirt()] function.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso [irtQ::est_irt()], [irtQ::shape_df()], [irtQ::shape_df_fipc()],
#' [irtQ::getirt()]
#'
#' @references Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood
#' estimation of item parameters: Application of an EM algorithm.
#' *Psychometrika, 46*, 443-459.
#'
#' Bock, R. D., & Zimowski, M. F. (1997). Multiple group IRT. In W. J.
#' van der Linden & R. K. Hambleton (Eds.), *Handbook of modern item response theory*
#' (pp. 433-448). New York: Springer.
#'
#' Kim, S. (2006). A comparative study of IRT fixed parameter calibration
#' methods. *Journal of Educational Measurement, 43*(4), 355-381.
#'
#' Kim, S., & Kolen, M. J. (2016). Multiple group IRT fixed-parameter estimation
#' for maintaining an established ability scale. *Center for Advanced
#' Studies in Measurement and Assessment Report, 49.*
#'
#' Meilijson, I. (1989). A fast improvement to the EM algorithm on its own
#' terms. *Journal of the Royal Statistical Society: Series B
#' (Methodological), 51*, 127-138.
#'
#' Woods, C. M. (2007). Empirical histograms in item response theory with
#' ordinal data. *Educational and Psychological Measurement, 67*(1), 73-87.
#'
#' @examples
#' \donttest{
#' ## ------------------------------------------------------------------------------
#' # 1. MG calibration using the simMG data
#' #  - Details:
#' #    (a) Constrain common items between groups to have
#' #        identical item parameters (i.e., items C1I1–C1I12 between
#' #        Groups 1 and 2, and items C2I1–C2I10 between Groups 2 and 3).
#' #    (b) Freely estimate the means and variances of the ability
#' #        distributions for all groups except the reference group,
#' #        where the mean and variance are fixed to 0 and 1, respectively.
#' ## ------------------------------------------------------------------------------
#' # 1-(1). Freely estimate the means and variances of Groups 2 and 3
#' # Import the true item metadata for the three groups
#' x <- simMG$item.prm
#'
#' # Extract model, score category, and item ID information
#' # from the item metadata for the three groups
#' model <- list(x$Group1$model, x$Group2$model, x$Group3$model)
#' cats <- list(x$Group1$cats, x$Group2$cats, x$Group3$cats)
#' item.id <- list(x$Group1$id, x$Group2$id, x$Group3$id)
#'
#' # Import the simulated response data sets for the three groups
#' data <- simMG$res.dat
#'
#' # Import the group names for the three groups
#' group.name <- simMG$group.name
#'
#' # Specify Groups 2 and 3 as the free groups where the scale
#' # of the ability distributions will be freely estimated.
#' # Group 1 will serve as the reference group, where the scale
#' # of the ability distribution is fixed to the values specified
#' # via the 'group.mean' and 'group.var' arguments
#' free.group <- c(2, 3) # or use 'free.group <- group.name[2:3]'
#'
#' # Estimate IRT parameters:
#' # As long as common items across groups share the same item IDs,
#' # their item parameters will be constrained to be equal across groups
#' # unless FIPC is implemented.
#' fit.1 <-
#'   est_mg(
#'     data = data, group.name = group.name, model = model,
#'     cats = cats, item.id = item.id, D = 1, free.group = free.group,
#'     use.gprior = TRUE, gprior = list(dist = "beta", params = c(5, 16)),
#'     group.mean = 0, group.var = 1, EmpHist = TRUE, Etol = 0.001, MaxE = 500
#'   )
#'
#' # Summary of the estimation
#' summary(fit.1)
#'
#' # Extract the item parameter estimates
#' getirt(fit.1, what = "par.est")
#'
#' # Extract the standard error estimates
#' getirt(fit.1, what = "se.est")
#'
#' # Extract the group-level parameter estimates (i.e., scale parameters)
#' getirt(fit.1, what = "group.par")
#'
#' # Extract the posterior latent ability distributions for each group
#' getirt(fit.1, what = "weights")
#'
#' # 1-(2). Alternatively, the same parameter estimation can be performed by
#' # inserting a list of item metadata for the groups into the 'x' argument.
#' # If the item metadata contains item parameters to be used as starting values,
#' # set 'use.startval = TRUE'.
#' # Also, specify the groups in which the ability distribution scales
#' # will be freely estimated using their group names.
#' free.group <- group.name[2:3]
#' fit.2 <-
#'   est_mg(
#'     x = x, data = data, group.name = group.name, D = 1,
#'     free.group = free.group, use.gprior = TRUE,
#'     gprior = list(dist = "beta", params = c(5, 16)),
#'     group.mean = 0, group.var = 1, EmpHist = TRUE, use.startval = TRUE,
#'     Etol = 0.001, MaxE = 500
#'   )
#'
#' # Summary of the estimation
#' summary(fit.2)
#'
#' ## ------------------------------------------------------------------------------
#' # 2. MG calibration with FIPC using simMG data
#' #  - Details:
#' #    (a) Fix the parameters of the common items between the groups
#' #        (i.e., items C1I1–C1I12 between Groups 1 and 2, and
#' #        items C2I1–C2I10 between Groups 2 and 3)
#' #    (b) Freely estimate the means and variances of the ability
#' #        distributions for all three groups
#' ## ------------------------------------------------------------------------------
#' # 2-(1). Freely estimate the means and variances for all three groups
#' # Set all three groups as free groups in which the scales
#' # of the ability distributions are to be freely estimated
#' free.group <- 1:3 # or use 'free.group <- group.name'
#'
#' # Specify the locations of items to be fixed in each group's metadata
#' # For Group 1: C1I1–C1I12 are located in rows 1–10 and 49–50
#' # For Group 2: C1I1–C1I12 are in rows 1–12, and
#' #              C2I1–C2I10 are in rows 41–50
#' # For Group 3: C2I1–C2I10 are in rows 1–10
#' fix.loc <- list(
#'   c(1:10, 49:50),
#'   c(1:12, 41:50),
#'   c(1:10)
#' )
#'
#' # Estimate IRT parameters using MG-FIPC:
#' # When FIPC is implemented, item metadata for all groups
#' # must be provided via the 'x' argument.
#' # For fixed items, their item parameters must be specified
#' # in the metadata. For non-fixed items, any placeholder values
#' # can be used in the metadata.
#' # Also set fipc = TRUE and fipc.method = "MEM"
#' fit.3 <-
#'   est_mg(
#'     x = x, data = data, group.name = group.name, D = 1,
#'     free.group = free.group, use.gprior = TRUE,
#'     gprior = list(dist = "beta", params = c(5, 16)),
#'     EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE,
#'     fipc.method = "MEM", fix.loc = fix.loc
#'   )
#'
#' # Summary of the estimation
#' summary(fit.3)
#'
#' # Extract the item parameter estimates
#' getirt(fit.3, what = "par.est")
#'
#' # Extract the standard error estimates
#' getirt(fit.3, what = "se.est")
#'
#' # Extract the group parameter estimates (i.e., scale parameters)
#' getirt(fit.3, what = "group.par")
#'
#' # Extract the posterior latent ability distributions of the groups
#' getirt(fit.3, what = "weights")
#'
#' # 2-(2). Alternatively, MG-FIPC can be implemented by specifying the
#' # IDs of the items to be fixed using the 'fix.id' argument.
#' # Provide a character vector of fixed item IDs to 'fix.id'
#' fix.id <- c(paste0("C1I", 1:12), paste0("C2I", 1:10))
#' fit.4 <-
#'   est_mg(
#'     x = x, data = data, group.name = group.name, D = 1,
#'     free.group = free.group, use.gprior = TRUE,
#'     gprior = list(dist = "beta", params = c(5, 16)),
#'     EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE,
#'     fipc.method = "MEM", fix.id = fix.id
#'   )
#'
#' # Summary of the estimation
#' summary(fit.4)
#'
#' ## ------------------------------------------------------------------------------
#' # 3. MG calibration with FIPC using simMG data
#' #    (Estimate group parameters only)
#' #  - Details:
#' #    (a) Fix all item parameters across all three groups
#' #    (b) Freely estimate the means and variances of the ability
#' #        distributions for all three groups
#' ## ------------------------------------------------------------------------------
#' # 3-(1). Freely estimate the means and variances for all three groups
#' # Set all three groups as free groups in which the scales
#' # of the ability distributions will be freely estimated
#' free.group <- 1:3 # or use 'free.group <- group.name'
#'
#' # Specify the locations of all fixed items in each group's metadata
#' fix.loc <- list(1:50, 1:50, 1:38)
#'
#' # Estimate group parameters only using MG-FIPC
#' fit.5 <-
#'   est_mg(
#'     x = x, data = data, group.name = group.name, D = 1,
#'     free.group = free.group, use.gprior = TRUE,
#'     gprior = list(dist = "beta", params = c(5, 16)),
#'     EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE,
#'     fipc.method = "MEM", fix.loc = fix.loc
#'   )
#'
#' # Summary of the estimation
#' summary(fit.5)
#'
#' # Extract the group parameter estimates (i.e., scale parameters)
#' getirt(fit.5, what = "group.par")
#'
#' ## ------------------------------------------------------------------------------
#' # 4. MG calibration with FIPC using simMG data
#' #    (Fix only the unique items of Group 1)
#' #  - Details:
#' #    (a) Fix item parameters of the unique items in Group 1 only
#' #    (b) Constrain the common items across groups to have
#' #        the same item parameters (i.e., C1I1–C1I12 between
#' #        Groups 1 and 2, and C2I1–C2I10 between Groups 2 and 3)
#' #    (c) Freely estimate the means and variances of the ability
#' #        distributions for all three groups
#' ## ------------------------------------------------------------------------------
#' # 4-(1). Freely estimate the means and variances for all three groups
#' # Set all three groups as free groups in which the scales
#' # of the ability distributions will be freely estimated
#' free.group <- group.name # or use 'free.group <- 1:3'
#'
#' # Specify the item IDs of the unique items in Group 1 to be fixed using
#' # the `fix.id` argument.
#' fix.id <- paste0("G1I", 1:38)
#'
#' # Alternatively, use the 'fix.loc' argument as
#' # 'fix.loc = list(11:48, NULL, NULL)'
#'
#' # Estimate IRT parameters using MG-FIPC
#' fit.6 <-
#'   est_mg(
#'     x = x, data = data, group.name = group.name, D = 1,
#'     free.group = free.group, use.gprior = TRUE,
#'     gprior = list(dist = "beta", params = c(5, 16)),
#'     EmpHist = TRUE, Etol = 0.001, MaxE = 500, fipc = TRUE,
#'     fipc.method = "MEM", fix.loc = NULL, fix.id = fix.id
#'   )
#'
#' # Summary of the estimation
#' summary(fit.6)
#'
#' # Extract the group parameter estimates (i.e., scale parameters)
#' getirt(fit.6, what = "group.par")
#'
#' }
#'
#'
#' @export
est_mg <- function(x = NULL,
                   data,
                   group.name = NULL,
                   D = 1,
                   model = NULL,
                   cats = NULL,
                   item.id = NULL,
                   free.group = NULL,
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
                   Etol = 1e-03,
                   MaxE = 500,
                   control = list(eval.max = 200, iter.max = 200),
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
    est_par <- est_mg_em(
      x = x, data = data, group.name = group.name, D = D, model = model, cats = cats, item.id = item.id,
      free.group = free.group, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g, a.val.1pl = a.val.1pl,
      a.val.gpcm = a.val.gpcm, g.val = g.val, use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
      aprior = aprior, bprior = bprior, gprior = gprior, missing = missing, Quadrature = Quadrature, weights = weights,
      group.mean = group.mean, group.var = group.var, EmpHist = EmpHist, use.startval = use.startval, Etol = Etol, MaxE = MaxE,
      control = control, se = se, verbose = verbose
    )
  } else {
    # implement FIPC method
    est_par <- est_mg_fipc(
      x = x, data = data, group.name = group.name, D = D, model = model, cats = cats, item.id = item.id,
      free.group = free.group, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g, a.val.1pl = a.val.1pl,
      a.val.gpcm = a.val.gpcm, g.val = g.val, use.aprior = use.aprior, use.bprior = use.bprior, use.gprior = use.gprior,
      aprior = aprior, bprior = bprior, gprior = gprior, missing = missing, Quadrature = Quadrature, weights = weights,
      group.mean = group.mean, group.var = group.var, EmpHist = EmpHist, use.startval = use.startval, Etol = Etol, MaxE = MaxE,
      control = control, fipc = fipc, fipc.method = fipc.method, fix.loc = fix.loc, fix.id = fix.id, se = se, verbose = verbose
    )
  }

  # return the estimation results
  class(est_par) <- "est_mg"
  est_par$call <- cl
  est_par
}


# multiple group calibration via EM
#' @import dplyr
est_mg_em <- function(x = NULL,
                      data,
                      group.name = NULL,
                      D = 1,
                      model = NULL,
                      cats = NULL,
                      item.id = NULL,
                      free.group = NULL,
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
                      Etol = 1e-03,
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

  # start parsing inputs
  if (verbose) {
    cat("Parsing input...", "\n")
  }

  # count the number of groups
  ngroup <- length(data)
  if (ngroup == 1) {
    stop("Use the est_irt() function when a single-group item calibration needs to be implemented.", call. = FALSE)
  }

  # create a group name vector when group.name = NULL
  if (is.null(group.name)) {
    group.name <- paste0("g", 1:ngroup)
  }

  # test length of each group
  nitem.gr <- purrr::map_dbl(.x = data, ncol)
  names(nitem.gr) <- group.name

  # sample size of each group
  nstd.gr <- purrr::map_dbl(.x = data, nrow)
  names(nstd.gr) <- group.name

  # extract information about the number of score categories and models
  if (!is.null(x)) {
    # confirm and correct all item metadata information
    x.gr <- purrr::map(.x = x, ~ {
      confirm_df(.x)
    })
    names(x.gr) <- group.name

    # item, cats, model information for each group
    item.id <- purrr::map(.x = x.gr, ~ {
      .x$id
    })
    cats.gr <- purrr::map(.x = x.gr, ~ {
      .x$cats
    })
    model.gr <- purrr::map(.x = x.gr, ~ {
      .x$model
    })
    if (!use.startval) {
      x.gr <-
        purrr::map(
          .x = x.gr,
          ~ {
            startval_df(cats = .x$cats, model = .x$model, item.id = .x$id)
          }
        )
    }
  } else {
    # check if the item IDs are specified when x = NULL
    if (is.null(item.id)) {
      stop("A list containing the vectors of the item IDs across all groups must be specified in the argument 'item.id'.",
        call. = FALSE
      )
    }

    # assign group name to the item id list
    names(item.id) <- group.name

    # IRT models for each group test form
    model.gr <-
      purrr::map2(
        .x = model, .y = nitem.gr,
        .f = function(x, y) {
          if (length(x) == 1) {
            rep(x, y)
          } else {
            x
          }
        }
      ) %>%
      purrr::map(~ {
        toupper(.x)
      })
    names(model.gr) <- group.name

    # score categories for each group test form
    cats.gr <-
      purrr::map2(
        .x = cats, .y = nitem.gr,
        .f = function(x, y) {
          if (length(x) == 1) {
            rep(x, y)
          } else {
            x
          }
        }
      )

    # create score category vectors when cats = NULL & all models are dichotomous item models
    if (is.null(cats)) {
      cats.gr <-
        purrr::map2(
          .x = nitem.gr, .y = model.gr,
          .f = function(x, y) {
            if (all(y %in% c("1PLM", "2PLM", "3PLM", "DRM"))) {
              rep(2, x)
            } else {
              stop("The number of score categories for the items should be specified in the argument 'cats'.", call. = FALSE)
            }
          }
        )
    }

    # create score category vectors when any of internal objects is Null in the cats argument &
    # all models are dichotomous item models
    cats.gr <-
      purrr::pmap(
        .l = list(x = cats.gr, y = nitem.gr, z = model.gr),
        .f = function(x, y, z) {
          if (is.null(x)) {
            if (all(z %in% c("1PLM", "2PLM", "3PLM", "DRM"))) {
              rep(2, y)
            } else {
              stop("The number of score categories for the items should be specified in the argument 'cats'.", call. = FALSE)
            }
          } else {
            x
          }
        }
      )
    names(cats.gr) <- group.name

    # create the item metadata for each group test form containing starting values
    x.gr <-
      purrr::pmap(
        .l = list(x = cats.gr, y = model.gr, z = item.id),
        ~ {
          startval_df(cats = ..1, model = ..2, item.id = ..3)
        }
      )
    names(x.gr) <- group.name
  }

  # copy data to data.gr object
  data.gr <- data
  names(data.gr) <- group.name

  # create combined response data sets across all groups
  # and assign item IDs
  data <-
    purrr::map(.x = data.gr, data.frame) %>%
    purrr::map2(
      .y = item.id,
      .f = function(x, y) {
        dplyr::rename_all(.tbl = x, .funs = function(x) y)
      }
    ) %>%
    dplyr::bind_rows()

  # a vector of item ID for the combined data set
  id <- colnames(data)

  # create an item metadata for the combined data est
  x <-
    dplyr::bind_rows(x.gr) %>%
    dplyr::distinct(.data$id, .keep_all = TRUE)
  if (!all(id == x$id)) {
    x <- x[match(id, x$id), ]
  }
  rownames(x) <- 1:nrow(x)

  # create two vectors of model & score categories for the combined data set
  model <- x$model
  cats <- x$cats

  # recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # transform a data set to matrix
  data <- data.matrix(data)

  # check the total number of examinees for the combined data
  nstd <- nrow(data)

  # check the total number of items in the combined data
  nitem <- ncol(data)

  # designate a reference group whose group parameters will be fixed
  if (is.null(free.group)) {
    ref.group <- 1:ngroup
  } else {
    if (is.character(free.group)) {
      free.group <- which(group.name %in% free.group)
      if (length(free.group) == 0) {
        stop("The group(s) in which abiltiy distribution(s) is(are) freely estimated do(does) not exist in the group.name argument.", call. = FALSE)
      }
      ref.group <- c(1:ngroup)[-free.group]
    } else {
      if (sum(1:ngroup %in% free.group) == 0) {
        stop("The group(s) in which abiltiy distribution(s) is(are) freely estimated do(does) not exist in the group name list.", call. = FALSE)
      }
      ref.group <- c(1:ngroup)[-free.group]
    }
  }

  # find the groups that have no common items with other groups
  pair <-
    t(utils::combn(group.name, m = 2)) %>%
    data.frame()
  overlap <- c()
  for (i in 1:nrow(pair)) {
    freq.tmp <-
      unlist(item.id[unlist(pair[i, ])]) %>%
      table()
    overlap[i] <- any(freq.tmp > 1)
  }
  link.gr <-
    dplyr::mutate(pair, overlap = overlap) %>%
    dplyr::filter(.data$overlap == TRUE) %>%
    dplyr::select(1:2) %>%
    unlist() %>%
    unique()
  nolink.gr <- c(1:ngroup)[!group.name %in% link.gr]

  # check if there exist free groups which do not share common items with other groups
  nolink.free <- free.group[free.group %in% nolink.gr]

  # when there exist free groups that do not share any common items,
  # stop the further process
  if (length(nolink.free) > 0) {
    warning.memo <-
      paste0(
        paste(group.name[nolink.free], collapse = ", "),
        " group(s) whose ability distribution(s) is(are) freely estimated do(does) share common items with other groups. \n",
        "Please specify the freely estimated group(s) correctly."
      )
    stop(paste0(warning.memo, " \n"), call. = FALSE)
  }

  # check the number of item responses across all items
  n.resp <- Rfast::colsums(!is.na(data))

  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)
  if (length(loc_allmiss) > 0L) {
    allmiss_item <- id[loc_allmiss]
    memo2 <- paste0(paste0(allmiss_item, collapse = ", "), " has/have no item response data across all groups. \n")
    stop(memo2, call. = FALSE)
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
    # create quadrature points
    quadpt <- seq(-Quadrature[2], Quadrature[2], length.out = Quadrature[1])

    # create the data.frame containing the quadrature points and weights
    weights <- gen.weight(dist = "norm", mu = group.mean, sigma = sqrt(group.var), theta = quadpt)
    n.quad <- length(quadpt)
  } else {
    quadpt <- weights[, 1]
    n.quad <- length(quadpt)
    moments.tmp <- cal_moment(node = quadpt, weight = weights[, 2])
    group.mean <- moments.tmp[1]
    group.var <- moments.tmp[2]
  }

  # a list containing the weights and densities for each group
  weights.gr <- replicate(n = ngroup, expr = weights, simplify = FALSE)
  names(weights.gr) <- group.name

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

  # create indices of samples (examinees) who belong to each group
  idx.tmp2 <- as.numeric(cumsum(nstd.gr))
  idx.tmp1 <- dplyr::lag(idx.tmp2, default = 0) + 1
  idx.std <- purrr::map2(.x = idx.tmp1, .y = idx.tmp2, ~ {
    .x:.y
  })
  names(idx.std) <- group.name

  # estimation
  if (verbose) {
    cat("Estimating item parameters...", "\n")
  }

  # implement EM algorithm
  time1 <- Sys.time()
  for (r in 1:MaxE) {
    # implement E-step
    estep <- Estep(
      elm_item = elm_item, idx.drm = idx.drm, idx.prm = idx.prm,
      data_drm = data_drm, data_prm = data_prm, data_all = data_all,
      weights = weights.gr, D = D, idx.std = idx.std
    )

    # implement M-step
    mstep <- Mstep(
      estep = estep, id = id, cats = cats, model = model, quadpt = quadpt, n.quad = n.quad,
      D = D, cols.item = cols.item, loc_1p_const = loc_1p_const, loc_else = loc_else, idx4est = idx4est,
      n.1PLM = n.1PLM, EmpHist = EmpHist, weights = weights.gr, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g,
      a.val.1pl = a.val.1pl, a.val.gpcm = a.val.gpcm, g.val = g.val, use.aprior = use.aprior, use.bprior = use.bprior,
      use.gprior = use.gprior, aprior = aprior, bprior = bprior, gprior = gprior, group.mean = group.mean,
      group.var = group.var, nstd = nstd.gr, Quadrature = Quadrature, control = control,
      iter = r, fipc = FALSE, reloc.par = param_loc$reloc.par, ref.group = ref.group,
      free.group = free.group, parbd = parbd
    )

    # compute the difference between previous and updated item parameter estimates
    diff_par <- mstep$elm_item$pars - elm_item$pars
    max.diff <- abs(max(diff_par, na.rm = TRUE))

    # loglikelihood value
    llike <- do.call(what = "sum", args = mstep$loglike)

    # print
    if (verbose) {
      cat("\r", paste0(
        "EM iteration: ", r, ", Loglike: ",
        format(round(llike, 4), nsmall = 4), ", Max-Change: ",
        format(round(max.diff, 6), nsmall = 5)
      ))
    }

    # check the convergence of EM algorithm
    converge <- max.diff <= Etol

    # extract the updated item parameter estimates
    elm_item$pars <- mstep$elm_item$pars

    # the quadrature points and the corresponding weights of the prior population density
    weights.gr <- mstep$weights

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
    weights = weights.gr, D = D, idx.std = idx.std
  )

  # compute the final log of marginal likelihood
  likehd.gr <- purrr::map(.x = idx.std, ~ {
    estep$likehd[.x, ]
  })
  llike.gr <- purrr::map2(
    .x = likehd.gr, .y = weights.gr,
    .f = ~ {
      sum(log(.x %*% matrix(.y[, 2])))
    }
  )
  llike <- do.call(what = "sum", args = llike.gr)

  # compute the mean and variance of the estimated density distributions
  pop_moments <-
    purrr::map(.x = weights.gr, ~ {
      cal_moment(node = quadpt, weight = .x[, 2])
    }) %>%
    purrr::map(~ {
      c(.x, sigma = sqrt(as.numeric(.x[2])))
    })
  if (length(ref.group) > 0) {
    pop_moments[ref.group] <-
      purrr::map(
        .x = 1:length(ref.group),
        ~ {
          c(mu = group.mean, sigma2 = group.var, sigma = sqrt(group.var))
        }
      )
  }

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

  # divide the parameter estimation results to each group
  full_all_df_gr <-
    purrr::map(
      .x = item.id,
      .f = function(x) {
        df.tmp <-
          full_all_df %>%
          dplyr::filter(.data$id %in% x)
        df.tmp <- df.tmp[match(x, df.tmp$id), ]
        rownames(df.tmp) <- 1:nrow(df.tmp)
        df.tmp
      }
    )
  names(full_all_df_gr) <- group.name
  par_df_gr <-
    purrr::map(
      .x = item.id,
      .f = function(x) {
        df.tmp <-
          par_df %>%
          dplyr::filter(.data$id %in% x)
        df.tmp <- df.tmp[match(x, df.tmp$id), ]
        rownames(df.tmp) <- 1:nrow(df.tmp)
        df.tmp
      }
    )
  names(par_df_gr) <- group.name
  se_df_gr <-
    purrr::map(
      .x = item.id,
      .f = function(x) {
        df.tmp <-
          se_df %>%
          dplyr::filter(.data$id %in% x)
        df.tmp <- df.tmp[match(x, df.tmp$id), ]
        rownames(df.tmp) <- 1:nrow(df.tmp)
        df.tmp
      }
    )
  names(se_df_gr) <- group.name

  # population density parameters
  group.par <- vector("list", ngroup)
  for (i in 1:ngroup) {
    if (i %in% ref.group) {
      moments.se <- rep(NA_real_, 3)
      group.par.tmp <- data.frame(rbind(pop_moments[[i]], moments.se))
    } else {
      moments.est <- pop_moments[[i]]
      se.mu <- moments.est[3] / sqrt(nstd.gr[i])
      se.sigma2 <- moments.est[2] * sqrt(2 / (nstd.gr[i] - 1))
      se.sigma <- (1 / (2 * moments.est[3])) * se.sigma2 # using Delta method
      moments.se <- c(se.mu, se.sigma2, se.sigma)
      group.par.tmp <- data.frame(rbind(moments.est, moments.se))
    }
    colnames(group.par.tmp) <- c("mu", "sigma2", "sigma")
    rownames(group.par.tmp) <- c("estimates", "se")
    group.par[[i]] <- group.par.tmp
  }
  names(group.par) <- group.name

  # prior information
  if (use.aprior) aprior.dist <- aprior else aprior.dist <- NULL
  if (use.bprior) bprior.dist <- bprior else bprior.dist <- NULL
  if (use.gprior) gprior.dist <- gprior else gprior.dist <- NULL

  # statistics based on the loglikelihood of the fitted model:
  npar.est <- length(param_loc$reloc.par) + length(free.group) * 2
  neg2llke <- -2 * llike
  aic <- 2 * npar.est + neg2llke
  bic <- npar.est * log(nstd) + neg2llke

  ## ---------------------------------------------------------------
  # check end time
  end.time <- Sys.time()

  # record total computation time
  est_time3 <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 2)

  rst <- list(
    estimates = list(overall = full_all_df, group = full_all_df_gr),
    par.est = list(overall = par_df, group = par_df_gr),
    se.est = list(overall = se_df, group = se_df_gr),
    pos.par = loc_df, covariance = cov_mat, loglikelihood = list(overall = llike, group = llike.gr), aic = aic, bic = bic,
    group.par = group.par, weights = weights.gr, posterior.dist = post_dist, data = list(overall = data, group = data.gr),
    scale.D = D, ncase = list(overall = nstd, group = nstd.gr), nitem = list(overall = nitem, group = nitem.gr),
    Etol = Etol, MaxE = MaxE, aprior = aprior.dist, bprior = bprior.dist, gprior = gprior.dist, npar.est = npar.est,
    niter = r, maxpar.diff = max.diff, EMtime = est_time1, SEtime = est_time2, TotalTime = est_time3, test.1 = memo3,
    test.2 = memo4, var.note = memo5, fipc = FALSE, fipc.method = NULL, fix.loc = list(overall = NULL, group = NULL)
  )

  if (verbose) {
    cat("Estimation is finished in", est_time3, "seconds.", "\n")
  }
  return(rst)
}


# multiple group FIPC via EM
#' @import dplyr
est_mg_fipc <- function(x = NULL,
                        data,
                        group.name = NULL,
                        D = 1,
                        model = NULL,
                        cats = NULL,
                        item.id = NULL,
                        free.group = NULL,
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

  # start parsing inputs
  if (verbose) {
    cat("Parsing input...", "\n")
  }

  # count the number of groups
  ngroup <- length(data)
  if (ngroup == 1) {
    stop("Use the est_irt() function when a single-group item calibration is implemented.", call. = FALSE)
  }

  # create a group name vector when group.name = NULL
  if (is.null(group.name)) {
    group.name <- paste0("g", 1:ngroup)
  }

  # test length of each group
  nitem.gr <- purrr::map_dbl(.x = data, ncol)
  names(nitem.gr) <- group.name

  # sample size of each group
  nstd.gr <- purrr::map_dbl(.x = data, nrow)
  names(nstd.gr) <- group.name

  # confirm and correct all item metadata information
  x.gr <- purrr::map(.x = x, ~ {
    confirm_df(.x)
  })
  names(x.gr) <- group.name

  # item, cats, model information for each group
  if (is.null(item.id)) {
    item.id <- purrr::map(.x = x.gr, ~ {
      .x$id
    })
  } else {
    item.id <- item.id
  }
  cats.gr <- purrr::map(.x = x.gr, ~ {
    .x$cats
  })
  model.gr <- purrr::map(.x = x.gr, ~ {
    .x$model
  })

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
    fix.item <- unique(fix.id)
    fix.item.gr <- purrr::map(.x = item.id, ~ {
      .x[.x %in% fix.item]
    })
    fix.loc.gr <- purrr::map(.x = item.id, ~ {
      which(.x %in% fix.item)
    })
    names(fix.loc.gr) <- group.name
  } else {
    fix.item.gr <- purrr::map2(.x = item.id, .y = fix.loc, ~ {
      .x[.y]
    })
    fix.item <- unique(unlist(fix.item.gr))
    fix.loc.gr <- fix.loc
    names(fix.loc.gr) <- group.name
  }

  # copy data to data.gr object
  data.gr <- data
  names(data.gr) <- group.name

  # create combined response data sets across all groups
  # and assign item IDs
  data <-
    purrr::map(.x = data.gr, data.frame) %>%
    purrr::map2(
      .y = item.id,
      .f = function(x, y) {
        dplyr::rename_all(.tbl = x, .funs = function(x) y)
      }
    ) %>%
    dplyr::bind_rows()

  # a vector of item ID for the combined data set
  id <- colnames(data)

  # find the locations of items that should be fixed in the combined data set
  fix.loc <- which(id %in% fix.item)

  # create an item metadata for the combined data est
  x <-
    dplyr::bind_rows(x.gr) %>%
    dplyr::distinct(.data$id, .keep_all = TRUE)
  if (!all(id == x$id)) {
    x <- x[match(id, x$id), ]
  }
  rownames(x) <- 1:nrow(x)

  # recode missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # transform a data set to matrix
  data <- data.matrix(data)

  # check the total number of examinees for the combined data
  nstd <- nrow(data)

  # check the total number of items in the combined data
  nitem.all <- ncol(data)

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

  # clean the two item metadata sets
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

  # designate a reference group whose group parameters will be fixed
  if (is.null(free.group)) {
    stop("The free.group argument should not be NULL when FIPC = TRUE.", call. = FALSE)
  }
  if (is.character(free.group)) {
    free.group <- which(group.name %in% free.group)
    if (length(free.group) == 0) {
      stop("The group(s) in which abiltiy distribution(s) is(are) freely estimated do(does) not exist in the group.name argument.", call. = FALSE)
    }
    ref.group <- c(1:ngroup)[-free.group]
  } else {
    if (sum(1:ngroup %in% free.group) == 0) {
      stop("The group(s) in which abiltiy distribution(s) is(are) freely estimated do(does) not exist in the group name list.", call. = FALSE)
    }
    ref.group <- c(1:ngroup)[-free.group]
  }

  # find the groups that have no common items with other groups
  pair <-
    t(utils::combn(group.name, m = 2)) %>%
    data.frame()
  overlap <- c()
  for (i in 1:nrow(pair)) {
    freq.tmp <-
      unlist(item.id[unlist(pair[i, ])]) %>%
      table()
    overlap[i] <- any(freq.tmp > 1)
  }
  link.gr <-
    dplyr::mutate(pair, overlap = overlap) %>%
    dplyr::filter(.data$overlap == TRUE) %>%
    dplyr::select(1:2) %>%
    unlist() %>%
    unique()
  nolink.gr <- c(1:ngroup)[!group.name %in% link.gr]

  # check if there exist free groups which do not share common items with other groups
  nolink.free <- free.group[free.group %in% nolink.gr]

  # count the number of fixed items in each group
  n.fix.gr <- purrr::map_dbl(.x = fix.item.gr, ~ {
    length(.x)
  })

  # when there exist groups that are specified as free groups but not have fixed items,
  # stop the further process
  if (!all(nolink.free %in% which(n.fix.gr != 0))) {
    warning.memo <-
      paste0(
        paste(group.name[nolink.free[nolink.free %in% which(n.fix.gr != 0)]], collapse = ", "),
        " group(s) in which ability distribution(s) is(are) freely estimated do(does) not share common items with other groups \n",
        "nor include fixed items. Please specify the freely estimated group(s) correctly."
      )
    stop(paste0(warning.memo, " \n"), call. = FALSE)
  }

  # check the number of item responses across all items
  n.resp <- Rfast::colsums(!is.na(data))

  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)
  if (length(loc_allmiss) > 0L) {
    allmiss_item <- id[loc_allmiss]
    memo2 <- paste0(paste0(allmiss_item, collapse = ", "), " has/have no item response data across all groups. \n")
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
  # create initial weights of prior ability distribution when it is not specified
  if (is.null(weights)) {
    # create quadrature points
    quadpt <- seq(-Quadrature[2], Quadrature[2], length.out = Quadrature[1])

    # create the data.frame containing the quadrature points and weights
    weights <- gen.weight(dist = "norm", mu = group.mean, sigma = sqrt(group.var), theta = quadpt)
    n.quad <- length(quadpt)
  } else {
    quadpt <- weights[, 1]
    n.quad <- length(quadpt)
    moments.tmp <- cal_moment(node = quadpt, weight = weights[, 2])
    group.mean <- moments.tmp[1]
    group.var <- moments.tmp[2]
  }

  # a list containing the weights and densities for each group
  weights.gr <- replicate(n = ngroup, expr = weights, simplify = FALSE)
  names(weights.gr) <- group.name

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
    mmt_dist_old <- matrix(rep(c(group.mean, group.var), ngroup), ncol = ngroup)
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

  # create indices of samples (examinees) who belong to each group
  idx.tmp2 <- as.numeric(cumsum(nstd.gr))
  idx.tmp1 <- dplyr::lag(idx.tmp2, default = 0) + 1
  idx.std <- purrr::map2(.x = idx.tmp1, .y = idx.tmp2, ~ {
    .x:.y
  })
  names(idx.std) <- group.name

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
          data_all1 = data.new, weights = weights.gr, D = D, idx.std = idx.std
        )
      } else {
        estep <- Estep_fipc(
          elm_item1 = elm_item_new, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
          idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
          data_all1 = data.new, weights = weights.gr, D = D, idx.std = idx.std
        )
      }
    } else {
      estep <- Estep_fipc(
        elm_item1 = elm_item_all, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
        idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
        data_all1 = data.all, weights = weights.gr, D = D, idx.std = idx.std
      )
      estep$elm_item <- NULL
    }

    # implement M-step
    mstep <- Mstep(
      estep = estep, id = id, cats = cats, model = model, quadpt = quadpt, n.quad = n.quad,
      D = D, cols.item = cols.item, loc_1p_const = loc_1p_const, loc_else = loc_else, idx4est = idx4est,
      n.1PLM = n.1PLM, EmpHist = EmpHist, weights = weights.gr, fix.a.1pl = fix.a.1pl, fix.a.gpcm = fix.a.gpcm, fix.g = fix.g,
      a.val.1pl = a.val.1pl, a.val.gpcm = a.val.gpcm, g.val = g.val, use.aprior = use.aprior, use.bprior = use.bprior,
      use.gprior = use.gprior, aprior = aprior, bprior = bprior, gprior = gprior, group.mean = group.mean,
      group.var = group.var, nstd = nstd.gr, Quadrature = Quadrature, control = control,
      iter = r, fipc = TRUE, reloc.par = param_loc$reloc.par, ref.group = ref.group,
      free.group = free.group, parbd = parbd
    )

    if (!is.null(x_new)) {
      # compute the difference between previous and updated item parameter estimates
      diff_par <- mstep$elm_item$pars - elm_item_new$pars
      max.diff <- abs(max(diff_par, na.rm = TRUE))
    } else {
      # compute the mean and sd of the updated prior distribution
      mmt_dist_new <-
        as.matrix(purrr::map_dfc(.x = mstep$weights, ~ {
          cal_moment(node = .x$theta, weight = .x$weight)
        }))
      diff_par <- mmt_dist_new - mmt_dist_old
      max.diff <- abs(max(diff_par, na.rm = TRUE))
    }

    # loglikelihood value
    llike <- do.call(what = "sum", args = mstep$loglike)

    # print
    if (verbose) {
      cat("\r", paste0(
        "EM iteration: ", r, ", Loglike: ",
        format(round(llike, 4), nsmall = 4), ", Max-Change: ",
        format(round(max.diff, 6), nsmall = 5)
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
    weights.gr <- mstep$weights

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

  # conduct one more E-step to update the posterior distribution using the final item parameter estimates
  if (!is.null(x_new)) {
    estep <- Estep_fipc(
      elm_item1 = elm_item_new, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
      idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
      data_all1 = data.new, weights = weights.gr, D = D, idx.std = idx.std
    )
  } else {
    estep <- Estep_fipc(
      elm_item1 = elm_item_all, elm_item2 = elm_item_all, idx.drm2 = idx.drm.all,
      idx.prm2 = idx.prm.all, data_drm2 = data_drm.all, data_prm2 = data_prm.all,
      data_all1 = data.all, weights = weights.gr, D = D, idx.std = idx.std
    )
  }

  # compute the final log of marginal likelihood
  likehd.gr <- purrr::map(.x = idx.std, ~ {
    estep$likehd[.x, ]
  })
  llike.gr <- purrr::map2(
    .x = likehd.gr, .y = weights.gr,
    .f = ~ {
      sum(log(.x %*% matrix(.y[, 2])))
    }
  )
  llike <- do.call(what = "sum", args = llike.gr)

  # compute the mean and variance of the estimated density distributions
  pop_moments <-
    purrr::map(.x = weights.gr, ~ {
      cal_moment(node = quadpt, weight = .x[, 2])
    }) %>%
    purrr::map(~ {
      c(.x, sigma = sqrt(as.numeric(.x[2])))
    })
  if (length(ref.group) > 0) {
    pop_moments[ref.group] <-
      purrr::map(
        .x = 1:length(ref.group),
        ~ {
          c(mu = group.mean, sigma2 = group.var, sigma = sqrt(group.var))
        }
      )
  }

  ## ---------------------------------------------------------------------
  # estimates the information matrix and standard errors of item parameter estimates
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

  # divide the parameter estimation results to each group
  full_all_df_gr <-
    purrr::map(
      .x = item.id,
      .f = function(x) {
        df.tmp <-
          full_all_df %>%
          dplyr::filter(.data$id %in% x)
        df.tmp <- df.tmp[match(x, df.tmp$id), ]
        rownames(df.tmp) <- 1:nrow(df.tmp)
        df.tmp
      }
    )
  names(full_all_df_gr) <- group.name
  par_df_gr <-
    purrr::map(
      .x = item.id,
      .f = function(x) {
        df.tmp <-
          x_all %>%
          dplyr::filter(.data$id %in% x)
        df.tmp <- df.tmp[match(x, df.tmp$id), ]
        rownames(df.tmp) <- 1:nrow(df.tmp)
        df.tmp
      }
    )
  names(par_df_gr) <- group.name
  se_df_gr <-
    purrr::map(
      .x = item.id,
      .f = function(x) {
        df.tmp <-
          se_all_df %>%
          dplyr::filter(.data$id %in% x)
        df.tmp <- df.tmp[match(x, df.tmp$id), ]
        rownames(df.tmp) <- 1:nrow(df.tmp)
        df.tmp
      }
    )
  names(se_df_gr) <- group.name

  # population density parameters
  group.par <- vector("list", ngroup)
  for (i in 1:ngroup) {
    if (i %in% ref.group) {
      moments.se <- rep(NA_real_, 3)
      group.par.tmp <- data.frame(rbind(pop_moments[[i]], moments.se))
    } else {
      moments.est <- pop_moments[[i]]
      se.mu <- moments.est[3] / sqrt(nstd.gr[i])
      se.sigma2 <- moments.est[2] * sqrt(2 / (nstd.gr[i] - 1))
      se.sigma <- (1 / (2 * moments.est[3])) * se.sigma2 # using Delta method
      moments.se <- c(se.mu, se.sigma2, se.sigma)
      group.par.tmp <- data.frame(rbind(moments.est, moments.se))
    }
    colnames(group.par.tmp) <- c("mu", "sigma2", "sigma")
    rownames(group.par.tmp) <- c("estimates", "se")
    group.par[[i]] <- group.par.tmp
  }
  names(group.par) <- group.name

  # prior information
  if (use.aprior) aprior.dist <- aprior else aprior.dist <- NULL
  if (use.bprior) bprior.dist <- bprior else bprior.dist <- NULL
  if (use.gprior) gprior.dist <- gprior else gprior.dist <- NULL

  # statistics based on the loglikelihood of the fitted model:
  if (!is.null(x_new)) {
    npar.est <- length(param_loc$reloc.par) + length(free.group) * 2
  } else {
    npar.est <- length(free.group) * 2
  }
  neg2llke <- -2 * llike
  aic <- 2 * npar.est + neg2llke
  bic <- npar.est * log(nstd) + neg2llke

  ## ---------------------------------------------------------------
  # check end time
  end.time <- Sys.time()

  # record total computation time
  est_time3 <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 2)

  rst <- list(
    estimates = list(overall = full_all_df, group = full_all_df_gr),
    par.est = list(overall = x_all, group = par_df_gr),
    se.est = list(overall = se_all_df, group = se_df_gr),
    pos.par = loc_all_df, covariance = cov_mat, loglikelihood = list(overall = llike, group = llike.gr), aic = aic, bic = bic,
    group.par = group.par, weights = weights.gr, posterior.dist = post_dist, data = list(overall = data_all, group = data.gr),
    scale.D = D, ncase = list(overall = nstd, group = nstd.gr), nitem = list(overall = nitem.all, group = nitem.gr),
    Etol = Etol, MaxE = MaxE, aprior = aprior.dist, bprior = bprior.dist, gprior = gprior.dist, npar.est = npar.est,
    niter = r, maxpar.diff = max.diff, EMtime = est_time1, SEtime = est_time2, TotalTime = est_time3, test.1 = memo3,
    test.2 = memo4, var.note = memo5, fipc = TRUE, fipc.method = fipc.method, fix.loc = list(overall = fix.loc, group = fix.loc.gr)
  )

  if (verbose) {
    cat("Estimation is finished in", est_time3, "seconds.", "\n")
  }
  return(rst)
}
