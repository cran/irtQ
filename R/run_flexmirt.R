#' Run flexMIRT from Within R
#'
#' This function runs flexMIRT (Cai, 2017) from within R by executing a model
#' specified in a flexMIRT syntax file (i.e., *.flexmirt). To use this function,
#' the flexMIRT software must be installed on your system. This interface is
#' especially useful for conducting simulation studies or automating batch
#' analyses involving flexMIRT.
#'
#' @param file.syntax A single string or character vector specifying the path(s)
#'   to one or more flexMIRT syntax files (with extension *.flexmirt) to be
#'   executed. For example: `"C:/Users/Data/irtmodel.flexmirt"`.
#' @param dir.flex A character string specifying the directory where flexMIRT is
#'   installed. The folder name typically includes "flexMIRT" (e.g.,
#'   "flexMIRT3", "flexMIRT 3.6"). If set to `NULL`, the function searches for
#'   flexMIRT in `"C:/Program Files"` and uses a default path if found (e.g.,
#'   `"C:/Program Files/flexMIRT3"`).
#' @param show.output.on.console Logical. If `TRUE`, the output of the system
#'   command is printed to the R console. Default is `FALSE`. See
#'   [base::system()].
#' @param ... Additional arguments passed to [base::system()].
#'
#' @details When using a version of flexMIRT earlier than 3.6, the directory
#'   specified in `dir.flex` must contain the following six files:
#' \itemize{
#'   \item \code{WinFlexMIRT.exe}
#'   \item \code{FlexMIRT_x64.exe}
#'   \item \code{FlexMIRT_x86.exe}
#'   \item \code{vpg.dll}
#'   \item \code{vpg.licensing.client.dll}
#'   \item \code{vpg.licensing.dll}
#' }
#'
#'   For flexMIRT version 3.6 or later, the directory must include the following
#'   five files:
#' \itemize{
#'   \item \code{WinFlexMIRT.exe}
#'   \item \code{vpg.dll}
#'   \item \code{vpg.licensing.client.dll}
#'   \item \code{vpg.licensing.dll}
#'   \item \code{VPGLicenseClientNet.dll}
#' }
#'   along with a subdirectory named \code{Resources} that contains the
#'   following two files:
#' \itemize{
#'   \item \code{flexMIRT_x64_AVX.exe}
#'   \item \code{flexMIRT_x86_AVX.exe}
#' }
#'
#' @return Output files generated by flexMIRT.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional
#'   item analysis and test scoring (Computer Software). Chapel Hill, NC: Vector
#'   Psychometric Group.
#'
#' @examples
#' # Examples below will run if the flexMIRT software is installed
#' # in the default directory "C:/Program Files/flexMIRT3".
#' # Otherwise, specify the directory where flexMIRT is installed
#' # using the 'dir.flex' argument.
#'
#' \dontrun{
#' # (1) Run a single syntax file
#' # Load an example flexMIRT syntax file for estimating item parameters using the 2PL model
#' file.syntax <- system.file("extdata", "2PLM_example.flexmirt", package = "irtQ")
#'
#' # Run flexMIRT to estimate item parameters for the 2PL model
#' run_flexmirt(file.syntax = file.syntax, dir.flex = NULL, show.output = TRUE)
#'
#' # Check the output file
#' out.file <- system.file("extdata", "2PLM_example-prm.txt", package = "irtQ")
#' bring.flexmirt(out.file, type = "par")
#'
#' # (2) Run multiple syntax files
#' # Load two example flexMIRT syntax files
#' file.syntax1 <- system.file("extdata", "2PLM_example.flexmirt", package = "irtQ")
#' file.syntax2 <- system.file("extdata", "3PLM_example.flexmirt", package = "irtQ")
#'
#' # Run flexMIRT to estimate item parameters for both models
#' run_flexmirt(file.syntax = c(file.syntax1, file.syntax2), dir.flex = NULL, show.output = FALSE)
#'
#' # Check the output files
#' out.file1 <- system.file("extdata", "2PLM_example-prm.txt", package = "irtQ")
#' out.file2 <- system.file("extdata", "3PLM_example-prm.txt", package = "irtQ")
#' bring.flexmirt(out.file1, type = "par")
#' bring.flexmirt(out.file2, type = "par")
#' }
#'
#' @export
run_flexmirt <- function(file.syntax,
                         dir.flex = NULL,
                         show.output.on.console = FALSE,
                         ...) {

  if (is.null(dir.flex)) {
    # find a default directory where flexMIRT is installed exists
    flex_exist <- dir(path = "C:/Program Files", pattern = "flexMIRT", full.names = TRUE)

    # provide a warning message when a default directory where flexMIRT is installed does not exist
    if (length(flex_exist) == 0) {
      stop("A default directory of flexMIRT in 'C:/Program Files/' does not exist.
            Please provide a full path of directory where flexMIRT is installed.", call. = TRUE)
    }

    # a path of flexMIRT execution file
    file.flex <- file.path(flex_exist, "WinFlexMIRT.exe")
  } else {
    # check if the provided directory where flexMIRT is installed exists
    flex_exist <- dir.exists(paths = dir.flex)

    # provide a warning message when the provided directory where flexMIRT is installed does not exist
    if (!flex_exist) {
      stop(paste0("The provided directory of '", dir.flex, "' does not exist.
           Please provide a full path of directory where flexMIRT is installed."), call. = TRUE)
    }

    # a path of flexMIRT execution file
    file.flex <- file.path(dir.flex, "WinFlexMIRT.exe")
  }

  # run flexMIRT
  # when a single syntax files are provided
  if (length(file.syntax) == 1) {
    system(
      command = paste("\"", file.flex, "\"", " -r ", "\"", file.syntax, "\"", sep = ""),
      show.output.on.console = show.output.on.console, ...
    )
  }

  # when multiple syntax files are provided
  if (length(file.syntax) > 1) {
    # start the progress bar
    pb <- utils::txtProgressBar(min = 0, max = 100, initial = 0, width = 50, style = 3)
    for (r in 1:length(file.syntax)) {
      system(
        command = paste("\"", file.flex, "\"", " -r ", "\"", file.syntax[r], "\"", sep = ""),
        show.output.on.console = show.output.on.console, ...
      )

      # modify the progress bar
      info <- sprintf("Parameter Estimation Progress: %d%%", round((r / length(file.syntax)) * 100))
      utils::setTxtProgressBar(pb, value = r / length(file.syntax) * 100, title = info, label = paste0("Estimation ", r, " is done"))
    }

    # closing the progress Bar
    close(pb)
  }
}
