#' Generate a correction matrix for cytometry data analysis
#'
#'
#' This function aids the correctUnmix function, to create a symmetrical
#' correction matrix that should be used together with a flowframe to correct
#' the errors of unmixing.
#' @param specMat The spectral matrix used to unmix the dataset of interest.
#' @return A symmetrical matrix of zeros with the right row- and column names.
#' @seealso \code{\link{correctUnmix}}
#' @importFrom BiocGenerics ncol colnames
#' @examples
#' # Load uncompensated data
#' data(fullPanel)
#'
#' # Load the spectral unmixing matrix generated with controls from the same
#' # experiment. These can be generated using the specMatCalc function.
#' data(specMat)
#'
#' # And now unmix
#' fullPanelUnmix <- specUnmix(fullPanel, specMat)
#'
#'
#' # Create an empty unmixinng matrix
#' corrMat <- corrMatCreate(specMat)
#'
#' @export corrMatCreate
corrMatCreate <- function(specMat) {
    corrNames <- rownames(specMat)
    correcMat <- matrix(0,
        nrow = length(corrNames),
        ncol = length(corrNames),
        dimnames =
            list(corrNames, corrNames)
    )
    diag(correcMat) <- 1
    return(correcMat)
}
