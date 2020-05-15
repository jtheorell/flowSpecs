#' Efficient inverse hyperbolic cosine transformation
#'
#' This is a simple wrapper function for the base asinh function, that is useful
#' for flowFrames and flowSets. It also allows for reversing the transformation
#' with the argument "unTrans".
#' @importFrom BiocGenerics colnames
#' @importFrom flowCore exprs fsApply
#' @param flowObj The fcs object to be transformed. Both flowFrames and flowSets
#' are accepted.
#' @param transNames The variables that should be normalized.
#' @param transCoFacs This value or vector of values define the values for the
#' transformation during the normalization. In the "default" case, the function
#' defines the object as a CyTOF object if >5 percent of the values are 0, and
#' applies the transformation value 5. Otherwise, the value 400 is applied.
#' @param unTrans If the reverse action should be taken, i.e. if an already
#' transformed dataset should be un-transformed. NB! It is of great importance
#' that the same transformation factors are used!
#' @return A flow object containing the transformed data, and with all metadata
#' left untouched.
#' @examples
#' # Import some data and the spectral matrix. The latter can be generated using
#' # specMatCalc
#' data(fullPanel)
#' data(specMat)
#'
#' fullPanelUnmixed <- specUnmix(fullPanel, specMat)
#'
#' # Identify the columns that should be transformed
#' colnames(fullPanelUnmixed)
#' # The time and scatter parameters should not, but apart from that, all should
#' # be included.
#' transNames <- colnames(fullPanelUnmixed)[seq(6,18)]
#'
#' # ow, transform this file, with the default transformation factor of 400.
#' # NB! It is alway advisable to check the data for the most optimal
#' # transformation factors. For flow cytometry data, this can potentially be
#' # done using the flowVS package:
#' # https://www.bioconductor.org/packages/release/bioc/html/flowVS.html
#'
#' fullPanelTrans <- arcTrans(fullPanelUnmixed, transNames)
#'
#' @export arcTrans
arcTrans <- function(flowObj, transNames, transCoFacs = "default",
                     unTrans = FALSE) {
    if (inherits(flowObj, "flowSet")) {
        focusFrame <- flowObj[[1]]
    } else if (inherits(flowObj, "flowFrame")) {
        focusFrame <- flowObj
    } else {
        stop("The flowObj needs to be either a flowSet or a flowFrame")
    }

    if (length(transCoFacs) == 1) {
        if (transCoFacs == "default") {
            transCoFacs <- massOrFlowTrans(
                focusFrame = focusFrame,
                transNames = transNames
            )
        } else if (is.numeric(transCoFacs)) {
            transCoFacs <- rep(transCoFacs, length(transNames))
            names(transCoFacs) <- transNames
        }
    } else {
        names(transCoFacs) <- transNames
    }

    if (inherits(flowObj, "flowSet")) {
        resultObj <- fsApply(flowObj, arcTransCoFunction,
            transCoFacs = transCoFacs,
            unTrans = unTrans
        )
    } else {
        resultObj <- arcTransCoFunction(
            focusFrame = flowObj,
            transCoFacs = transCoFacs,
            unTrans = unTrans
        )
    }
    return(resultObj)
}

arcTransCoFunction <- function(focusFrame, transCoFacs, unTrans) {
    focusFrameResult <- focusFrame
    if (unTrans) {
        for (i in names(transCoFacs)) {
            exprs(focusFrameResult)[, i] <-
                sinh(exprs(focusFrame)[, i]) * transCoFacs[i]
        }
    } else {
        for (i in names(transCoFacs)) {
            exprs(focusFrameResult)[, i] <-
                asinh(exprs(focusFrame)[, i] / transCoFacs[i])
        }
    }

    message("Transformation complete")

    return(focusFrameResult)
}
