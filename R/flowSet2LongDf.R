#' Convert a flowSet to one long dataframe with all identifiers as separate
#' #columns.
#'
#'
#' This function is mainly used for compatibility with matrix-based clustering
#' algorithms, such as depeche in the DepecheR package.
#'
#' @importFrom flowCore fsApply exprs flowSet sampleNames description
#' @param flowObj The flowSet or flowFrame to be converted to a dataframe.
#' @param idInfo A list of any number of characteristics that can be derived
#' from the file names. For each entry, a gsub specification of where to find
#' the information in the file name should be added, such as id=""..._|..."".
#' @return A long data frame with one column per PMT/APD (or fluorochrome,
#' depending on the state of the imported files), one for the acquisition date
#' (for fcs files) and one colum for each specified slot above. If no
#' gsub-pattern is provided, only a single column with the full file name will
#' be used to separate the observations from each file.
#' @seealso \code{\link[DepecheR]{depeche}}
#' @examples
#' #' # Load uncompensated data
#' data(fullPanel)
#'
#' # Load the spectral unmixing matrix generated with controls from the same
#' # experiment. These can be generated using the specMatCalc function.
#' data(specMat)
#'
#' # Now unmix
#' fullPanelUnmix <- specUnmix(fullPanel, specMat)
#'
#' # Transform all fluorescent channels
#' fullPanelTrans <- arcTrans(fullPanelUnmix,
#'     transNames = colnames(fullPanelUnmix)[6:18])
#'
#' # This function is primarily meant to be used with flowSets.
#' # If we had only one flowFrame, we could just extract the data by
#' # the use of the flowCore function exprs(), so we will convert the data to a
#' # flowSet now.
#' library(flowCore)
#' fullPanelFs <- flowSet(fullPanelTrans)
#'
#' # Before converting to a dataframe it is important to get an idea of the
#' # structure of the names, to be able to extract meaningful parts of the name.
#' # Here, we have an exceptional case again, as the flowSet has just been
#' # created, so there is actually no meaningful name of the flowFrame inside
#' # it. So for example reasons, we will give it one now:
#' sampleNames(fullPanelFs) <- "PBMC_full_panel_d1.fcs"
#'
#' # And now, we generate the dataframe:
#' fullPanelDf <- flowSet2LongDf(fullPanelFs, idInfo =
#'                                  list("Tissue" = "|_full_panel_..\\.fcs",
#'                                      "Donor" = "...._full_panel_|\\.fcs"))
#' # This is the result
#' str(fullPanelDf)
#'
#' @export flowSet2LongDf
flowSet2LongDf <- function(flowObj, idInfo) {
    if (inherits(flowObj, "flowFrame")) {
        flowObj <- flowSet(flowObj)
    }
    flowSetExprs <- data.frame(fsApply(flowObj, exprs))
    nameVector <- sampleNames(flowObj)

    if (missing(idInfo)) {
        flowSetExprs$names <-
            unlist(retrieveFlowSetNames(
                nameVector = nameVector,
                specFlowSet = flowObj, gsubpattern = ""
            ))
    } else {
        for (i in seq_along(idInfo)) {
            flowSetExprs$names <-
                unlist(retrieveFlowSetNames(
                    nameVector = nameVector,
                    specFlowSet = flowObj,
                    gsubpattern = idInfo[[i]]
                ))
            colnames(flowSetExprs)[which(colnames(flowSetExprs) == "names")] <-
                names(idInfo)[i]
        }
    }

    dateVector <- as.vector(fsApply(
        flowObj,
        function(x) description(x)$`$DATE`
    ))
    flowSetExprs$acqDate <-
        unlist(retrieveFlowSetNames(
            nameVector = dateVector,
            specFlowSet = flowObj,
            gsubpattern = ""
        ))

    return(flowSetExprs)
}

retrieveFlowSetNames <- function(nameVector, specFlowSet, gsubpattern) {
    lengthList <- as.list(fsApply(specFlowSet, function(x) nrow(exprs(x))))
    nameVectorShort <- as.list(gsub(pattern = gsubpattern, "\\1", nameVector))
    return(as.vector(mapply(
        function(x, y) rep(x, times = y), nameVectorShort,
        lengthList
    )))
}
