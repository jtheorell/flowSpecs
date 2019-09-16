#' This function is a co-function to the oneVsAllPlot function, to reduce
#' computational time.
#' @param flowObj The object to downsample
#' @param nRows The number of rows in the downsampled sample.
#' @return The downsampled flowObj.
#' @importFrom BiocGenerics nrow
plotDownSample <- function(flowObj, nRows = 10000) {
    if (inherits(flowObj, "flowSet")) {
        if (sum(fsApply(flowObj, BiocGenerics::nrow)) > nRows) {
            numIncluded <- floor(nRows / length(flowObj))
            plotExprs <- fsApply(flowObj, function(x)
                return(exprs(x)[sample(
                    seq_len(BiocGenerics::nrow(x)),
                    numIncluded
                ), ]))
        } else {
            plotExprs <- fsApply(flowObj, function(x) return(exprs(x)))
        }
    } else if (BiocGenerics::nrow(flowObj) > nRows) {
        plotExprs <-
            exprs(flowObj)[sample(
                seq_len(BiocGenerics::nrow(flowObj)),
                nRows
            ), ]
    } else {
        plotExprs <- exprs(flowObj)
        message(
            "The number of rows in the dataset was only ",
            nrow(plotExprs),
            ", so the output will be restricted to this"
        )
    }
    return(plotExprs)
}
