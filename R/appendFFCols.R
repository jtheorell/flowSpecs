# This function is a slightly adjusted copy of the fr_append_cols
# function in flowCore. It mainly allows for the presence of additional
# pData variables
# @param focusFrame The flowFrame that will get appended.
# @param newCols The new columns that should be added.
# @return The focusFrame, with the newCols appended.
#' @importFrom flowCore pData exprs parameters
appendFFCols <- function(focusFrame, newCols) {
    pd <- flowCore::pData(flowCore::parameters(focusFrame))
    cn <- colnames(newCols)
    new_pid <- max(as.integer(gsub("\\$P", "", rownames(pd)))) +
        1
    new_pid <- seq(new_pid, length.out = ncol(newCols))
    new_pid <- paste0("$P", new_pid)

    new_pd <- do.call(rbind, lapply(cn, function(i) {
        vec <- newCols[, i]
        rg <- range(vec)
        new_pd <- pd[1, ]
        new_pd[1, ] <- NA
        new_pd$name <- i
        new_pd$range <- diff(rg) + 1
        new_pd$minRange <- rg[1]
        new_pd$maxRange <- rg[2]
        return(new_pd)
    }))
    rownames(new_pd) <- new_pid
    pd <- rbind(pd, new_pd)
    newExprs(focusFrame) <- cbind(exprs(focusFrame), newCols)
    flowCore::pData(flowCore::parameters(focusFrame)) <- pd
    focusFrame
}

#' This function lets us exchange a flowframes exprs portion to an unrelated one
#' It is solely meant to be used within appendFFCols.
#' @param x A flowFrame
#' @param value A matrix suitable to be an exprs object.
#' @return A new flowFrame with "value" as the exprs portion.
#' @importFrom methods setGeneric setMethod
#' @keywords internal
setGeneric("newExprs<-", function(x, value) standardGeneric("newExprs<-"))

setMethod("newExprs<-", "flowFrame", function(x, value) {
    x@exprs <- value
    x
})
