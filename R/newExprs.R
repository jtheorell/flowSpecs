#' This function lets us exchange a flowframes exprs portion to an unrelated one
#' It is solely meant to be used internally, as it is a strange practice.
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
