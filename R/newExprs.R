#' This function lets us exchange a flowframes exprs portion to an unrelated one

setGeneric("newExprs<-", function(x, value) standardGeneric("newExprs<-"))

setMethod("newExprs<-", "flowFrame", function(x, value) {
    x@exprs <- value
    x
})
