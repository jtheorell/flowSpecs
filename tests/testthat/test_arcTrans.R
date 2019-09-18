# arcTrans
context("arcTrans")
value <- 5
coFactor <- 5
result <- asinh(value/coFactor)
valueFrame <- flowCore::flowFrame(matrix(value, 1, 1, dimnames = list(NA, "Val")))
transFrame <- arcTrans(valueFrame, transNames = colnames(valueFrame),
                       transCoFacs = coFactor)
test_that("arcTransResult", {
    expect_equal(as.numeric(exprs(transFrame)), result)
})

#Untransform the data again
unTransFrame <- arcTrans(transFrame, transNames = colnames(valueFrame),
                         transCoFacs = coFactor, unTrans = TRUE)
test_that("arcUnTransResult", {
    expect_equal(as.numeric(exprs(unTransFrame)), value)
})
