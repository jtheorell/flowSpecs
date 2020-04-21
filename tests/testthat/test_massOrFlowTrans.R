# correctUnmix
context("massOrFlowTrans")
data(fullPanel)

testFrame <- flowFrame(matrix(c(rep(0, 50), 1),
                              nrow= 51, ncol = 1, dimnames  =
                                  list(rep(NA, 51), "V1")))
as.numeric(flowSpecs:::massOrFlowTrans(testFrame, transNames = "V1"))

test_that("massOrFlowTransResult", {
    expect_equal(5, as.numeric(
        flowSpecs:::massOrFlowTrans(testFrame, transNames = "V1")))
})


