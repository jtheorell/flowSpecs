# correctUnmix
context("specMatCalc")
data(unmixCtrls)
testUnmixSet <- unmixCtrls[c(1,2,12,15)]

testSpecMat <- specMatCalc(testUnmixSet, groupNames = c("Beads_"),
                       autoFluoName = "PBMC_unstained.fcs")

test_that("massOrFlowTransResult", {
    expect_equal(1, testSpecMat[1,"R2-A"])
})


