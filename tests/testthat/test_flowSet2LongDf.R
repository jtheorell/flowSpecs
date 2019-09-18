# correctUnmix
context("flowSet2LongDf")
library(flowCore)
data(fullPanel)

testFrame <- fullPanel[seq(1,2),1]
exprs(testFrame)[,1] <- c(1,2)
testSet <- flowSet(testFrame[1,], testFrame[2,])
sampleNames(testSet) <- c("Test_1.fcs", "Test_2.fcs")
testDf <- flowSet2LongDf(testSet, idInfo = list("Number" = "Test_|\\.fcs"))

test_that("longDfRes1.1", {
    expect_equal(1, testDf[1,1])
})

test_that("longDfRes1.2", {
    expect_equal("1", testDf[1,2])
})

test_that("longDfRes2.1", {
    expect_equal(2, testDf[2,1])
})

test_that("longDfRes2.2", {
    expect_equal("2", testDf[2,2])
})

test_that("longDfResDate", {
    expect_equal("25-Oct-2018", testDf[1,3])
})

