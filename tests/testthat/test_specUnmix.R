# correctUnmix
context("specUnmix")

testSpecMat <- matrix(c(1,0.5,0,1), ncol = 2, nrow = 2,
                      dimnames = list(c("Res1", "Res2"),
                                      c("Start1", "Start2")))
testFlowFrame <- flowCore::flowFrame(matrix(c(1,1,1,1),
                                            nrow = 2, ncol = 2,
                                            dimnames = list(c(NA, NA),
                                                            c("Start1",
                                                              "Start2"))))
resultExprs <- exprs(specUnmix(testFlowFrame, testSpecMat))

test_that("specUnmixResult1", {
    expect_equal("Res1", colnames(resultExprs)[1])
})

test_that("specUnmixResult2", {
    expect_equal(0.5, resultExprs[1,1])
})

test_that("specUnmixResult3", {
    expect_equal(1, resultExprs[1,2])
})


