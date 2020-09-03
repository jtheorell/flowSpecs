# correctUnmix
context("correctUnmix")
valueFrame <-
    flowCore::flowFrame(matrix(c(1,1, 1, 1), 2, 2,
                     dimnames = list(c(NA, NA), c("V1", "V2"))))
corrMat <- matrix(c(1, -0.5, 0, 1), nrow = 2, ncol = 2,
                  dimnames = list(c("V1", "V2"), c("V1", "V2")))
locRes <- as.numeric(exprs(correctUnmix(
    valueFrame, corrMat, transCoFacs = FALSE))[1,1])
test_that("correctUnmixResult", {
    expect_equal(locRes, 1.5)
})


