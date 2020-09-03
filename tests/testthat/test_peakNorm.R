data(fullPanel)
extract <- fullPanel[,c(2,20)]

extractPlus1000 <- flowFrame(exprs(extract)+1000)

#Now normalize the new one to the old.
normPanel1000 <- peakNorm(flowSet(extractPlus1000), 1, extract)

#And now check the new result

test_that("massOrFlowTransResult", {
    expect_equal(range(exprs(extract)), range(exprs(normPanel1000[[1]])))
})
