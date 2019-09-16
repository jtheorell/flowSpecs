# The point of this function is to enable autodetection of the "optimal" value
# for the arcsinh cofactor in the transformation of the data. This is because
# CyTOF and flow cytometry data have very different needs, but the data can also
# easily be differentiated, as CyTOF data is much sparser.
# focusFrame = the frame that will be transformed.
# transNames = the names of the columns in the focusFrame that should be
# transformed.

massOrFlowTrans <- function(focusFrame, transNames) {
    if (unname(sum(exprs(focusFrame) == 0) / (dim(focusFrame)[1] *
        dim(focusFrame)[2])) > 0.05) {
        transCoFacs <- rep(5, times = length(transNames))
    } else {
        transCoFacs <- rep(400, times = length(transNames))
    }
    names(transCoFacs) <- transNames

    return(transCoFacs)
}
