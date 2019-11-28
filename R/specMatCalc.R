#' Calculating the matrix used for spectral unmixing
#'
#'
#' This algoritm takes a flowSet containing single-stained controls and
#' negative controls, including an autofluorescence control and estimates the
#' unmixing for all fluorescent variables.
#' @param unmixCtrls A flowSet containing all the single stained and
#' unstained files necessary to create an spectral unmixing matrix.
#' @param groupNames A character vector containing strings common to the groups
#' of non-autofluoresence unmixCtrls that could be present. If for example
#' all antibodies single stains are anti-mouse bead-based the dead cell marker
#' is stained PBMC, and the files congruently either have a prefix containing
#' "Bead" or "PBMC", then the vector should be c("Bead", "PBMC"). The system
#' is not case specific.
#' @param autoFluoName The sample name of the autofluorescence control.
#' @return A data frame with each row representing a fluorochrome or
#' or autofluorescence and each column representing a detector.
#' @importFrom BiocGenerics colnames ncol
#' @examples
#' # Load suitable unmixing controls. NB! If these originate from different
#' # sample types, such as beads and PBMC, there should be a negative control
#' # for each group and the names should reflect this, so that all PBMC samples
#' # would be called PBMC_unstained, PBMC_DCM, etc.
#' data(unmixCtrls)
#'
#' # If the dataset contains cell controls, make sure that the cell population
#' # interest dominates FSC-A, as the data highest peak in this channel will be
#' # used.
#'
#' # And run the function
#' specMat <- specMatCalc(unmixCtrls, groupNames = c("Beads_", "Dead_"),
#' autoFluoName = "PBMC_unstained.fcs")
#' @export specMatCalc
specMatCalc <- function(unmixCtrls, groupNames, autoFluoName) {


    # The spectrum for each file is calculated

    specCalcMat <- fsApply(unmixCtrls, specCalc)

    # Now, the samples are categorized into groups depending on their sample
    # type reflected in the names of the samples. If any samples are singlets,
    # then they are put to the side. The most likely reason for having one
    # singlet is that all unmixing controls have been acquired with beads,
    # but that the autofluorescence control is unstained cells.

    singleStainGroupsList <- lapply(groupNames, function(x)
        return(specCalcMat[which(grepl(x, row.names(specCalcMat))), ]))

    # Now, it is checked that if there is only one group, this group contains
    # more than samples, to make sure that unmixing is not attempted with one
    # color only, as this leads to some downstream negative effects, and also
    # is more or less meaningless.

    if(length(singleStainGroupsList) == 1 &&
       nrow(singleStainGroupsList[[1]]) < 3){
           stop("It seems like unmixing of one color is attempted, which is
           not meaningful. Please try again with another set of controls,
           or name them differently, if you believe that you have more than
           one color. ")
       }


    # Now, in each matrix in the list, the row with the lowest sum
    # is identified as the unstained
    negCtrlRows <- lapply(
        singleStainGroupsList,
        function(x) which.min(rowSums(x))
    )

    # If the autoFluoName is not in the

    # Here, the subtractions are made
    rawSpecMatList <- lapply(seq_along(negCtrlRows), function(x) {
        localSpecMat <- apply(singleStainGroupsList[[x]], 1, function(y)
            y - singleStainGroupsList[[x]][negCtrlRows[[x]], ])
    })

    # Here, the data is coerced into a matrix
    rawSpecMat <- do.call("cbind", rawSpecMatList)

    # Now, the unstained controls are removed
    specMatNoUnstain <- rawSpecMat[, -which(colSums(rawSpecMat) == 0)]

    # Here, the column names are cleaned up.
    specMatColNamesRaw1 <- colnames(specMatNoUnstain)
    specMatColNamesRaw2 <- gsub("|\\.fcs", "", specMatColNamesRaw1)
    specMatColNames <- vector()
    for (i in specMatColNamesRaw2) {
        for (j in groupNames) {
            if (grepl(j, i)) {
                specMatColNames[i] <- gsub(paste0(j, "|"), "", i)
            }
        }
    }

    colnames(specMatNoUnstain) <- specMatColNames

    # Now, the autofluorescence medians are added
    specMat <- cbind(specMatNoUnstain,
        "Autofluo" = specCalcMat[autoFluoName, ]
    )

    specMatFrac <- t(apply(specMat, 2, function(x) x / max(x)))

    # And finally, all negative values resulting from minor errors in detection,
    # are removed.
    specMatFrac[which(specMatFrac < 0)] <- 0

    return(specMatFrac)
}

specCalc <- function(flowFrame) {
    focusColNames <- BiocGenerics::colnames(flowFrame)

    # First, a gate is applied to FSC.A, to simplify work with cells
    fscVar <- which(grepl("FSC", focusColNames) &
        grepl("A", focusColNames))

    fscGatedFrame <- madFilter(flowFrame, gateVar = fscVar, nMads = 1.5)
    fscFilteredFrame <- filterOut(fscGatedFrame,
        filterName = "FSC-A_auto_filter"
    )[
        , seq(1, BiocGenerics::ncol(flowFrame))
    ]

    # Now, a gate is applied to ssc, to clean up all files.
    sscVar <- which(grepl("SSC", focusColNames) &
        grepl("A", focusColNames))[1]

    sscGatedFrame <- madFilter(fscFilteredFrame, gateVar = sscVar, nMads = 1.5)
    sscFilteredFrame <- filterOut(sscGatedFrame,
        filterName = "SSC-A_auto_filter"
    )[
        , seq(1, BiocGenerics::ncol(flowFrame))
    ]

    # Here, all non-fluorescent channels are excluded
    fluoFrame <- sscFilteredFrame[, -which(grepl("ime", focusColNames) |
        grepl("SC", focusColNames))]
    # Then the median is calculated for all fluorescence channels on this
    # filtered population
    fluoColNames <- BiocGenerics::colnames(fluoFrame)

    rawMedVals <- vapply(fluoColNames, function(x)
        median(exprs(fluoFrame[, x])), 1)


    # Now, the highest peak is identified, and the data further gated on this
    # variable, to reduce the variance further
    maxMedVar <- which.max(rawMedVals)

    # This channel is now produced separately, to increase computational
    # speed
    maxMedVarFrame <- fluoFrame[, maxMedVar]

    # And here, this channel is transformed for the madFilter to work correctly
    maxMedVarFrameTrans <- arcTrans(maxMedVarFrame,
        transNames = colnames(maxMedVarFrame),
        transCoFacs = 400
    )

    # Here, the data is gated
    maxGatedFrame <- madFilter(maxMedVarFrameTrans,
        gateVar = 1,
        nMads = 1.5, returnSepFilter = TRUE
    )
    maxFilteredFrame <- fluoFrame[which(maxGatedFrame == 1), ]

    # And finally, the median procedure is repeated for this final population
    resultMedVals <- vapply(fluoColNames, function(x)
        median(exprs(maxFilteredFrame[, x])), 1)
}
