#' Normalize batch differences in intensities by aligning peaks
#'
#'
#' This function is intended to be used on standardized controls,
#' preferrably one standard that has been acquired with every batch.
#' This function needs to be applied separately for each batch, with the same
#' standard.
#'
#' @param fs The flowSet to be normalized
#' @param ctrlPos The position in the flowSet that contains the internal
#' control
#' @param standFF The external standard to which all batches should be
#' normalized
#' @param transNames If parameters are not transformed prior to running this,
#' internal transformation is necessary to identify peaks, so specify which
#' variables that need transformation here. The data for these variables
#' are untransformed at the end, so the data will have the same scale in and
#' out.
#' @param transCoFacs Also inherited from arcTrans. Low values (2-10) for CyTOF
#' data, high values (200-2000) for flow data, all depending on the number of
#' input channels.
#' @param volThresh The threshold for how small the area of a peak can be compared to
#' the largest peak, and still count.
#' @param normOrNot This needs to be a logical vector with the same length as the
#' number of columns in the fs and standFF. If it is known that a certain
#' variable should not be normalized, that information can be specified
#' here.
#' @return A normalized flowSet. In addition, output to console, to clarify
#' #for which markers the same number of peaks was detected in the control and
#' the standard, as normalization will only be applied for these.
#' @importFrom flowCore sampleNames<-
#' @examples
#' #Load uncompensated data and spectral matrix.
#' data(fullPanel)
#' data(specMat)
#'
#' # And now unmix
#' fullPanelUnmix <- specUnmix(fullPanel, specMat)
#'
#'  #Create a new file with the value 1000 added to all values
#'  library(flowCore)
#'  fullPanelPlus1000 <- flowFrame(exprs(fullPanelUnmix)+1000)
#' # Check how they differ.
#' range(exprs(fullPanelUnmix)[,1])
#' # 143 1187733
#' range(exprs(fullPanelPlus1000)[,1])
#' # 1143 1188733
#'
#' #Now normalize the new one to the old. NB! Here we will only
#' normPanel1000 <- peakNorm(flowSet(fullPanelPlus1000), 1, fullPanelUnmix,
#' transNames = colnames(fullPanelUnmix)[6:18], transCoFacs = 500)
#'
#' #And now check the new result
#' range(exprs(normPanel1000[[1]] )[,1])
#' # 143 1187733
#' @export peakNorm
#'
peakNorm <- function(fs, ctrlPos, standFF, transNames = FALSE, transCoFacs,
                     volThresh = 0.005, normOrNot = rep(TRUE, ncol(standFF))){

    #First, the data is centered to the median in the standFF or the ctrlPos
    standMed <- apply(exprs(standFF), 2, median)
    ctrlMed <- apply(exprs(fs[[ctrlPos]]), 2, median)
    standFFCent <- medCent(standFF, standMed)
    fsCent <- fsApply(fs, medCent, ctrlMed)
    message("median centering done")

    #Now we make a rough arcTrans of the data, to get peaks
    if(length(transNames) == 1 && transNames == FALSE){
        standTrans <- standFFCent
        fsTrans <- fsCent
    } else {
        standTrans <- arcTrans(standFFCent, transNames, transCoFacs)
        fsTrans <- arcTrans(fsCent, transNames, transCoFacs)
    }

    #And now we normalize data. First, we identify the values in the control
    #frame, and then we apply them to all values in the flowset.
    peakVals <- lapply(
        list(exprs(standTrans),
             exprs(fsTrans[[ctrlPos]])), function(x){
                 locRes <- apply(x, 2, function(y){
                     locPeak <- peakIdenti(y, volThresh,
                                                       returnStats = TRUE)
                     if(length(locPeak$PeakPos) == 1){
                         peakPoss <- 1
                     } else {
                         peakPoss <- c(1,2)
                     }
                     locPeak$Median <- unlist(lapply(peakPoss, function(z){
                         median(y[
                             which(y > locPeak$Width[[z]][1] &
                                       y < locPeak$Width[[z]][2])])
                     }))
                     locPeak
                 })
             })

    fsFinal <- flowSet(lapply(seq_along(fsTrans), function(x){
        peakNormCoFunc(fsTrans[[x]], peakVals[[2]], peakVals[[1]],
                       standMed, inFF = fs[[x]], normOrNot,
                       transNames, transCoFacs)
    }))
    sampleNames(fsFinal) <- sampleNames(fs)
    fsFinal
}
#########
medCent <- function(FF, cent){
    locExprs <- exprs(FF)
    locResMedCentExprs <-
        do.call("cbind", lapply(seq_len(ncol(locExprs)),
                                function(y){
                                    yRes <- locExprs[,y] - cent[y]
                                }))
    colnames(locResMedCentExprs) <- colnames(locExprs)
    exprs(FF) <- locResMedCentExprs
    return(FF)
}
#########
peakNormCoFunc <- function(FF, ctrlPeakVals, standPeakVals, standMed, inFF,
                           normOrNot, transNames, transCoFacs){
    locRes <- lapply(seq_len(ncol(exprs(FF))), function(x){
        locDat <- exprs(FF)[,x]
        ctrlPeaks <- ctrlPeakVals[[x]]
        standPeaks <- standPeakVals[[x]]
        peaksMatch <- TRUE
        if(length(ctrlPeaks$Median) == 1 &&
           length(standPeaks$Median) == 1){
            subtrDat <- locDat - ctrlPeaks$Median
            normDat <- subtrDat + standPeaks$Median

        } else if(length(ctrlPeaks$Median) == 2 &&
                  length(standPeaks$Median) == 2){
            subtrDat <- (locDat - ctrlPeaks$Median[1])/
                (ctrlPeaks$Median[2] - ctrlPeaks$Median[1])
            normDat <- (subtrDat * (standPeaks$Median[2] -
                                        standPeaks$Median[1])) +
                standPeaks$Median[1]

        } else {
            normDat <- locDat
            peaksMatch <- FALSE
        }
        return(list(normDat, peaksMatch))

    })

    locResMat <- do.call("cbind", lapply(locRes, "[[", 1))
    colnames(locResMat) <- colnames(FF)
    locNormOrNot <- unlist(lapply(locRes, "[[", 2))
    for(i in seq_along(normOrNot)){
        if(locNormOrNot[i] == FALSE){
            normOrNot[i] <- FALSE
        }
    }
    names(normOrNot) <- colnames(FF)
    message(paste0(names(normOrNot), "=", normOrNot, " "))

    exprs(FF) <- locResMat
    #Now we can untransform the data again and bring it back to its
    #original values.
    if(length(transNames) == 1 && transNames == FALSE){
        FFUntrans <- FF
    } else {
        FFUntrans <- arcTrans(FF, transNames, transCoFacs, unTrans = TRUE)
    }

    #And now finally we reverse the median subtraction from the first step.
    FFOrig <- medCent(FFUntrans, -standMed)

    #Now, for the columns where there was a mismatch in the number of peaks
    #between the standard and the internal control, the original, untouched
    #data is returned
    exprs(FFOrig)[,which(normOrNot == FALSE)] <-
        exprs(inFF)[,which(normOrNot == FALSE)]

    return(FFOrig)
}
