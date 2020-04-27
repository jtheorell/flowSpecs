#' Peak identification for higher-level functions.
#'
#' This function is primarily thought to be used internally to define peaks
#' in data.
#' One function is borrowed from package vulcan, namely the vulcan::densityauc,
#' which is neat, but the package is large and significantly increases the
#' installation time, and the importing is thus discarded.
#' @param markerData The data that the peaks should be identified for
#' @param volThresh The cutoff ratio of the volume for each secondary peak,
#' under which it is not considered to be a peak
#' @param distThresh The cutoff under which two peaks are considered one, as
#' they are too close to each other. This value between 0 and 1 corresponds to a
#' fraction from the 10th to the 90th percentile of the data range that the
#' peaks must be separated by to count. Defaults to 0.1 or 10 percent of the
#' distance.
#' @param adjust The value deciding the accuracy of the density calculation. The
#' higher the value, the lower the sensitivity for small aberrations in the
#' density.
#' @param nPeaks The number of peaks that should be exported. If n+1
#' fulfilling the volRatio criterion are found, the peaks most separated in
#' space are chosen.
#' @param returnStats Should the deflection points defining the peaks, the peak
#' hight and the lowest deflection point between the two most extreme peaks
#' be included in the export?
#' @return The information about the peaks in question. Depending on if
#' returnStats is TRUE or not and the number of peaks, it will change in
#' complexity.
#' @seealso \code{\link[vulcan]{densityauc}}
#' @keywords internal
#' @importFrom zoo rollmean
peakIdenti <- function(markerData, volThresh = 0.05, distThresh = 0.1,
                       adjust = 2, nPeaks = 2, returnStats = FALSE) {
    #First, we remove outliers, as they can have unintended effects
    #especially on the distThresh, as true clusters can look very close if a
    #few extreme values are included.
    quantile1_99 <- quantile(markerData, c(0.001, 0.999))
    markerData <- markerData[which(markerData > quantile1_99[1] &
                                       markerData < quantile1_99[2])]
    Da <- density(markerData, adjust = adjust, na.rm = TRUE)
    DeltaY <- diff(Da$y)
    Turns <- which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1

    # In the case that the peak is located at the lowest position
    # in the density Y vector, it is included in the turns vector here
    if (which.max(Da$y) == 1) {
        Turns <- c(1, Turns)
    }

    # Here, the top peak decides if the peaks are odd or even numbers in Turns.
    maxPeakPos <- which.max(Da$y[Turns])
    if (maxPeakPos %% 2 == 0) {
        allPeakMaxima <- Turns[seq_along(Turns) %% 2 == 0]
    } else {
        allPeakMaxima <- Turns[seq_along(Turns) %% 2 != 0]
    }
    # Either, there is only one peak, and in that case, no further analyses are
    # needed. If not, a set of functions will find which peaks to choose
    if (length(allPeakMaxima) == 1) {
        if (returnStats) {
            return(list(
                "PeakPos" = Da$x[allPeakMaxima],
                "Width" = list(c(min(markerData), max(markerData))),
                "DensVol" = NA,
                "Height" = NA, "LowVertex" = NA
            ))
        } else {
            return(Da$x[allPeakMaxima])
        }
    }

    # Now, the volume for each peak is calculated.
    # Each peak is defined as the region from the preceding to the succeding
    # deflection point
    densVols <- vector()
    peakRegVals <- list()
    peakHeights <- vector()
    for (i in seq_along(allPeakMaxima)) {
        focusTurnPos <- which(Turns == allPeakMaxima[i])
        if (focusTurnPos == 1) {
            peakRegion <- c(1, Turns[2])
        } else if (focusTurnPos == length(Turns)) {
            peakRegion <- c(Turns[focusTurnPos - 1], length(Da$x))
        } else {
            peakRegion <- c(Turns[focusTurnPos - 1], Turns[focusTurnPos + 1])
        }
        localPeakRegVals <- c(Da$x[peakRegion[1]], Da$x[peakRegion[2]])

        # Now calculate the volume of the region in question
        xt <- diff(Da$x[Da$x > localPeakRegVals[1] & Da$x <
            localPeakRegVals[2]])
        yt <- rollmean(
            Da$y[Da$x > localPeakRegVals[1] & Da$x <
                localPeakRegVals[2]],
            2
        )
        densVols[i] <- sum(xt * yt)
        peakRegVals[[i]] <- localPeakRegVals
        peakHeights[i] <- Da$y[allPeakMaxima[i]]
    }

    # Here, the volumes are sorted by size
    allPeakMaximaSorted <- allPeakMaxima[order(densVols, decreasing = TRUE)]
    peakRegValsSorted <- peakRegVals[order(densVols, decreasing = TRUE)]
    densVolsOrdered <- densVols[order(densVols, decreasing = TRUE)]
    peakHeightsSorted <- peakHeights[order(densVols, decreasing = TRUE)]

    # Now, all peaks with a volume smaller than volThresh are excluded
    peakMaximaReal <- allPeakMaximaSorted[which(densVolsOrdered > volThresh)]
    peakRegValsReal <- peakRegValsSorted[which(densVolsOrdered > volThresh)]
    densVolsReal <- densVolsOrdered[which(densVolsOrdered > volThresh)]
    peakHeightsReal <- peakHeightsSorted[which(densVolsOrdered > volThresh)]

    # Now, if this returns only one value, this value is returned here
    if (length(peakMaximaReal) == 1 || nPeaks == 1) {
        if (returnStats) {
            return(list(
                "PeakPos" = Da$x[peakMaximaReal[1]],
                "Width" = peakRegValsReal[1],
                "DensVol" = densVolsReal[1],
                "Height" = peakHeightsReal[1],
                "LowVertex" = NA
            ))
        } else {
            return(Da$x[peakMaximaReal][1])
        }
    }

    # Here, it is investigaded if the peaks are separated by more than
    # a set fraction of the total width of the data. This is done
    # iteratively, first considering the largest and the second largest
    # peak, then considering if the third peak is separated from both the
    # first and the second peak, and so on. The intention with this
    # algorithm is to avoid detecting a peak with a minor aberration at the top
    # to be detected as two peaks.

    threshDist <- length(Da$x) * distThresh
    peakMaximaFinal <- peakMaximaReal[1]
    j <- 2
    for (i in seq(2, length(peakMaximaReal))) {
        if (all(vapply(
            seq(1, i - 1), function(x)
                abs(peakMaximaReal[x] - peakMaximaReal[i]) > threshDist,
            TRUE
        ))) {
            peakMaximaFinal[j] <- peakMaximaReal[i]
            j <- j + 1
        }
    }

    peakRegValsFinal <-
        peakRegValsReal[which(peakMaximaReal %in% peakMaximaFinal)]
    densVolsFinal <- densVolsReal[which(peakMaximaReal %in% peakMaximaFinal)]
    peakHeightsFinal <-
        peakHeightsReal[which(peakMaximaReal %in% peakMaximaFinal)]

    # Here, we are identifying the lowest vertex present between any of the
    # peaks
    dayPeakSection <- Da$y[seq(
        min(peakMaximaFinal),
        max(peakMaximaFinal)
    )]

    lowVertex <- dayPeakSection[which.min(dayPeakSection)]

    # Now, the exports of more complex solutions is generated
    if (nPeaks < length(peakMaximaFinal)) {
        mostSepVals <- vapply(peakMaximaFinal, function(x)
            sum(abs(x - peakMaximaFinal)), 1)

        reducedPositions <- sort(peakMaximaFinal[order(mostSepVals,
            decreasing = TRUE
        )][seq_len(nPeaks)])

        peakMaximaFinalReduced <- peakMaximaFinal[which(peakMaximaFinal
        %in% reducedPositions)]

        # And here comes the reduction of all the metadata
        peakRegValsFinalReduced <-
            peakRegValsFinal[which(peakMaximaFinal %in% reducedPositions)]
        densVolsFinalReduced <- densVolsFinal[which(peakMaximaFinal
        %in% reducedPositions)]
        peakHeightsFinalReduced <- peakHeightsFinal[
            which(peakMaximaFinal %in% reducedPositions)
        ]

        # Followed by sorting into the natural order
        peakMaximaFinalSorted <-
            peakMaximaFinalReduced[order(peakMaximaFinalReduced)]

        peakRegValsFinalSorted <-
            peakRegValsFinalReduced[order(peakMaximaFinalReduced)]

        densVolsFinalSorted <-
            densVolsFinalReduced[order(peakMaximaFinalReduced)]

        peakHeightsFinalSorted <-
            peakHeightsFinalReduced[order(peakMaximaFinalReduced)]
    } else {
        peakMaximaFinalSorted <- peakMaximaFinal[order(peakMaximaFinal)]
        peakRegValsFinalSorted <- peakRegValsFinal[order(peakMaximaFinal)]
        densVolsFinalSorted <- densVolsFinal[order(peakMaximaFinal)]
        peakHeightsFinalSorted <- peakHeightsFinal[order(peakMaximaFinal)]
    }

    if (returnStats) {
        return(list(
            "PeakPos" = Da$x[peakMaximaFinalSorted],
            "Width" = peakRegValsFinalSorted,
            "DensVol" = densVolsFinalSorted,
            "Height" = peakHeightsFinalSorted,
            "LowVertex" = lowVertex
        ))
    } else {
        return(Da$x[peakMaximaFinalSorted])
    }
}
