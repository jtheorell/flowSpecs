% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correctUnmix.R
\name{correctUnmix}
\alias{correctUnmix}
\title{Correct defects in spectral unmixing by compensation}
\usage{
correctUnmix(unmixFlowObj, corrMat, transCoFacs = 400)
}
\arguments{
\item{unmixFlowObj}{A flowframe or flowset post unmixing.}

\item{corrMat}{A correction matrix. If this is the first round, the
executionof this function needs to be preceeded by the generation of this
matrix, for example by using the \code{\link{corrMatCreate}} function.}

\item{transCoFacs}{If transformation should be performed, the transformation
cofactors can be added here. Three possible inputs: a vector with specific
cofactors for each variable, a set value that will be used for all variables,
and FALSE.
Note: It might be good to set this to FALSE in the final round, to optimize
the transoformations externally.}
}
\value{
The unmixed flow object, now corrected with the values from the
corrMat.
}
\description{
This function provides a way to reduce the defects in the spectral unmixing,
by creating a secondary correction matrix, which is symmetrical.
}
\examples{
# Load uncompensated data
data(fullPanel)

# Load the spectral unmixing matrix generated with controls from the same
# experiment. These can be generated using the specMatCalc function.
data(specMat)

# And now unmix
fullPanelUnmix <- specUnmix(fullPanel, specMat)


# Create an empty unmixinng matrix
corrMat <- corrMatCreate(specMat)

# Now correct the data with this. In the first instance, this will of course
# not have any effect, more  than transformation, as the corrMat is empty.
fullPanelCorr <- correctUnmix(fullPanelUnmix, corrMat)

# This now needs to be investigated, to identify any possible compensation
# defects. This is most easily done with the oneVsAllPlot executed in the
# following way:
\dontrun{
oneVsAllPlot(fullPanelCorr)
}
# One obvoius defect that shows when doing this is between CD56 and IgM:
oneVsAllPlot(fullPanelCorr, "BV650_CD56", saveResult = FALSE)

# This is corrcted the following way:
corrMat["BV650_CD56", "AF647_IgM"] <- -0.03
fullPanelCorr <- correctUnmix(fullPanelUnmix, corrMat)
oneVsAllPlot(fullPanelCorr, "BV650_CD56", saveResult = FALSE)

# This process is iterated until there are no remaining artifacts. Good help
# to do this is a set of fluorescence-minus-one controls. If that is not
# available, a rule of thumb is that if the signal in marker x is
# strongly negatively correlated to marker y, so that highly
# single-x-posisive values are below zero, then this is with all likelihood
# an artifact. The situation becomes more complicated with strong positive
# correlations, as they can occur in biology, so there one has to take more
# care and keep the marker biology in mind.

}
\seealso{
\code{\link{specUnmix}}, \code{\link{arcTrans}},
\code{\link{corrMatCreate}}
}
