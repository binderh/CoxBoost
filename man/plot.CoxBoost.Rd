\name{plot.CoxBoost}
\alias{plot.CoxBoost}
\title{Plot coefficient paths from CoxBoost fit}
\description{
Plots coefficient paths, i.e. the parameter estimates plotted against the boosting steps as obtained from a CoxBoost object fitted by \code{\link{CoxBoost}}.
}
\usage{
\method{plot}{CoxBoost}(x,line.col="dark grey",label.cex=0.6,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,\ldots)
}
\arguments{
\item{x}{fitted CoxBoost object from a \code{\link{CoxBoost}} call.}
\item{line.col}{color of the lines of the coefficient path}
\item{label.cex}{scaling factor for the variable labels.}
\item{xlab}{label for the x axis, default label when \code{NULL}.}
\item{ylab}{label for the y axis, default label when \code{NULL}.}
\item{xlim,ylim}{plotting ranges, default label when \code{NULL}.}
\item{\dots}{miscellaneous arguments, passed to the plot routine.}
}
\value{
No value is returned, but a plot is generated.
}
\author{
Harald Binder \email{binderh@uni-mainz.de}
}
