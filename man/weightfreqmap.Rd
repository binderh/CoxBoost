\name{weigthfreqmap}
\alias{weigthfreqmap}
\title{Plots weight-frequency maps from resCoxBoost fits}
\description{
Plots the weight-frequency map, i.e. a heat map for the different weights and different stably selected covariates  as obtained from a resCoxBoost object fitted by \code{\link{resCoxBoost}}.
}
\usage{
weightfreqmap<-function(RIF,mix.list=c(0.001, 0.01, 0.05, 0.1, 0.25, 0.35, 0.5, 0.7, 0.9, 0.99),
                        plotmix=c(0.001, 0.01, 0.05, 0.1, 0.25, 0.35, 0.5, 0.7, 0.9, 0.99),
                        lowerRIFlimit=0.5,method="complete")
}
\arguments{
\item{RIF}{list obtained from a \code{\link{resample.CoxBoost}} call.}
\item{mix.list}{vector of weights which were also entered in the \code{\link{resample.CoxBoost}} call.}
\item{plotmix}{vector of weights which should be plotted in the stability trajectory plot.}
\item{lowerRIFlimit}{cutoff RIF value: All covariates which have an resampling inclusion frequency (RIF) greater or equal this lowerRIFlimit value at any weight are presented in the stability trajectory plot.}
\item{method}{agglomeration method to be used for performing a hierarchical cluster analysis to get groups of covariates which are associated with different patient subgroups. See the \code{hclust} function for further details of the different methods. Default is \code{method="complete"}}
}
\details{The \code{weightfreqmap} function generates a visualization tool (weight-frequency map), which is a heatmap that plots the RIFs for stably selected covariates selected via \code{resample.CoxBoost}. This visualization tool is a useful tool for identifying groups of covariates which depend on the subgroup that is in the focus, on the other subgroup or on both subgroups. In this tool light gray color stands for a good RIF whereas darker gray color indicates a lower RIF and therefore decreased variable selection stability. The number of plotted covariates in this plot depends on \code{lowerRIFlimit}. How many weights are plotted which are element of \code{mix.list} can be changed with \code{plotmix}.
}

\value{
No value is returned, but the weight-frequency map is generated.
}
\author{
Veronika Weyer \email{weyer@uni-mainz.de} and Harald Binder \email{binderh@uni-mainz.de}
}