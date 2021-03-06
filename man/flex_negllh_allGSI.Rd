% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flex_negllh_allGSI.R
\name{flex_negllh_allGSI}
\alias{flex_negllh_allGSI}
\title{calculate the negative log-likelihood with only PBT and GSI (optional) information}
\usage{
flex_negllh_allGSI(params, nPBT, nGSI, ohnc, t, utGSI, ohnc_gsi)
}
\arguments{
\item{params}{vector of paramaters to optimize}

\item{nPBT}{number of PBT groups to estimate}

\item{nGSI}{number of GSI groups to estimate}

\item{ohnc}{vector of number of observed (PBT assigned) hatchery fish in each PBT and GSI group (gsi groups should be 0)}

\item{t}{vector of tag rates for all PBT and GSI groups (gsi groups should be 0)}

\item{utGSI}{vector of number of un-PBT assigned fish in each GSI group}

\item{ohnc_gsi}{matrix of counts of fish GSI assigned to various groups}
}
\description{
I expect that most users interested in MLE of composition will be best served by using \code{MLEwrapper},
which preps data, calls this function, and then returns the results in a nice format. This function calculates the
likelihood for a given set of parameters, and may be useful to users who want MLE but want more control over how the
data is used and/or the optimization procedure. This function does not allow another variable. for that,
see \code{flex_negllh_var}. This calculates a simplified version of the negative
log-likelihood (constants removed) for computational efficiency.
}
