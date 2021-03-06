% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mwlsr.R
\name{mwlsr.Ftest}
\alias{mwlsr.Ftest}
\title{mwlsr.Ftest}
\usage{
mwlsr.Ftest(fit)
}
\arguments{
\item{fit}{mwlsr fit object}
}
\value{
data.frame with F-test results
}
\description{
Calculates F-statistic and p-values for all models in the fit. Returns 
a table of the results.  This is only sensible if your design included
an intercept.
}
