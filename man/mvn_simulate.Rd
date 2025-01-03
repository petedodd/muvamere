% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{mvn_simulate}
\alias{mvn_simulate}
\title{mvn_simulate}
\usage{
mvn_simulate(X, Beta, Sigma)
}
\arguments{
\item{X}{multivariable explanatory data - each line a set of predictors for a patient}

\item{Beta}{matrix of regression coefficients}

\item{Sigma}{covariance matrix for response}
}
\value{
a matrix of responses Y, with nrow(Y)==nrow(X) and ncol(Y)==ncol(Sigma)
}
\description{
MVN sample for a single study
}
\details{
MVN sample for a single study using explanatory variables X, regression coefficient matrix Beta, and covariance matrix Sigma
}
\author{
Pete Dodd
}
