% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{mvn_simulate_studies}
\alias{mvn_simulate_studies}
\title{mvn_simulate_studies}
\usage{
mvn_simulate_studies(
  Xlist,
  betag,
  sigb,
  rhoA,
  rhoB,
  taug,
  sigt,
  lkj_local,
  lkj_global
)
}
\arguments{
\item{Xlist}{list of covariate data matrices}

\item{betag}{matrix global regression parameters mean}

\item{sigb}{matrix global regression parameters SD}

\item{rhoA}{beta prior parameter for interpolant between local and global correlations}

\item{rhoB}{beta prior parameter for interpolant between local and global correlations}

\item{taug}{global tau mean(s)}

\item{sigt}{global tau SD(s)}

\item{lkj_local}{local correlation LKJ prior parameter}

\item{lkj_global}{global correlation LKJ prior parameter}
}
\value{
a data frame with responses for all studies, with studyno and obsno columns appended
}
\description{
Simulate responses for a number of studies
}
\details{
Given a list of covariate data inputs, a corresponding set of responses is generated.
}
\author{
Pete Dodd
}
