% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DAMisc_functions.R
\name{simPredpolr}
\alias{simPredpolr}
\title{Calculate Predictions for Proportional Odds Logistic Regression}
\usage{
simPredpolr(object, coefs, n.coef)
}
\arguments{
\item{object}{An object of class \code{polr}.}

\item{coefs}{A vector of coefficients where elements 1 to \code{n.coef} give
model coefficients and elements \code{n.coef}+1 to k have intercepts.}

\item{n.coef}{Number of coefficients (minus intercepts) for the \code{polr}
model.}
}
\value{
An n x m-category matrix of predicted probabilities
}
\description{
Calculates predicted probabilities from models of class \code{polr} from a
model object and a vector of coefficient values.  This is an auxiliary
function used in \code{pre} if \code{sim=TRUE}.
}
\author{
Dave Armstrong
}
