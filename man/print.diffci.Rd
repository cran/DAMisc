% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DAMisc_functions.R
\name{print.diffci}
\alias{print.diffci}
\title{Print Confidence Intervals for Predicted Probabilities and Their Differences}
\usage{
\method{print}{diffci}(
  x,
  type = c("pr", "pw"),
  ...,
  digits = 4,
  filter = NULL,
  const = NULL,
  onlySig = FALSE
)
}
\arguments{
\item{x}{A object of class \code{diffci} produced by \code{\link{probci}}.}

\item{type}{Which kind of result to print - predictions (\code{pr}) or 
pairwise differences in predictions (\code{pw}).}

\item{...}{Other arguments to be passed down to print, currently
unimplemented.}

\item{digits}{How many digits to round output.}

\item{filter}{A named list of values where the names indicate the variable
to be filtered and the values in the vector indicate the values to include
for the filtering variable.}

\item{const}{A string identifying the name of the variable to be held
constant across comparisons.  Only applies if \code{type = "pw"}.}

\item{onlySig}{Logical indicating whether all differes should be displayed
or only those significant at the 95\% two-tailed level.}
}
\value{
An data frame with the following variables: \item{variables}{The
variables and the values at which they are held constant.  For example,
\code{tmp1} would be the first value of \code{tmp} used in the probability
calculation and \code{tmp2} would be the second value of \code{tmp} used in
the probability calculation. } \item{pred_prob}{The difference in predicted
probability given the following change in \code{X}:
\code{tmp2}-\code{tmp1}.} \item{lower, upper}{The lower and upper 95\%
confidence bounds.}
}
\description{
Print method for output from the \code{probci} function.
}
\examples{

data(france)
left.mod <- glm(voteleft ~ male + age + retnat + 
	poly(lrself, 2, raw=TRUE), data=france, family=binomial)
data(france)
left.mod <- glm(voteleft ~ male + age + retnat + 
                  poly(lrself, 2, raw=TRUE), data=france, family=binomial)
out <- probci(left.mod, france, numQuantVals=3, 
              changeX=c("retnat", "lrself"), calcPW = TRUE)
print(out, filter=list(retnat=c("Better", "Worse")))
print(out, type="pw", 
      filter=list(retnat=c("Better", "Worse")),
      const="lrself")
}
\author{
Dave Armstrong
}
