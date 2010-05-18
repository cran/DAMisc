\encoding{UTF-8}
\name{factorplot}
\alias{factorplot}
\title{Exact Factor Contrast Differences}
\description{
This function calculates and presents graphically all simple contrasts for factor covariates in linear models and GLMs.  
}
\usage{
factorplot(obj, factor.variable=NULL, pval = 0.05, two.sided = TRUE, bonferroni = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object of class 'lm' or 'glm'}
  \item{factor.variable}{String containing the name of the factor to be plotted}
  \item{pval}{The (uncorrected) Type I error probability required, defaults to 0.05}
  \item{two.sided}{Logical indicating whether two-sided hypothesis tests should be conducted.  If FALSE (not the default), one-sided tests in the observed direction are computed.}
  \item{bonferroni}{Should a bonferroni adjsustment to the p-value be employed to correct for simultaneous testing.  If TRUE (not the default), the p-value is changed from p to p/k, where k is the number of contrasts being tested. }
}
\details{
This function calculates and presents graphically all simple contrasts for factor covariates in linear models and GLMs.  Specifically, it calculates:\cr\cr
 \deqn{\frac{b_{row}-b_{col}}{\sqrt{V(b_{row}-b_{col})}}}{(b[row]-b[col])/V(b[row]-b[col])^(1/2)}\cr where \deqn{V(b_{row}-b_{col}) = V(b_{row}) + V(b_{col}) + 2(1)(-1)V(b_{row},b_{col})}{V(b[row]-b[col]) = V(b[row]) + V(b[col]) + 2(1)(-1)V(b[row],b[col])}\cr\cr  This function is an alternative to Firth's quasi-variances.  The beauty of the quasi-variances is that they can be used to calculate any contrast.  Drawbacks include a presentation that still requires calculation to perform hypothesis tests and more problematically, the possibility that the underlying model produces negative quasi-variances.  While the function here is less flexible, it provides exact hypothesis tests directly for simple contrasts, even when the quasi-variances are negative. \cr\cr
The function automatically scales the resulting graph to the appropriate size, but will be less useful as the number of contrasts gets large (e.g., > 30).  If more than one factor covariate is present, the function generates a text-based menu in the R GUI that will allow the users to pick the term for which they want a graph.    
}
\value{
\item{graph}{A graph is produced, though no other output is currently saved.}
}
\references{
Firth, D. (2000) Quasi-variances in Xlisp-Stat and on the web. Journal of Statistical Software 5.4, 1–13. \cr
\cr
Firth, D. (2003) Overcoming the reference category problem in the presentation of statistical models. Sociological Methodology 33, 1–18.}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}
\examples{
data(FactorplotEx)
mod <- lm(y ~ x1 + x2, data=FactorplotEx)
factorplot(mod, factor.variable='x1', 0.05, TRUE, FALSE)
}
