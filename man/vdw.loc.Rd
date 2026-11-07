\name{vdw.loc}
\alias{vdw.loc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Van der Waerden Estimator of Location}
\description{
  Computes the location estimator based on van der Waerden scores (sometimes also referred to as normal scores).
  }
\usage{
vdw.loc(x, int.diff = 10, maxiter = 1000, na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector.}
  \item{int.diff}{no longer used, kept for compatibility with earlier versions. }
  \item{maxiter}{no longer used, kept for compatibility with earlier versions.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}
\details{
  The estimator is the Hodges-Lehmann type estimator based on the signed-rank statistic with
  van der Waerden scores
  \deqn{T(\mu) = \sum_{i=1}^n sign(x_i - \mu) \Phi^{-1}\left(\frac{1}{2} + \frac{R_i^+(\mu)}{2(n+1)}\right),}{T(mu) = sum_i sign(x_i - mu) qnorm(1/2 + R_i(mu) / (2(n+1))),}
  where \eqn{R_i^+(\mu)}{R_i(mu)} is the rank of \eqn{|x_i - \mu|}{|x_i - mu|} among
  \eqn{|x_1 - \mu|, \ldots, |x_n - \mu|}{|x_1 - mu|, ..., |x_n - mu|}. The function \eqn{T}{T} is nonincreasing
  and constant between neighbouring Walsh averages \eqn{(x_i + x_j)/2}{(x_i + x_j)/2}, \eqn{i \le j}{i <= j},
  and the estimate is
  \deqn{\hat\mu = \frac{1}{2}\left(\sup\{\mu : T(\mu) > 0\} + \inf\{\mu : T(\mu) < 0\}\right).}{mu.hat = (sup{mu : T(mu) > 0} + inf{mu : T(mu) < 0}) / 2.}
  It is computed exactly by a bisection over the sorted Walsh averages.

  With Wilcoxon scores instead of van der Waerden scores the same definition gives the
  Hodges-Lehmann estimator \code{\link{hl.loc}}.

  Up to version 1.1-3 the estimate was obtained by a linear interpolation between two Walsh averages,
  which gave only an approximation of the estimate and could fail.
}
\value{
  the van der Waerden score estimator of location.
}
\references{
\cite{Hodges, J.L. and Lehmann, E.L. (1963), Estimates of location based on rank tests, \emph{The Annals of Mathematical Statistics}, \bold{34}, 598--611.}

\cite{Hettmansperger, T.P. and McKean, J.W. (1998), \emph{Robust Nonparametric Statistical Methods},
London, Arnold. }}
\author{Klaus Nordhausen}

\seealso{\code{\link{hl.loc}}}

\examples{
set.seed(1)
x <- rt(100, df = 3)
vdw.loc(x)
rm(.Random.seed)
}
\keyword{nonparametric}
\keyword{robust}
\keyword{univar}
