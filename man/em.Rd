\name{em}
\alias{em}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EM algorithm}
\description{
EM algorithm used in the mixture cure model by assuming a latent indicator of uncure. Detailed estimation methods can be found in the reference section}
\usage{
em(Time, Status, X, Z, offsetvar, b, beta, model, link, emmax, eps)
}
\arguments{
  \item{Time}{
time to event of interest}
  \item{Status}{
status indicator, 0=alive, 1=dead
}
  \item{X}{
a vector or matrix of covariates corresponding to latency part
}
  \item{Z}{
a vector or matrix of covariates corresponding to incidence part
}
  \item{offsetvar}{
offset variable }
  \item{b}{
initial value for parameter b
}
  \item{beta}{
initial value for parameter beta
}
  \item{model}{
either "ph" or "aft"}
  \item{link}{
specifies the link in incidence part. The "logit", "probit" or complementary loglog ("cloglog") links are available. By default link = "logit".}
  \item{emmax}{
specifies the maximum iteration number. If the convergence criterion is not met, the EM iteration will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default emmax = 100.}
  \item{eps}{
sets the convergence criterion. The default is eps = 1e-7. The iterations are considered to be converged when the maximum relative change in the parameters and likelihood estimates between iterations is less than the value specified.
}
}
