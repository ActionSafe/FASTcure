\name{predictsmcure}
\alias{predictsmcure}
\title{
prediction of semicure model}
\description{
Prediction of semicure model}
\usage{
predictsmcure(object, newX, newZ, model = c("ph", "aft"), ...)
}
\arguments{
  \item{object}{
an object of smcure
}
  \item{newX}{
new value(s) of X}
  \item{newZ}{
new value(s) of Z}
  \item{model}{
either 'ph' or 'aft'}
 \item{\dots}{
further arguments to be passed to the predictsmcure function}
}
\details{
Predicted population survival function can be calculated by the following equation
\deqn{S_{pop}(t)=\pi+(1-\pi)*S(t)
 }
}