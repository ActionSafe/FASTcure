\name{plotpredictsmcure}
\alias{plotpredictsmcure}
\title{
Plot predicted smcure object
}
\description{
plot predicted survival curve(s) from the estimated mixture cure model
}
\usage{
plotpredictsmcure(object, type = "S", xlab = "Time",
ylab = "Predicted Survival Probability",
model = c("ph", "aft"), ...)
}
\arguments{
   \item{object}{
an object of the predictsmcure function
}
  \item{type}{
type of plot. "S" means steps plot.
}
  \item{xlab}{
a label for the x axis
}
  \item{ylab}{
a label for the y axis
}
  \item{model}{
either "ph" or "aft"}
 \item{\dots}{
Further options in plotfunction can be passed to the plotpredictsmcure function}
}
