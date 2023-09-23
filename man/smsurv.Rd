\name{smsurv}
\alias{smsurv}
\title{
Estimation of the baseline survival
}
\description{
This R-program uses the Breslow method to estimate baseline survival of PH mixture cure model and AFT mixture cure model.
}
\usage{
smsurv(Time, Status, X, beta, w, model)
}
\arguments{
  \item{Time}{
this is the follow up time for "ph" model. If model is "aft", then this is residual.
}
  \item{Status}{
The status indicator, normally 0=alive, 1=dead
}
  \item{X}{
effects of covariates of uncured patients
}
  \item{beta}{
initial beta from coxph
}
  \item{w}{
conditional probability of the ith individual remains uncured at the mth iteration. We use Status as initial value
}
  \item{model}{
specifies your model, it can be  "ph" or "aft"
}
}

