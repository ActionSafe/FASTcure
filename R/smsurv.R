haz = function(Time, Status, death_point, w, coxexp, model){

  event <- numeric(length(death_point))
  s_time = subset(Time, Status == 1)
  wexpb = w * coxexp
  orders = order(Time)
  wexpb = wexpb[orders]
  Time = Time[orders]
  event = table(sort(s_time))
  lambda = haz_(Time, event, death_point, wexpb)
  return(lambda)
}

smsurv <- function(Time, Status, X, beta, w, model) {
  death_point <- sort(unique(subset(Time, Status == 1)))

  if (model == 'ph') {
    coxexp <- exp(as.matrix(X[, -1]) %*% beta)
  }

  lambda = haz(Time, Status, death_point, w, coxexp, model)
  HHazard = cumhaz(Time, lambda, death_point)

  survival = exp(-HHazard)
  return(list(survival = survival))
}
