#' entry point of fastcure package
#'
#' @param formula a formula object
#' @param cureform specifies the variables in the incidence
#' @param offset variable(s) with coefficient 1 in PH model or AFT model
#' @param data a data.frame in which to interpret the variables named in the formula and cureform
#' @param na.action a missing-data filter function. By default na.action = na.omit
#' @param model specifies your model ph or aft
#' @param link incidence part
#' @param Var By default Var = TRUE
#' @param emmax maximum iteration number
#' @param eps convergence criterion
#' @param nboot number of bootstrap sampling
#'
#' @return a smcure object
#' @importFrom survival coxph survreg
#' @importFrom stats optim pnorm model.extract model.frame model.matrix na.omit optim var
#' @importFrom graphics matplot lines
#' @export
#'
#' @examples
#' data(e1684)
#' pd <- smcure(Surv(FAILTIME,FAILCENS)~TRT+SEX+AGE,
#' cureform=~TRT+SEX+AGE,data=e1684,model="ph",
#' Var = FALSE)
#' printsmcure(pd,Var = FALSE)

smcure <- function(formula,cureform,offset=NULL,data,na.action=na.omit,
                   model= c("aft", "ph"),link="logit", Var=TRUE,emmax=100,
                   eps=1e-7,nboot=100, mc.cores = 1, silence = TRUE, robust = TRUE)
{
  call <- match.call()
  model <- match.arg(model)
  if(!silence)message(paste("FASTcure: program is running with parameters: \n",
                            "bootstrap =",Var, "\n", "convergence =",eps, "\n",
                            "max_iter =",emmax, "\n", "cluster_size =",mc.cores, "\n"))
  ## prepare data
  data <- na.action(data)
  n <- dim(data)[1]
  mf <- model.frame(formula,data)
  # cvars <- all.vars(cureform)
  # Z <- as.matrix(cbind(rep(1,n),data[,cvars]))
  # colnames(Z) <- c("(Intercept)",cvars)
  Z <- model.matrix(cureform, data)

  if(!is.null(offset)) {
    offsetvar <- all.vars(offset)
    offsetvar<-data[,offsetvar]}
  else offsetvar <- NULL
  Y <- model.extract(mf,"response")
  X <- model.matrix(attr(mf,"terms"), mf)
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  Time <- Y[,1]
  Status <- Y[,2]
  bnm <- colnames(Z)
  nb <- ncol(Z)
  if(model == "ph") {
    betanm <- colnames(X)[-1]
    nbeta <- ncol(X)-1}
  if(model == "aft"){
    betanm <- colnames(X)
    nbeta <- ncol(X)}
  ## initial value
  w <- Status
  b <- eval(parse(text = paste("glm", "(", "w~Z[,-1]",",family = quasibinomial(link='", link, "'",")",")",sep = "")))$coef
  if(model=="ph") beta <- coxph(Surv(Time, Status)~X[,-1]+offset(log(w)), subset=w!=0, method="breslow")$coef
  if(model=="aft") beta <- survreg(Surv(Time,Status)~X[,-1])$coef
  ## do EM algo
  emfit <- em(Time,Status,X,Z,offsetvar,b,beta,model,link,emmax,eps)
  if(emfit$tau > eps) {
    warning(paste("\nFASTcure: fail to converge after", emmax, "Iterations !"))
  } else {
    message("\nFASTcure: EM algorithm converges successfully.")
  }
  b <- emfit$b
  beta <- emfit$latencyfit
  s <- emfit$Survival
  logistfit <- emfit$logistfit
  converge_flag = TRUE
  results = NULL
  if(Var){
    if(!silence)message(paste("FASTcure: Starting bootstrap jobs with ", mc.cores, "cpu cores..."))
    if(model=="ph") {b_boot<-matrix(rep(0,nboot*nb), nrow=nboot)
    beta_boot<-matrix(rep(0,nboot*(nbeta)), nrow=nboot)
    iter <- matrix(rep(0,nboot),ncol=1)}

    if(model=="aft") {b_boot<-matrix(rep(0,nboot*nb), nrow=nboot)
    beta_boot<-matrix(rep(0,nboot*(nbeta)), nrow=nboot)}
    tempdata <- cbind(Time,Status,X,Z)
    data1<-subset(tempdata,Status==1);data0<-subset(tempdata,Status==0)
    n1<-nrow(data1);n0<-nrow(data0)

    # enable parallel computing
    cl <- makeCluster(mc.cores)
    registerDoParallel(cl)
    suppressWarnings({
      results <- foreach(i = 1:nboot, .packages = c("FASTcure"),
                         .export = c("n1", "n0", "data1", "data0", "bnm", "model",
                                     "betanm", "offsetvar", "b", "beta", "link",
                                     "emmax", "eps", "boot_step", "em", "n")
      ) %dopar% {
        boot_step(n, n1, n0, data1, data0, bnm, model, betanm, offsetvar, b, beta, link, emmax, eps)
      }
    })
    stopCluster(cl)
    if(!silence) message(paste("FASTcure: bootstrap jobs finished, all worker stopped"))
    if(robust) {
      results = Filter(function(x) (x$convergence < eps), results)
      message(paste("FASTcure: Bootstrap convergence rate: ", length(results), "/", nboot))
      if(length(results) < 2) {
        warning("FASTcure: Bootstrap jobs fail to converge !")
        converge_flag = FALSE
      } else {
        converge_flag = TRUE
      }
    }
  }
  fit<-list()
  class(fit) <- c("smcure")
  fit$logistfit <- logistfit
  fit$b <- b
  fit$beta <- beta
  if(Var & converge_flag){
    for (i in 1:length(results)) {
      b_boot[i, ] <- results[[i]]$b
      beta_boot[i, ] <- results[[i]]$latencyfit
    }
    b_var <- apply(b_boot, 2, var)
    beta_var <- apply(beta_boot, 2, var)
    b_sd <- sqrt(b_var)
    beta_sd <- sqrt(beta_var)

    fit$b_var <- b_var
    fit$b_sd <- b_sd
    fit$b_zvalue <- fit$b/b_sd
    fit$b_pvalue <- (1-pnorm(abs(fit$b_zvalue)))*2
    fit$beta_var <- beta_var
    fit$beta_sd <- beta_sd
    fit$beta_zvalue <- fit$beta/beta_sd
    fit$beta_pvalue <- (1-pnorm(abs(fit$beta_zvalue)))*2	}
  fit$call <- call
  fit$bnm <- bnm
  fit$betanm <- betanm
  fit$s <- s
  fit$Time <- Time
  fit$convergence = emfit$tau
  if(model=="aft"){
    error <- drop(log(Time)-beta%*%t(X))
    fit$error <- error}
  # printsmcure(fit,Var)
  return(fit)
}


boot_step <- function(n, n1, n0, data1, data0, bnm, model, betanm, offsetvar, b, beta, link, emmax, eps) {
  id1 <- sample(1:n1, n1, replace = TRUE)
  id0 <- sample(1:n0, n0, replace = TRUE)
  bootdata <- rbind(data1[id1, ], data0[id0, ])
  bootZ <- bootdata[, bnm]

  if (model == "ph") {
    bootX <- as.matrix(cbind(rep(1, n), bootdata[, betanm]))
  }

  if (model == "aft") {
    bootX <- bootdata[, betanm]
  }

  bootfit <- em(bootdata[, 1], bootdata[, 2], bootX, bootZ, offsetvar, b, beta, model, link, emmax, eps)
  # convergence checking
  return(list(b = bootfit$b, latencyfit = bootfit$latencyfit, convergence = bootfit$tau))
}


