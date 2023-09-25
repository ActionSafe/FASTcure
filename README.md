# FASTcure <!-- omit in toc -->
An R-package for Estimating Semiparametric PH Cure Models. Optimized for big data.

**Table of Contents**
- [1. Motivation](#1-motivation)
- [2. Modification](#2-modification)
  - [2.1. Efficiency](#21-efficiency)
  - [2.2. Bug fix](#22-bug-fix)
    - [2.2.1. Problem with factor variable in incidence model](#221-problem-with-factor-variable-in-incidence-model)
    - [2.2.2. Problem with one covariate](#222-problem-with-one-covariate)
    - [2.2.3. Infinite bootstrap loop](#223-infinite-bootstrap-loop)
- [3. Usage](#3-usage)
- [4. Simulation result](#4-simulation-result)
- [5. Todos](#5-todos)


## 1. Motivation
Existing R-packages for estimating PH Cure Models are very inefficient and quite buggy, so I modify the original package [smcure](https://cran.r-project.org/web/packages/smcure/index.html) just to make it faster and more reliable (hopefully).

## 2. Modification
### 2.1. Efficiency
As you can see, when you feed too much data to smcure (in this case, I feed a `34, 008` row dataframe to it), it will consume a significant amount of memory (up to 14 GB !). Besides, `smcure` will only utilize one cpu core when performs bootstrap steps, even if your new PC or powerful server has many cpu cores !

![smcure](imgs/smcure.png)

After profiling smcure, I found this performance bottleneck function:

```R
smsurv <-
function(Time,Status,X,beta,w,model){    
    death_point <- sort(unique(subset(Time, Status==1)))
	if(model=='ph') coxexp <- exp((beta)%*%t(X[,-1]))  
    lambda <- numeric()
    event <- numeric()
      for(i in 1: length(death_point)){
       event[i] <- sum(Status*as.numeric(Time==death_point[i]))
                 if(model=='ph')  temp <- sum(as.numeric(Time>=death_point[i])*w*drop(coxexp))
       		if(model=='aft')  temp <- sum(as.numeric(Time>=death_point[i])*w)
                  temp1 <- event[i]
       lambda[i] <- temp1/temp
        }
    HHazard <- numeric()
    for(i in 1:length(Time)){
        HHazard[i] <- sum(as.numeric(Time[i]>=death_point)*lambda)
        if(Time[i]>max(death_point))HHazard[i] <- Inf
        if(Time[i]<min(death_point))HHazard[i] <- 0
        }
   survival <- exp(-HHazard)
   list(survival=survival)
}
```

There are two main reasons why this function is slow:

- Frequent GC: this expression `lambda <- numeric()` create an vector and subsequent code dynamically grows its size, which invokes too many memory allocation system call.
- Nested loop: when your data has many rows, the two `for` loops inside this function can be quite slow.

To fix these issues, I rewrite this funciton in C++ (`haz.cpp` and `cumhaz.cpp`) using [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) library and use [foreach](https://cran.r-project.org/web/packages/foreach/index.html) package to parallize bootstrap steps.

The following screenshot shows that after optimization, you can use all of your cpu cores and there's a significant reduction in memory usage (and GC frequence, which is not shown here)

![fastcure](imgs/fastcure.png)

### 2.2. Bug fix

#### 2.2.1. Problem with factor variable in incidence model

```R
cvars <- all.vars(cureform)
Z <- as.matrix(cbind(rep(1,n),data[,cvars]))
colnames(Z) <- c("(Intercept)",cvars)
```

In this code block, `as.matrix` will convert dataframe to matrix, but if your dataframe contains factor variables, `as.matrix` will convert the whole dataframe to string matrix (because matrix can hold exactly one type of data), this behaviour will cause subsequent code crashes:

```R
b <- eval(parse(text = paste("glm", "(", "w~Z[,-1]",",family = quasibinomial(link='", link, "'",")",")",sep = "")))$coef
```

So I use `model.matrix` to extract covariates.

#### 2.2.2. Problem with one covariate

When you specify only one covariate in submodels like the following code:

```R
smcure(formula = Surv(time, cancer) ~ stage)
```
smcure will crash with error: `non-conformable arguments`, this problem is due to the following expression in `smsurv` function:

```R
exp(X[, -1] %*% beta)
```

When X has only two columns, `X[, -1]` will produce an atomic vector instead of a matrix, this behaviour will cause `non-conformable arguments` error. So a `dim` attribute should be added to this vector to make it a matrix.

#### 2.2.3. Infinite bootstrap loop

```R
i<-1
while (i<=nboot){
  # some code
  if (bootfit$tau<eps) i<-i+1
}
```

When em does not converge (this is common when the dataset is relatively large), this loop will never finish.

## 3. Usage

```R
smcure(formula = lformula,
                cureform = iformula,
                data = d, model = "ph",emmax = 100, eps = 1e-06, nboot = 50, mc.cores = 10, Var = TRUE, robust = TRUE, silence = TRUE)
```

`mc.cores = 10`: create 10 sub-processes to run bootstrap steps

`robust = TRUE`: only use converged bootstrap step to caculate std.error

`silence=TRUE`: reduce the number of messages printed to console


## 4. Simulation result

Setting: `n=300, nsims = 500`

![sims](/imgs/sims.png)

## 5. Todos

- use better std.error estimation method
- testing AFT model
- update docs
- update plot function

