# Homework 3
Matthew Farr  
October 27, 2016  
FW 849
Fall 2014
HW3: Model comparison and model checking

The cpus dataset in R can be found using:
>library(MASS)
>data(cpus)
It includes a relative measure of performance and characteristics of 209 cpus. For more details about this dataset, please use the help command:
>?cpus

This dataset is often used to illustrate model selection in regression. We will use this data in this HW. HINT: use the provided R code to upload and standardize the data.

A preliminary analysis has found that the variables mmin (minimum main memory), mmax (maximum main memory), cach (cache size) and chmax (maximum number of chanels) allow to build a model with good predictive properties for the performance measure variable (perf). Please answer the following questions:

<a name="8"></a>

1)	Write a linear regression model that allows to predict performance using the 4 listed predictor variables (assume indepdendently and identically distributed errors), make sure to include all prior assumptions. (20 points)

<span style="color:red">See [BUGS Model](#1) for the linear regression.

2)	Fit the model using OpenBUGS. Provide convergence diagnostics (10 points) and posterior distribution of parameters of interest: linear coefficients and error variance (10). HINT: be careful with computational cost! Start with small number of iterations (5000) and perform diagnostics only on linear parameters and precision. Once you decide the number of iterations and burn in, you can start monitoring other quantities.

<span style="color:red">See [Model 1 Test](#2) for test run posterior distributions and  convergence diagnostics.

<span style="color:red">See [Model 1](#3) for full run posterior distributions and convergence diagnostics. Note that Gelman and Rubin diagnostic is contained within the summary output under rhat.

3)	Assess the posterior predictive distribution of FUTURE observations. Hint: a graphic of predicted versus observed or residuals versus observed could be used, BUT make sure you generate FUTURE observations. (20 points)

<span style="color:red">See [Future vs Observed Residuals](#4) to assess the model fit for FUTURE observations. The plot suggest that the model has good fit, and would be able to predict FUTURE observations.

4)	An alternative model to the presented model is the model that excludes mmin and chmax (model that predicts perf using only mmax and cach). Fit this model using OpenBUGS (20 points)

<span style="color:red">The [BUGS Model](#1) was reused, but the data was [reformatted](#5) to exlude "mmin" and "chmax." See [Model 2](#6) for full run posterior distributions and convergence diagnostics. Note that Gelman and Rubin diagnostic is contained within the summary output under rhat.

5)	If you have to choose between the model in 1) and the model in 4), which one would you choose and why? HINT: if you want to re-fit the model to extract predictive criteria, monitor the minimum number of parameters to minimize memory use and output.

<span style="color:red">The [DIC comparison](#7) shows that [Model 1](#3) better fits the data than [Model 2](#6). Thus the full model is better than the model when removing mmin and chmax and those 2 predictors are important in describing the performance of CPUs.

Load Packages and Data

```r
library(MASS)
library(coda)
library(jagsUI)
```

```
## Loading required package: lattice
```

```
## 
## Attaching package: 'jagsUI'
```

```
## The following object is masked from 'package:coda':
## 
##     traceplot
```

```r
data(cpus)

y<-cpus$perf
x<-as.matrix(cpus[,c("mmin","mmax","cach","chmax")])
x<-apply(x,2,scale)
head(x)
```

```
##            mmin       mmax      cach       chmax
## [1,] -0.6734091 -0.4942755 5.6805690  4.22089912
## [2,]  1.3231141  1.7229127 0.1672280  0.52821054
## [3,]  1.3231141  1.7229127 0.1672280  0.52821054
## [4,]  1.3231141  1.7229127 0.1672280  0.52821054
## [5,]  1.3231141  0.3584892 0.1672280 -0.08723756
## [6,]  1.3231141  1.7229127 0.9548481  0.52821054
```

```r
y<-(y-mean(y))/sd(y)

N <- length(y)
p <- length(x[1,])
```

Look at standardized Data

```r
hist(y)
```

![](HW3_files/figure-html/unnamed-chunk-2-1.png)

<a name="1"></a>BUGS Model

```r
cat("
  model{
  
  #Priors
  
  tau ~ dgamma(0.0001, 0.0001)
  for(i in 1:p){
  beta[i] ~ dnorm(0, 0.001)
  }

  #Likelihood
  for(i in 1:N){
  mu[i] <- inprod(x[i,], beta)
  y[i] ~ dnorm(mu[i], tau)
  yfuture[i] ~ dnorm(mu[i], tau)
  e[i] <- y[i] - mu[i]
  enew[i] <- yfuture[i] - mu[i]
  }
  fit <- sum(e[])
  fitnew <- sum(enew[])

  }",fill=TRUE,file="hw3_1.txt")
```

[Return to Top](#8)

Compile Data, Generate Inits, and Run Test Model

```r
data.regress<-list(y = y, N = N, p = p, x = x)

inits<-function(){
  list(tau=runif(1,0,10)
  )
}

params <-c("beta","tau")

out <- jags(data = data.regress, inits = inits, parameters.to.save = params, 
            model.file = "hw3_1.txt", n.chains = 1, n.iter = 5100, n.burnin = 100, n.thin = 1)
```

```
## 
## Processing function input....... 
## 
## Done. 
##  
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 209
##    Unobserved stochastic nodes: 214
##    Total graph size: 2074
## 
## Initializing model
## 
## Adaptive phase, 100 iterations x 1 chains 
## If no progress bar appears JAGS has decided not to adapt 
##  
## 
##  Burn-in phase, 100 iterations x 1 chains 
##  
## 
## Sampling from joint posterior, 5000 iterations x 1 chains 
##  
## 
## Calculating statistics.......
```

```
## Warning in process.output(samples, DIC = DIC, codaOnly, verbose = verbose):
## At least one Rhat value could not be calculated.
```

```
## 
## Done.
```

<a name="2"></a>Model 1 Test Posterior Distribution Summary and Convergence Diagnostics

```r
out
```

```
## JAGS output for model 'hw3_1.txt', generated by jagsUI.
## Estimates based on 1 chains of 5100 iterations,
## burn-in = 100 iterations and thin rate = 1,
## yielding 5000 total samples from the joint posterior. 
## MCMC ran for 0.005 minutes at time 2016-10-27 17:31:24.
## 
##             mean    sd    2.5%     50%   97.5% overlap0 f
## beta[1]    0.358 0.044   0.273   0.357   0.446    FALSE 1
## beta[2]    0.389 0.048   0.294   0.390   0.483    FALSE 1
## beta[3]    0.148 0.034   0.080   0.148   0.214    FALSE 1
## beta[4]    0.233 0.034   0.164   0.233   0.300    FALSE 1
## tau        7.009 0.690   5.710   6.984   8.433    FALSE 1
## deviance 186.998 3.242 182.686 186.344 194.929    FALSE 1
## 
## overlap0 checks if 0 falls in the parameter's 95% credible interval.
## f is the proportion of the posterior with the same sign as the mean;
## i.e., our confidence that the parameter is positive or negative.
## 
## DIC info: (pD = var(deviance)/2) 
## pD = 5.3 and DIC = 192.254 
## DIC is an estimate of expected predictive error (lower is better).
```

```r
nm.theta <- paste("beta[",1:p,"]",sep="")
beta <- out$samples[[1]][,nm.theta]
traceplot(beta)
```

![](HW3_files/figure-html/unnamed-chunk-5-1.png)![](HW3_files/figure-html/unnamed-chunk-5-2.png)![](HW3_files/figure-html/unnamed-chunk-5-3.png)![](HW3_files/figure-html/unnamed-chunk-5-4.png)

```r
autocorr.plot(beta)
```

![](HW3_files/figure-html/unnamed-chunk-5-5.png)

```r
nm.theta <- paste("tau",sep="")
tau <- out$samples[[1]][,nm.theta]
traceplot(tau)
```

![](HW3_files/figure-html/unnamed-chunk-5-6.png)

```r
autocorr.plot(tau)
```

![](HW3_files/figure-html/unnamed-chunk-5-7.png)

[Return to Top](#8)

Set Parameters and Run Full Model 1

```r
params2 <- c("beta", "tau", "fit", "fitnew")

out2 <- jags(data = data.regress, inits = inits, parameters.to.save = params2, 
model.file = "hw3_1.txt", n.chains = 3, n.iter = 15000, n.burnin = 5000, n.thin = 8)
```

```
## 
## Processing function input....... 
## 
## Done. 
##  
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 209
##    Unobserved stochastic nodes: 214
##    Total graph size: 2074
## 
## Initializing model
## 
## Adaptive phase, 100 iterations x 3 chains 
## If no progress bar appears JAGS has decided not to adapt 
##  
## 
##  Burn-in phase, 5000 iterations x 3 chains 
##  
## 
## Sampling from joint posterior, 10000 iterations x 3 chains 
##  
## 
## Calculating statistics....... 
## 
## Done.
```

<a name="3"></a>Full Model 1 Posterior Distribution Summary and Convergence Diagnostics

```r
out2[["summary"]][1:5,]
```

```
##              mean         sd       2.5%       25%       50%       75%
## beta[1] 0.3584187 0.04343105 0.27298087 0.3288986 0.3580263 0.3880736
## beta[2] 0.3888643 0.04820845 0.29583226 0.3556259 0.3886584 0.4221010
## beta[3] 0.1486809 0.03431302 0.08235207 0.1252666 0.1488715 0.1714051
## beta[4] 0.2321486 0.03450303 0.16548390 0.2093484 0.2323988 0.2547066
## tau     7.0112845 0.69510499 5.66525556 6.5317317 7.0034394 7.4761168
##             97.5%     Rhat n.eff overlap0 f
## beta[1] 0.4433627 1.000073  3750        0 1
## beta[2] 0.4859403 1.000113  3750        0 1
## beta[3] 0.2163549 1.000498  3750        0 1
## beta[4] 0.2999942 1.000203  3750        0 1
## tau     8.3766959 1.000272  2985        0 1
```

```r
traceplot(out2, "beta")
```

![](HW3_files/figure-html/unnamed-chunk-7-1.png)![](HW3_files/figure-html/unnamed-chunk-7-2.png)![](HW3_files/figure-html/unnamed-chunk-7-3.png)![](HW3_files/figure-html/unnamed-chunk-7-4.png)

[Return to Top](#8)

<a name="4"></a>Future Residuals vs Observed Residuals

```r
plot(out2$sims.list$fit, out2$sims.list$fitnew)
abline(0,1, lwd = 2, col = "black")
```

![](HW3_files/figure-html/unnamed-chunk-8-1.png)

[Return to Top](#8)

<a name="5"></a>Model 2 Reformat and Recompile Data

```r
x1<-as.matrix(cpus[,c("mmax","cach")])
x1<-apply(x1,2,scale)
head(x1)
```

```
##            mmax      cach
## [1,] -0.4942755 5.6805690
## [2,]  1.7229127 0.1672280
## [3,]  1.7229127 0.1672280
## [4,]  1.7229127 0.1672280
## [5,]  0.3584892 0.1672280
## [6,]  1.7229127 0.9548481
```

```r
p1 <- length(x1[1,])

data.regress1<-list(y = y, N = N, p = p1, x = x1)
```

[Return to Top](#8)

Run Full Model 2

```r
out3 <- jags(data = data.regress1, inits = inits, parameters.to.save = params, 
model.file = "hw3_1.txt", n.chains = 3, n.iter = 15000, n.burnin = 5000, n.thin = 8)
```

```
## 
## Processing function input....... 
## 
## Done. 
##  
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 209
##    Unobserved stochastic nodes: 212
##    Total graph size: 1538
## 
## Initializing model
## 
## Adaptive phase, 100 iterations x 3 chains 
## If no progress bar appears JAGS has decided not to adapt 
##  
## 
##  Burn-in phase, 5000 iterations x 3 chains 
##  
## 
## Sampling from joint posterior, 10000 iterations x 3 chains 
##  
## 
## Calculating statistics....... 
## 
## Done.
```

[Return to Top](#8)

<a name="6"></a>Full Model 2 Posterior Distribution Summary and Convergence Diagnostics

```r
out3[["summary"]][1:3,]
```

```
##              mean         sd      2.5%       25%       50%       75%
## beta[1] 0.7127143 0.03724684 0.6404161 0.6876106 0.7123626 0.7379061
## beta[2] 0.2793986 0.03680358 0.2072741 0.2540355 0.2790458 0.3048541
## tau     4.9754032 0.48292019 4.0883494 4.6426904 4.9565175 5.2882545
##             97.5%      Rhat n.eff overlap0 f
## beta[1] 0.7862381 1.0001106  3750        0 1
## beta[2] 0.3519813 0.9999853  3750        0 1
## tau     5.9712714 1.0000377  3750        0 1
```

```r
traceplot(out3, "beta")
```

![](HW3_files/figure-html/unnamed-chunk-11-1.png)![](HW3_files/figure-html/unnamed-chunk-11-2.png)

```r
traceplot(out3, "tau")
```

![](HW3_files/figure-html/unnamed-chunk-11-3.png)

[Return to Top](#8)

<a name="7"></a>Compare DIC between Model 1 and Model 2

Model1

```r
out2$DIC
```

```
## [1] 192.1735
```
Model2

```r
out3$DIC
```

```
## [1] 261.8292
```

[Return to Top](#8)
