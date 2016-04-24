script.start = start.time()
setwd("/Users/SeanOlson/Desktop/bayes_project/")
waic <- function(stanfit){
   log_lik <- extract(stanfit, "log_lik")$log_lik
   dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
      c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
   S <- nrow(log_lik)
   n <- ncol(log_lik)
   lpd <- log(colMeans(exp(log_lik)))
   p_waic <- colVars(log_lik)
   elpd_waic <- lpd - p_waic
   waic <- -2*elpd_waic
   loo_weights_raw <- 1/exp(log_lik-max(log_lik))
   loo_weights_normalized <- loo_weights_raw/
      matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
   loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))
   elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/
                      colMeans(loo_weights_regularized))
   p_loo <- lpd - elpd_loo
   pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
   total <- colSums(pointwise)
   se <- sqrt(n*colVars(pointwise))
   return(list(waic=total["waic"], elpd_waic=total["elpd_waic"],
               p_waic=total["p_waic"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
               pointwise=pointwise, total=total, se=se))
}

colVars <- function(a){
   n <-dim(a)[[1]];
   c <- dim(a)[[2]];
   return(.colMeans(((a - matrix(.colMeans(a, n, c), nrow = n, ncol = c, byrow = TRUE)) ^ 2), n, c) * n / (n - 1))
}



model <- "
// negative binomial parameterized as eta (log(mu)) and dispersion (phi)
// see p286 in stan-reference-2.4.0.pdf
data {
   int<lower=1> n; // number of observations
   int<lower=1> p; // number of cols in model matrix
   int <lower=1> nBlocks; // number of blocks
   int y[n]; // the response
   matrix[n,p] X; // the model matrix
   int<lower=1,upper=nBlocks> blk[n];
}
parameters {
   vector[p] beta; // treatment effects
   real<lower=0.1> phi; // dispersion parameter
   vector[nBlocks] b;
   real<lower=0> tau_b;
   //vector[n] mu;
}
transformed parameters{
   real sig2b;
   //vector[n] mu;
   sig2b <- 1/tau_b;
   //mu<-exp(X*beta+b[blk]*sqrt(sig2b));
}
model {
   vector[n] mu; // linear predictor
   real sigb;
   //vector[n] log_lik;
   phi ~ gamma(0.5,0.5);
   tau_b ~ gamma(0.001,0.001);
   beta ~ normal(0,10);
   sigb<-1/sqrt(tau_b);
   b~normal(0,0.5);
   mu<-exp(X*beta+b[blk]*sigb);
   y ~ neg_binomial_2(mu, phi);
   //log_lik <- neg_binomial_2_log(y, mu, phi);

}

generated quantities{
   real log_lik;
   log_lik <- neg_binomial_2_log(y, exp(X*beta+b[blk]*sqrt(sig2b)), phi);
}
"
write(model, "nbmodel.stan")
require(rstan)

x <- c(1,2)
N <- length(x)

fits <- list()
mod <- stan_model("gp-sim.stan")
for(i in 1:100)
{
   fits[i] <- sampling(mod, data=list(x=x,N=N), iter=1, chains=1)
}
mod = stan_model("nbmodel.stan")
nb.fit = list()
for(i in 1:1000){
   cat(paste0("\n RUN ", i, "\n\n"))
   nb.data = nb.data1[which(nb.data1$expt==i),]
   n=dim(nb.data)[1];
   X=model.matrix(~factor(nb.data$treatment)-1,data=nb.data);
   Z = model.matrix(~factor(nb.data$Block)-1, data=nb.data);
   p=dim(X)[2]


   nb.dat=list(n=n,
               p=p,
               nBlocks=max(nb.data$Block),
               y=nb.data$count,
               X=X,
               blk=nb.data$Block);
   nb.fit[[i]] = sampling(mod, data = nb.dat, iter = 50000, chains = 4)
}
nb.run = function(i, ...){
   cat(paste0("\n RUN ", i, "\n\n"))
   nb.data = nb.data1[which(nb.data1$expt==i),]
   n=dim(nb.data)[1];
   X=model.matrix(~factor(nb.data$treatment)-1,data=nb.data);
   Z = model.matrix(~factor(nb.data$Block)-1, data=nb.data);
   p=dim(X)[2]


   nb.dat=list(n=n,
               p=p,
               nBlocks=max(nb.data$Block),
               y=nb.data$count,
               X=X,
               blk=nb.data$Block);
    stan("nbmodel.stan", data =nb.dat,
                       iter = 50000, chains = 4)
}

start = proc.time()
library(rstan)
#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())
nb.data1 = read.csv("nega.csv")
nb.data1 = nb.data1[,-16]

nb.fit = list()
for(i in 1:50){
   nb.fit[[i]] <- nb.run(i)
}
for(i in 51:100){
   nb.fit[[i]] <- nb.run(i)
}
for(i in 101:150){
   nb.fit[[i]] <- nb.run(i)
}
for(i in 151:200){
   nb.fit[[i]] <- nb.run(i)
}

for(i in 201:250){
   nb.fit[[i]] <- nb.run(i)
}

for(i in 251:300){
   nb.fit[[i]] <- nb.run(i)
}

for(i in 301:350){
   nb.fit[[i]] <- nb.run(i)
}

for(i in 351:400){
   nb.fit[[i]] <- nb.run(i)
}


for(i in 401:450){
   nb.fit[[i]] <- nb.run(i)
}

for(i in 451:500){
   nb.fit[[i]] <- nb.run(i)
}


for(i in 501:550){
   nb.fit[[i]] <- nb.run(i)
}

for(i in 551:600){
   nb.fit[[i]] <- nb.run(i)
}


for(i in 601:650){
   nb.fit[[i]] <- nb.run(i)
}

for(i in 651:700){
   nb.fit[[i]] <- nb.run(i)
}


for(i in 701:750){
   nb.fit[[i]] <- nb.run(i)
}

for(i in 751:800){
   nb.fit[[i]] <- nb.run(i)
}


for(i in 801:850){
   nb.fit[[i]] <- nb.run(i)
}

for(i in 851:900){
   nb.fit[[i]] <- nb.run(i)
}


for(i in 901:950){
   nb.fit[[i]] <- nb.run(i)
}

for(i in 951:1000){
   nb.fit[[i]] <- nb.run(i)
}



for(i in 1:1000){
   cat(paste0("\n RUN ", i, "\n\n"))
   nb.data = nb.data1[which(nb.data1$expt==i),]
   n=dim(nb.data)[1];
   X=model.matrix(~factor(nb.data$treatment)-1,data=nb.data);
   Z = model.matrix(~factor(nb.data$Block)-1, data=nb.data);
   p=dim(X)[2]


   nb.dat=list(n=n,
               p=p,
               nBlocks=max(nb.data$Block),
               y=nb.data$count,
               X=X,
               blk=nb.data$Block);
   nb.fit[[i]] <- stan("nbmodel.stan", data =nb.dat,
               iter = 50000, chains = 4)
}

end = proc.time()
(end - start) / 60


######################################


model <- "
// negative binomial parameterized as eta (log(mu)) and dispersion (phi)
// see p286 in stan-reference-2.4.0.pdf
data {
int<lower=1> n; // number of observations
int<lower=1> p; // number of cols in model matrix
int <lower=1> nBlocks; // number of blocks
int <lower=1> nUnits; // no. units w/i blocks
int y[n]; // the response
  matrix[n,p] X; // the model matrix
   int<lower=1,upper=nBlocks> blk[n];
   int<lower=1, upper=p> unit[n];
}
parameters {
vector[p] beta; // treatment effects
real<lower=0.1> phi; // dispersion parameter
vector[nBlocks] b;
real<lower=0> tau_b;
vector[n] u;
real<lower=0> tau_u;
//vector[n] mu;
}
transformed parameters{
real sig2b;
real sig2u;
sig2b <- 1/tau_b;
sig2u <- 1/tau_u;
//mu<-exp(X*beta+b[blk]*sqrt(sig2b));
}
model {
vector[n] mu; // linear predictor
real sigb;
real sigu;
phi ~ gamma(0.5,0.5);
tau_b ~ gamma(0.001,0.001);
tau_u ~ gamma(0.001, 0.001);
beta ~ normal(0,10);
sigb<-1/sqrt(tau_b);
sigu <- 1/sqrt(tau_u);
b~normal(0,0.5);
u ~ normal(0,1);
mu<-exp(X*beta+b[blk]*sigb + u[unit]*sigu);
y ~ poisson(mu);

}

generated quantities{
real log_lik;
log_lik <- poisson_log(y, exp(X*beta+b[blk]*sqrt(sig2b) + u[unit]*sqrt(sig2u)));
}
"
write(model, "pmodel.stan")




start = proc.time()
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
p.data1 = read.csv("poiunit.csv")

p.fit = list()
for(i in 1:500){
   cat(paste0("\n RUN ", i, "\n\n"))
   p.data = p.data1[which(p.data1$expt==i),]
   n=dim(p.data)[1];
   X=model.matrix(~factor(p.data$treatment)-1,data=p.data);
   Z = model.matrix(~factor(p.data$Block)-1, data=p.data);
   p=dim(X)[2]


   p.dat=list(n=n,
               p=p,
               nBlocks=max(p.data$Block),
               y=p.data$count,
               X=X,
               blk=p.data$Block,
              nUnits=4,
              unit = p.data$treatment);
   p.fit[[i]] <- stan("pmodel.stan", data =p.dat,
                       iter = 50000, chains = 4)
}

end = proc.time()
(end - start) / 60

script.end = proc.time()

(script.start - script.end)/60

