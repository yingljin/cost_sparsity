# this script defines a function that fit model for 1 data set 
# and saves the fitted model
# three types of data are generated
# each type contains 1000 simulated data set

#library(ncvreg)
library(grpreg)
library(hdrm)
library(dplyr)
library(reshape2)

set.seed(12345)

load("Data/sim_data.RData")

##### function to fit model on each data set #####

sim_fit <- function(sim_list, c_grp, grp, rs, family = "gaussian"){
  
  # set up
  X_train <- sim_list$train$X
  X_test <- sim_list$test$X
  y_train <- sim_list$train$y
  y_test <- sim_list$test$y
  
  # model 1: lasso
  fit1 <- ncvreg(X_train, y_train, penalty = "lasso")
  
  # model 1.5: cost-ranked lasso
  fit1.5 <- ncvreg(X_train, y_train, penalty = "lasso", penalty.factor = c_grp[grp])
  
  # model 2: group lasso (GL)
  fit2 <- grpreg(X_train, y_train, group = grp, penalty = "grLasso", family = family)
  
  # model 3: group lasso weighted by cost (CGL)
  fit3 <- grpreg(X_train, y_train, group = grp, penalty = "grLasso",
                 family = family, group.multiplier = c_grp)
  
  # model 3: group lasso weighted by cost^r (rCGL)
  fit4_list <- list()
  
  ## choose the best l for each r
  for(j in 1:length(rs)) {
    fit4_list[[j]] <- grpreg(X_train, y_train, group = grp, penalty = "grLasso",
                             family = family, group.multiplier = c_grp^rs[j])
  }
  return(list(lasso = fit1, crl=fit1.5, gl = fit2, cgl = fit3, rcgl = fit4_list))
  
}

# generate cost
cost <- rnorm(10, 15, runif(1, 0.25, 5))

grp_ind <- rep(c(1:10), each = 5)

gamma <- c(0, .25, .5, .75, 1)

#### case 1: no cor #####

this_fit1 <- sim_fit(sim_list = sim_data1[[1]], grp = grp_ind, c_grp = cost, rs = gamma)
this_fit2 <- sim_fit(sim_list = sim_data2[[1]], grp = grp_ind, c_grp = cost, rs = gamma)
this_fit3 <- sim_fit(sim_list = sim_data3[[1]], grp = grp_ind, c_grp = cost, rs = gamma)

save(this_fit1, this_fit2, this_fit3, file = "Data/sim_model_once.RData")
