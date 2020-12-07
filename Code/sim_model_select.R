# this script fit model on simulated data set
# simulated data set is from sim_data.R
# data is saved in Data/sim_data.RData

load("Data/sim_data.RData")

library(ncvreg)
library(grpreg)
library(hdrm)
library(dplyr)
library(reshape2)

set.seed(12345)

# sim_list <- sim_data1[[1]]
# c_grp <-  rnorm(10, 15, runif(1, 0.25, 5)
# grp <- grp_ind
# rs <- gamma

##### function to fit model on each data set #####

sim_fit <- function(sim_list, c_grp, grp, rs, family = "gaussian"){
  
  # set up
  X_train <- sim_list$train$X
  X_test <- sim_list$test$X
  y_train <- sim_list$train$y
  y_test <- sim_list$test$y
  
  # model 1: lasso
  fit1 <- ncvreg(X_train, y_train, penalty = "lasso")
  l1 <- fit1$lambda[which.min(BIC(fit1))]
  pred1 <- predict(fit1, X = X_test, lambda = l1)                     
  pe1 <- mean((pred1 - y_test)^2)        
  grp1 <- unique(grp[coef(fit1, lambda = l1)[-1] != 0])
  cpp1 <- sum(c_grp[grp1])
  
  # model 1.5: cost-ranked lasso
  fit1.5 <- ncvreg(X_train, y_train, penalty = "lasso", penalty.factor = c_grp[grp])
  l1.5 <- fit1.5$lambda[which.min(BIC(fit1.5))]
  pred1.5 <- predict(fit1.5, X = X_test, lambda = l1.5)                     
  pe1.5 <- mean((pred1.5 - y_test)^2)        
  grp1.5 <- unique(grp[coef(fit1.5, lambda = l1.5)[-1] != 0])
  cpp1.5 <- sum(c_grp[grp1.5])
  
  # model 2: group lasso (GL)
  fit2 <- grpreg(X_train, y_train, group = grp, penalty = "grLasso", family = family)
  l2 <- fit2$lambda[which.min(BIC(fit2))]
  pred2 <- predict(fit2, X = X_test, lambda = l2)                     
  pe2 <- mean((pred2 - y_test)^2)  
  grp2 <- as.numeric(predict(fit2, type = "groups", lambda = l2))
  cpp2 <- sum(c_grp[grp2])
  
  # model 3: group lasso weighted by cost (CGL)
  fit3 <- grpreg(X_train, y_train, group = grp, penalty = "grLasso",
                           family = family, group.multiplier = c_grp)
  l3 <- fit3$lambda[which.min(BIC(fit3))]
  pred3 <- predict(fit3, X = X_test, lambda = l3)                     
  pe3 <- mean((pred3 - y_test)^2)                     
  grp3 <- as.numeric(predict(fit3, type = "groups", lambda = l3))
  cpp3 <- sum(c_grp[grp3])
  
  # model 3: group lasso weighted by cost^r (rCGL)
  fit4_list <- list()
  fit4_l <- numeric(length(rs))
  fit4_bic <- numeric(length(rs))
  
  ## choose the best l for each r
  for(j in 1:length(rs)) {
    fit4_list[[j]] <- grpreg(X_train, y_train, group = grp, penalty = "grLasso",
                             family = family, group.multiplier = c_grp^rs[j])
    fit4_l[[j]] <- fit4_list[[j]]$lambda[which.min(BIC(fit4_list[[j]]))]
    fit4_bic[j] <- min(BIC(fit4_list[[j]]))
  }
  ## Pick best r
  fit4 <- fit4_list[[which.min(fit4_bic)]]
  l4 <- fit4_l[which.min(fit4_bic)]
  pred4 <- predict(fit4, X = X_test, lambda = l4)                     
  pe4 <- mean((pred4 - y_test)^2)                     
  grp4 <- as.numeric(predict(fit4, type = "groups", lambda = l4))
  cpp4 <- sum(c_grp[grp4])
  
  pe <- c("Lasso" = pe1, "CRL" = pe1.5, "GL" = pe2, "CGL" = pe3, "rCGL" = pe4)
  cpp <- c("Lasso" = cpp1,"CRL" = cpp1.5,  "GL" = cpp2, "CGL" = cpp3, "rCGL" = cpp4)
  
  return(list(pe = pe, cpp = cpp))
  
}

##### test #####
grp_ind <- rep(c(1:10), each = 5)

gamma <- c(0, .25, .5, .75, 1)

# case 1: no cor
sim_pe1 <- list()
sim_cpp1 <- list()

for(l in 1:length(sim_data1)){
  this_fit <- sim_fit(sim_list = sim_data1[[l]], grp = grp_ind,
                      c_grp = rnorm(10, 15, runif(1, 0.25, 5)), rs = gamma)
  sim_pe1[[l]] <- this_fit$pe 
  sim_cpp1[[l]] <- this_fit$cpp
}

sim_pe1 %>% bind_rows()
sim_cpp1 %>% bind_rows()

# case 2: within cor
sim_pe2 <- list()
sim_cpp2 <- list()

for(l in 1:length(sim_data2)){
  this_fit <- sim_fit(sim_list = sim_data2[[l]], grp = grp_ind, 
                      c_grp = rnorm(10, 15, runif(1, 0.25, 5)), rs = gamma)
  sim_pe2[[l]] <- this_fit$pe 
  sim_cpp2[[l]] <- this_fit$cpp
}


# case 3: both cor
sim_pe3 <- list()
sim_cpp3 <- list()

for(l in 1:length(sim_data3)){
  this_fit <- sim_fit(sim_list = sim_data3[[l]], grp = grp_ind, 
                      c_grp = rnorm(10, 15, runif(1, 0.25, 5)), rs = gamma)
  sim_pe3[[l]] <- this_fit$pe 
  sim_cpp3[[l]] <- this_fit$cpp
}

save(sim_pe1, sim_pe2, sim_pe3, sim_cpp1, sim_cpp2, sim_cpp3,
     file = "Data/sim_model_select.RData")

