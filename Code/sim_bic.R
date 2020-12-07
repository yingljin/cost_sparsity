# this script generates simulations 
# and select models with BIC

##### generate data #####

# true covariate
true_beta <- c(1, 2, 3, 4, 5,
               0, 0, 0, 0, 0,
               1, 2, 3, 4, 5,
               rep(0, 35))

# group index
grp <- rep(c(1:10), each = 5)
table(grp)

# all data
x_grp_all <- genDataGrp(n = 400, J = 10, K = 5, beta = true_beta,
                        family = "gaussian", SNR = 1, 
                        rho = 0.5, rho.g = 0.5)

table(x_grp_all$group)
x_grp_all$beta
cor(x_grp_all$X)

# training set
x_grp <- list(X = x_grp_all$X[1:200,], y = x_grp_all$y[1:200])

# testing set
x_grp_test <- list(X = x_grp_all$X[201:400, ], y = x_grp_all$y[201:400])

# generate group-wise cost
c_grp <- rnorm(10, 15, runif(1, 0.25, 5))
w_grp <- c_grp/sum(c_grp)
sum(w_grp)

##### Lasso #####
# penalize each covariate equally
lasso.fit <- ncvreg(x_grp$X, x_grp$y, penalty = "lasso")

# select best model by BIC
l1 <- lasso.fit$lambda[which.min(BIC(lasso.fit))]

# prediction error
pred1 <- predict(lasso.fit, X = x_grp_test$X, lambda = l1)                     
pe1 <- mean((pred1 - x_grp_test$y)^2)                     

# cost 
cpp1 <- c_grp[unique(grp[coef(lasso.fit, lambda = l1)[-1] != 0])]
cpp1 <- sum(cpp1)

##### Group lasso #####
# penaliza each group equally
grp.lasso.fit <- grpreg(x_grp$X, x_grp$y, group = grp, penalty = "grLasso",
                           family = "gaussian")


# select best model by BIC
l2 <- grp.lasso.fit$lambda[which.min(BIC(grp.lasso.fit))]

# prediction error
pred2 <- predict(grp.lasso.fit, X = x_grp_test$X, lambda = l2)                     
pe2 <- mean((pred2 - x_grp_test$y)^2)                     

# cost
grp2 <- as.numeric(predict(grp.lasso.fit, type = "groups", lambda = l2))
cpp2 <- sum(c_grp[grp2])

##### Group lasso with cost #####
cgrp.lasso.fit <- grpreg(x_grp$X, x_grp$y, group = grp, penalty = "grLasso",
                        family = "gaussian", group.multiplier = w_grp)

# select best model by BIC
l3 <- cgrp.lasso.fit$lambda[which.min(BIC(cgrp.lasso.fit))]

# prediction error
pred3 <- predict(cgrp.lasso.fit, X = x_grp_test$X, lambda = l3)                     
pe3 <- mean((pred3 - x_grp_test$y)^2)                     

# cost
grp3 <- as.numeric(predict(cgrp.lasso.fit, type = "groups", lambda = l3))
cpp3 <- sum(c_grp[grp3])

##### Group lasso with cost tunned by gamma no scale #####

# gamma
rs <- c(0, .25, .5, .75, 1)

# select best gamma and lambda
cgrp_fits <- list()
cgrp_ls <- numeric(length(rs))
cgrp_bics <- numeric(length(rs))

for(i in 1:length(rs)) {
  cgrp_fits[[i]] <- grpreg(x_grp$X, x_grp$y, group = grp, penalty = "grLasso",
                           family = "gaussian", group.multiplier = c_grp^rs[i])
  cgrp_ls[[i]] <- cgrp_fits[[i]]$lambda[which.min(BIC(cgrp_fits[[i]]))]
  cgrp_bics[i] <- min(BIC(cgrp_fits[[i]]))
}

# Pick best model
best_cgrp <- cgrp_fits[[which.min(cgrp_bics)]]
best_l <- cgrp_ls[which.min(cgrp_bics)]

## prediction error
pred4 <- predict(best_cgrp, X = x_grp_test$X, lambda = best_l)                     
pe4 <- mean((pred4 - x_grp_test$y)^2)                     

## cost
grp4 <- as.numeric(predict(best_cgrp, type = "groups", lambda = best_l))
cpp4 <- sum(c_grp[grp4])


##### Results #####

pred_error <- c("lasso" = pe1, "grp_lasso" = pe2, "cost_grep_lasso" = pe3, "wcost_grp_lasso" = pe4)
cpp <- c("lasso" = cpp1, "grp_lasso" = cpp2, "cost_grep_lasso" = cpp3, "wcost_grp_lasso" = cpp4)
