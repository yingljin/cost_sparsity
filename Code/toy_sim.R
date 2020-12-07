
##### Generate data #####

true_beta <- c(1, 2, 3, 4, 5,
               0, 0, 0, 0, 0,
               1, 2, 3, 4, 5,
               rep(0, 35))

grp <- rep(c(1:10), each = 5)
table(grp)

# generate grouped covariates
x_grp_all <- genDataGrp(n = 400, J = 10, K = 5, beta = true_beta,
                    family = "gaussian", SNR = 1, 
                    rho = runif(1))

# training set
x_grp <- list(X = x_grp_all$X[1:200,], y = x_grp_all$y[1:200])

# testing set
x_grp_test <- list(X = x_grp_all$X[201:400, ], y = x_grp_all$y[201:400])

# ssy <- sum((x_grp_test$y-mean(x_grp_test$y))^2)
# generate group-wise cost
c_grp <- rnorm(10, 15, runif(1, 0.25, 5))
w_grp <- c_grp/sum(c_grp)
sum(w_grp)

##### Fit regular LASSO #####
# without accounting for different cost
# why use BIC for regularization parameter

lasso.fit <- cv.glmnet(x_grp$X, x_grp$y)
# plot(lasso.fit)
# coef(lasso.fit)

lasso_grp <- lapply(list("G01", "G02", "G03", "G04", "G05",
                         "G06", "G07", "G08", "G09", "G10" ), grep, 
                    x =  colnames(x_grp$X)[which(coef(lasso.fit)[-1]!=0)])
lasso_grp <- lapply(lasso_grp, function(x)sum(x)!=0)
lasso_grp <- unlist(lasso_grp)

cpp1 <- sum(c_grp[lasso_grp])
grp1 <- paste(c(1:10)[lasso_grp], collapse = ",")

pe1 <- mean((predict(lasso.fit, newx = x_grp_test$X) - x_grp_test$y)^2)
#lasso.r2 <- 1-sum((predict(lasso.fit, newx = x_grp_test$X) - x_grp_test$y)^2)/ssy

##### Regular group lasso #####

# without considering cost
# now all the coefficients are the same magnitude
# want to know what if high cost are slightly more correlated


grp.lasso.fit <- cv.grpreg(x_grp$X, x_grp$y, group = grp, penalty = "grLasso",
                           family = "gaussian")
# plot(grp.lasso.fit)
# coef(grp.lasso.fit)

pe2 <- mean((predict(grp.lasso.fit, X = x_grp_test$X, type = "response") - x_grp_test$y)^2)
cpp2 <- sum(c_grp[as.numeric(predict(grp.lasso.fit, type = "groups"))])
grp2 <- paste(predict(grp.lasso.fit, type = "groups"), collapse = ",")

#grp.lasso.r2 <- 1-sum((predict(grp.lasso.fit, X = x_grp_test$X, type = "response") - x_grp_test$y)^2)/ssy
##### group lasso considering cost#####

# now all the coefficients are the same magnitude
# want to know what if high cost are slightly more correlated

cgrp.lasso.fit <- cv.grpreg(x_grp$X, x_grp$y, group = grp, penalty = "grLasso",
                           family = "gaussian", group.multiplier  = w_grp)
# plot(cgrp.lasso.fit)
# coef(cgrp.lasso.fit)
pe3 <- mean((predict(cgrp.lasso.fit, X = x_grp_test$X, type = "response") - x_grp_test$y)^2)

cpp3 <- sum(c_grp[as.numeric(predict(cgrp.lasso.fit, type = "groups"))])
grp3 <- paste(predict(cgrp.lasso.fit, type = "groups"), collapse = ",")

#cgrp.lasso.r2 <- 1-sum((predict(cgrp.lasso.fit, X = x_grp_test$X, type = "response") - x_grp_test$y)^2)/ssy

##### Outputs #####

model_list <- list(lasso.fit, grp.lasso.fit, cgrp.lasso.fit)
pred_error <- c("lasso" = pe1,
                "group lasso" = pe2,
                "group lasso with cost" = pe3)

cpp <- c("lasso" = cpp1,
         "group lasso" = cpp2,
         "group lasso with cost" = cpp3)

# pe <- c("lasso" = lasso.pe,
#          "group lasso" = grp.lasso.pe,
#          "group lasso with cost" = cgrp.lasso.pe)

grp <- c("lasso" = grp1,
         "group lasso" = grp2,
         "group lasso with cost" = grp3)
