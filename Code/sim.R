## Ranked cost model

library(ncvreg)
library(mvtnorm)
library(visreg)

S <- 10000
set.seed(13)

AICc <- function(fit, eps = 1) { 
  ll <- logLik(fit)
  k <- attr(ll, "df")
  n <- attr(ll, "n")
  k_star <- min(k, n - eps - 1)
  AIC(fit) + (2 * k^2 + 2*k) / (n - k_star - 1)
}

n <- 200
p <- 20
beta <- c(1:5, rep(0, p - 5))
sigma_e <- 12
x_corr_s <- function() runif(1, 0, 1)
x_cost_var_s <- function() runif(1, 0.25,5)
pf_exp <- c(0, .25, .5, .75, 1)
block <- F

save_fname <- "M:/General Research Directory/Peterson Thesis/Ranked Skepticism/Ranked Cost/sim_results_rsq.RData"
pf_hat <- numeric(S)

rsq <- costs <- matrix(ncol = 2, nrow = S)
x_corr <- var_cost_var <- optimal_model_cost <- numeric(S)


pb <- dplyr::progress_estimated(S)
for(s in 1:S) {
  x_corr[s] <- x_corr_s()
  
  if(block) {
    Sigma11 <- pmin(diag(5*2) + x_corr[s], 1)
    Sigma12 <- matrix(0, nrow = 10, ncol = p-10)
    Sigma22 <- diag(p-5*2)
    Sigma21 <- matrix(0, nrow = p-10, ncol = 10)
    Sigma <- rbind(cbind(Sigma11, Sigma12), cbind(Sigma21, Sigma22))
    
  } else Sigma <- pmin(diag(p) + x_corr[s], 1)
  
  X <- rmvnorm(n = n, sigma = Sigma)
  var_cost_var[s] <- x_cost_var_s()
  var_cost <- pmax(0.1, rnorm(p, mean = 15, sd=  var_cost_var[s]))
  optimal_model_cost[s] <- (beta != 0) %*% var_cost
  y <- X %*% beta + rnorm(n, sd = sigma_e)
  
  fit0 <- ncvreg(X, y, penalty = "lasso", penalty.factor = rep(1, p))
  l_hat0 <- fit0$lambda[which.min(AICc(fit0))]
  s0 <- summary(fit0, lambda = l_hat0)
  
  fits1 <- s1 <- list()
  l_hats1 <- numeric(length(pf_exp))
  for(pf in 1:length(pf_exp)) {
    fits1[[pf]] <- ncvreg(X, y, penalty = "lasso", penalty.factor = (var_cost)^pf_exp[pf])
    l_hats1[[pf]] <- fits1[[pf]]$lambda[which.min(AICc(fits1[[pf]]))]
    s1[[pf]] <- summary(fits1[[pf]] , lambda = l_hats1[[pf]] )
  }
  
  ## Pick best SR model
  AICcs <- sapply(lapply(fits1, AICc), min)
  
  pf_hat[s] <- pf_exp[which.min(AICcs)]
  fit1 <- fits1[[which.min(AICcs)]]
  l_hat1 <- l_hats1[[which.min(AICcs)]]
  
  ## Gen new data
  X_new <- rmvnorm(n = 10000, sigma = pmin(diag(p) + x_corr[s], 1))
  y_new <- X_new %*% beta + rnorm(10000, 0, sigma_e)
  
  p0 <- predict(fit0, X = X_new, lambda = l_hat0)
  p1 <- predict(fit1, X = X_new, lambda = l_hat1)
  sst <- sum((y_new - mean(y_new))^2)
  rsq[s,] <- 1 - c(sum((y_new - p0)^2) / sst,  sum((y_new - p1)^2) / sst)
  
  costs[s,] <- c((coef(fit0, lambda = l_hat0) != 0)[-1] %*% var_cost,
                 (coef(fit1, lambda = l_hat1) != 0)[-1] %*% var_cost)
  pb$tick()$print()
}

boxplot(rsq)
boxplot(costs)
pf_hat

## Using paired t-test
t.test(log(rsq[,1]),log(rsq[,2]), paired = T)
t.test(log(costs[,1]), log(costs[,2]), paired = T)

## Models for R-squared
r2rat <- rsq[,2] / rsq[,1]


# Models for cost
crat <- costs[,2] / costs[,1]

par(mfrow = c(1,4), mar = c(5,4,4,2) + .1, font.main = 1)
colnames(costs) <- colnames(rsq) <- c("LS", "SRL")
boxplot(rsq, ylab = expression(R^2), main = expression(R^2))
points(1:2, apply(rsq,2, mean), pch = 4, col = 4)
boxplot(costs, ylab = "CPP", main = "CPP", font.main = 1)
points(1:2, apply(costs,2, mean), pch = 4, col = 4)

boxplot(cbind("R-squared" = r2rat, "CPP" = crat), log = "y", 
        ylab = "Ratio of SRL to LS",
        main = "Ratio of SRL to LS", font.main = 1)
points(1:2, apply(cbind("R-squared" = r2rat, "CPP" = crat),2, mean), 
       pch = 4, col = 4)

r2rat2 <- r2rat
crat2 <- crat
r2rat2[.9999 < r2rat & r2rat < 1.0001] <- NA
crat2[.9999 < crat & crat < 1.0001] <- NA

boxplot(cbind("R-squared" = r2rat2, "CPP" = crat2), log = "y", 
        ylab = "Ratio of SRL to LS(exact matches omitted)",
        main = "Ratio of SRL to LS\n(exact matches omitted)", font.main = 1)
points(1:2, apply(cbind("R-squared" = r2rat2, "CPP" = crat2),2, mean, na.rm = T), 
       pch = 4, col = 4)

## Try looking at loess plots by model, kappa, rho
get_loess <- function(x, y, xx) {
  fit_loess <- loess(y ~ x)
  pp_loess <- predict(fit_loess, newdata = xx)
  pp_loess
}

## Check out ratios
par(mfrow = c(2,2))
plot(var_cost_var, r2rat, pch = 20, main = "R-squared", ylab = "R-squared", xlab = expression(kappa), log = "y")
add_loess(var_cost_var, r2rat, col = 2, lwd = 2)
plot(x_corr, r2rat, pch = 20, main = "R-squared", ylab = "R-squared", xlab = expression(rho), log = "y")
add_loess(x_corr, r2rat, col = 2, lwd = 2)

plot(var_cost_var, crat, pch = 20, main = "Cost ratio", ylab = "Cost ratio", xlab = expression(kappa), log = "y")
add_loess(var_cost_var, crat, col = 2, lwd = 2)
plot(x_corr, crat, pch = 20, main = "Cost ratio", ylab = "Cost ratio", xlab = expression(rho), log = "y")
add_loess(x_corr, crat, col = 2, lwd = 2)

xx1 <- seq(.25, 5, length = 200)
rsq_loess <- list(
  r1 = get_loess(var_cost_var, r2rat, xx1),
  c1 = get_loess(var_cost_var, crat, xx1)
)

xx2 <- seq(0, 1, length = 200)
xcorr_loess <- list(
  r1 = get_loess(x_corr, r2rat, xx2),
  c1 = get_loess(x_corr, crat, xx2)
)

par(mfrow = c(2,2))
plot(xx1, rsq_loess$r1, main = "R-squared Ratio", 
     lwd = 2, type = "l", 
     ylab = "R-squared Ratio", xlab = expression(kappa))
plot(xx2, xcorr_loess$r1, main = "R-squared Ratio", 
     lwd = 2, type = "l", 
     ylab = "R-squared Ratio", xlab = expression(rho))

plot(xx1, rsq_loess$c1, main = "Cost Ratio", 
     ylim = range(rsq_loess$c1, 1, na.rm = T),
     lwd = 2, type = "l", 
     ylab = "Cost Ratio", xlab = expression(kappa))

plot(xx2, xcorr_loess$c1, main = "Costs", 
     ylim = range(xcorr_loess$c1, 1, na.rm = T),
     lwd = 2, type = "l", 
     ylab = "Cost/Optimal cost", xlab = expression(rho))


## GAM fits
rsq_mean <- apply(rsq, 1, mean)

bestfit_r2 <- mgcv::gam(log(r2rat) ~ s(x_corr) + s(var_cost_var))
bestfit_c <- mgcv::gam(log(crat) ~ s(x_corr) + s(var_cost_var))

summary(bestfit_r2)
summary(bestfit_c)

ylim1 <- c(.99, 1.01)
ylim2 <- c(.75, 1.05)

par(mfrow = c(2,2), mar = c(2, 5, 4, 1.2))
visreg(bestfit_r2, "var_cost_var", partial = F, rug = F, trans = exp, #ylim = ylim1, 
       ylab = expression(R[SRL]^2/R[LS]^2), xlab = "",
       main = expression(Varying~cost~variance*","~~kappa))
abline(h = 1, lty = 2)
par(mar = c(2, 3, 4, 2))
visreg(bestfit_r2, "x_corr", partial = F, rug = F, trans = exp, #ylim = ylim1, 
       ylab = "", xlab = "",
       main = expression(Varying~covariate~correlation*","~rho))
abline(h = 1, lty = 2)

par(mar = c(5, 5, 2, 1.2))
visreg(bestfit_c, "var_cost_var", partial = F, rug = F, trans = exp, #ylim = ylim2, 
       ylab = expression(CPP[SRL]/CPP[LS]), xlab = expression(kappa))
abline(h = 1, lty = 2)
par(mar = c(5, 3, 2, 2))
visreg(bestfit_c, "x_corr", partial = F, rug = F, trans = exp, #ylim = ylim2, 
       ylab = "", xlab = expression(rho))
abline(h = 1, lty = 2)

simfit_r2 <- lm(log(r2rat) ~ x_corr + var_cost_var)
simfit_c <- lm(log(crat) ~ var_cost_var + I(log(1-x_corr)))

summary(simfit_r2)
summary(simfit_c)

par(mfrow = c(2,2), mar = c(2, 5, 4, 1.2))
visreg(simfit_r2, "var_cost_var", partial = F, rug = F, trans = exp, #ylim = ylim1, 
       ylab = expression(R[SRL]^2/R[LS]^2), xlab = "",
       main = expression(Varying~cost~variance*","~~kappa))
abline(h = 1, lty = 2)
par(mar = c(2, 3, 4, 2))
visreg(simfit_r2, "x_corr", partial = F, rug = F, trans = exp, #ylim = ylim1, 
       ylab = "", xlab = "",
       main = expression(Varying~covariate~correlation*","~rho))
abline(h = 1, lty = 2)

par(mar = c(5, 5, 2, 1.2))
visreg(simfit_c, "var_cost_var", partial = F, rug = F, trans = exp, #ylim = ylim2, 
       ylab = expression(CPP[SRL]/CPP[LS]), xlab = expression(kappa))
abline(h = 1, lty = 2)
par(mar = c(5, 3, 2, 2))
visreg(simfit_c, "x_corr", partial = F, rug = F, trans = exp, #ylim = ylim2, 
       ylab = "", xlab = expression(rho))
abline(h = 1, lty = 2)

save(rsq, costs, optimal_model_cost, simfit_r2, simfit_c, x_corr, var_cost_var,
     bestfit_c, bestfit_r2, block, crat, r2rat,
     file = save_fname)

## Try a multivariate spline
# fitmv <- mgcv::gam(log(r2rat) ~ s(x_corr, var_cost_var))
# summary(fitmv)
# plot(fitmv)

## Plot both bestfit and simfit together?

par(mfrow = c(2,2), mar = c(2, 5, 4, 1.2))
visreg(simfit_r2, "var_cost_var", partial = F, rug = F, trans = exp, #ylim = ylim1, 
       ylab = expression(R[SRL]^2/R[LS]^2), xlab = "",
       main = expression(Varying~cost~variance*","~~kappa), 
       line.par = list(col = "tomato"),
       fill.par = list(col = rgb(255/255, 99/255,71/255,alpha = .5)), nn=301)
abline(h = 1, lty = 2)

v2 <- visreg(bestfit_r2, "var_cost_var", partial = F, trans = exp, plot = F, nn = 301)
lines(v2$fit$var_cost_var, v2$fit$visregFit, col = "slateblue", lwd = 2)
polygon(c(v2$fit$var_cost_var, rev(v2$fit$var_cost_var)), 
        c(v2$fit$visregLwr, rev(v2$fit$visregUpr)),
        col = rgb(106/255, 90/255, 205/255, alpha = .4), border = NA)

par(mar = c(2, 3, 4, 2))
visreg(simfit_r2, "x_corr", partial = F, rug = F, trans = exp, #ylim = ylim1, 
       ylab = "", xlab = "",
       main = expression(Varying~covariate~correlation*","~rho), 
       line.par = list(col = "tomato"),
       fill.par = list(col = rgb(255/255, 99/255,71/255,alpha = .5)), nn=301)
abline(h = 1, lty = 2)
v2 <- visreg(bestfit_r2, "x_corr", partial = F, trans = exp, plot = F, nn = 301)
lines(v2$fit$x_corr, v2$fit$visregFit, col = "slateblue", lwd = 2)
polygon(c(v2$fit$x_corr, rev(v2$fit$x_corr)), 
        c(v2$fit$visregLwr, rev(v2$fit$visregUpr)),
        col = rgb(106/255, 90/255, 205/255, alpha = .4), border = NA)

par(mar = c(5, 5, 2, 1.2))
visreg(simfit_c, "var_cost_var", partial = F, rug = F, trans = exp, #ylim = ylim2, 
       ylab = expression(CPP[SRL]/CPP[LS]), xlab = expression(kappa), 
       line.par = list(col = "tomato"),
       fill.par = list(col = rgb(255/255, 99/255,71/255,alpha = .5)), nn=301)
abline(h = 1, lty = 2)
v2 <- visreg(bestfit_c, "var_cost_var", partial = F, trans = exp, plot = F, nn = 301)
lines(v2$fit$var_cost_var, v2$fit$visregFit, col = "slateblue", lwd = 2)
polygon(c(v2$fit$var_cost_var, rev(v2$fit$var_cost_var)), 
        c(v2$fit$visregLwr, rev(v2$fit$visregUpr)),
        col = rgb(106/255, 90/255, 205/255, alpha = .4), border = NA)

par(mar = c(5, 3, 2, 2))
visreg(simfit_c, "x_corr", partial = F, rug = F, trans = exp, #ylim = ylim2, 
       ylab = "", xlab = expression(rho), 
       line.par = list(col = "tomato"),
       fill.par = list(col = rgb(255/255, 99/255,71/255,alpha = .5)), nn=301)
abline(h = 1, lty = 2)
v2 <- visreg(bestfit_c, "x_corr", partial = F, trans = exp, plot = F, nn = 301)
lines(v2$fit$x_corr, v2$fit$visregFit, col = "slateblue", lwd = 2)
polygon(c(v2$fit$x_corr, rev(v2$fit$x_corr)), 
        c(v2$fit$visregLwr, rev(v2$fit$visregUpr)),
        col = rgb(106/255, 90/255, 205/255, alpha = .4), border = NA)


## Table 

s1 <- summary(simfit_r2)
coefs <- s1$coefficients[,1]
ses <- s1$coefficients[,2]
ciub <- coefs + 1.96 * ses
cilb <- coefs - 1.96 * ses

prettyc_r2 <- paste0(
  round(exp(coefs),5), 
  " (", 
  round(exp(cilb), 5), 
  ", ",
  round(exp(ciub), 5), 
  ")"
)


s2 <- summary(simfit_c)
coefs <- s2$coefficients[,1]
ses <- s2$coefficients[,2]
ciub <- coefs + 1.96 * ses
cilb <- coefs - 1.96 * ses

prettyc_c <- paste0(
  round(exp(coefs),3), 
  " (", 
  round(exp(cilb), 3), 
  ", ",
  round(exp(ciub), 3), 
  ")"
)

data.frame(
  "$R^2 \\hat \\theta$ (95\\% CI)" = prettyc, 
  "$CPP \\hat \\theta$ (95\\% CI)" = prettyc_c, 
  check.names = FALSE         
)
