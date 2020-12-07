
# this script generates data for simulation
# three types of data are generated
# each type contains 1000 simulated data set

# data 1: cor within group = cor between group = 0 
# data 2: cor within group = 0.5, cor between group = 0
# data 3: cor within group = cor between group = 0.5
set.seed(12345)
true_beta <- c(1, 2, 3, 4, 5,
               0, 0, 0, 0, 0,
               1, 2, 3, 4, 5,
               rep(0, 35))

##### function to generate data #####

gen_data <- function(iter = 1000, n1 = 200, n2 = 200, J = 10, K = 5, true_beta,
                     family = "gaussian", SNR = 1, 
                     rho = 0, rho.g = 0){
  sim_data <- list()
  for(i in 1:iter){
    sim_data[[i]] <- list()
    this_sim1 <- genDataGrp(n = n1, J = J, K = K, beta = true_beta,
                            family = family, SNR = SNR, 
                            rho = rho, rho.g = rho.g)
    this_sim2 <- genDataGrp(n = n2, J = J, K = K, beta = true_beta,
                            family = family, SNR = SNR, 
                            rho = rho, rho.g = rho.g)
    sim_data[[i]]$train <- list(X = this_sim1$X, y = this_sim1$y)
    sim_data[[i]]$test <- list(X = this_sim2$X, y = this_sim2$y)
  }
  return(sim_data)
}

##### no corrlation #####

sim_data1 <- gen_data(true_beta = c(1, 2, 3, 4, 5,
                                    0, 0, 0, 0, 0,
                                    1, 2, 3, 4, 5,
                                    rep(0, 35)), rho = 0, rho.g = 0)


##### within correlation ####
sim_data2 <- gen_data(true_beta = c(1, 2, 3, 4, 5,
                                    0, 0, 0, 0, 0,
                                    1, 2, 3, 4, 5,
                                    rep(0, 35)), rho = 0, rho.g = 0.5)

##### within & between correlation ####
sim_data3 <- gen_data(true_beta = c(1, 2, 3, 4, 5,
                                    0, 0, 0, 0, 0,
                                    1, 2, 3, 4, 5,
                                    rep(0, 35)), rho = 0.5, rho.g = 0.5)

save(sim_data1, sim_data2, sim_data3, file = "Data/sim_data.RData")

