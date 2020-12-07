

library(ncvreg)
library(grpreg)
library(hdrm)
library(dplyr)

set.seed(123)


##### no correlations #####

pe_list1 <- list()
cpp_list1 <- list()

# rho = 0, rho.g = 0
for(iter in 1:1000){
  source("Code/sim_bic.R")
  pe_list1[[iter]] <- pred_error
  cpp_list1[[iter]] <- cpp
}

pe_df1 <- pe_list1 %>% bind_rows() 
cpp_df1 <- cpp_list1 %>% bind_rows() 

lapply(pe_df1, mean)
lapply(cpp_df1, mean)

save(pe_df1, cpp_df1, file = "Data/sim_bic_no_cor.RData")

##### within-group correlations #####

pe_list2 <- list()
cpp_list2 <- list()

# rho = 0, rho.g = 0.5
for(iter in 1:1000){
  source("Code/sim_bic.R")
  pe_list2[[iter]] <- pred_error
  cpp_list2[[iter]] <- cpp
}

pe_df2 <- pe_list2 %>% bind_rows() 
cpp_df2 <- cpp_list2 %>% bind_rows() 

lapply(pe_df2, mean)
lapply(cpp_df2, mean)

save(pe_df2, cpp_df2, file = "Data/sim_bic_within_cor.RData")

##### both between and within group correlations #####

pe_list3 <- list()
cpp_list3 <- list()

# rho = 0.5, rho.g = 0.5
for(iter in 1:1000){
  source("Code/sim_bic.R")
  pe_list3[[iter]] <- pred_error
  cpp_list3[[iter]] <- cpp
}

pe_df3 <- pe_list3 %>% bind_rows() 
cpp_df3 <- cpp_list3 %>% bind_rows() 

lapply(pe_df3, mean)
lapply(cpp_df3, mean)

save(pe_df3, cpp_df3, file = "Data/sim_bic_both_cor.RData")

##### one iteration #####
source("Code/sim_bic.R")

cor(x_grp_all$X)

grp1 <- which(sapply(split(coef(lasso.fit, lambda = l1)[-1]!=0, f = grp), sum)!=0)
grp1
grp2 
grp3
grp4
