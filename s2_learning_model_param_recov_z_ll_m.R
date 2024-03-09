# The current script corresponds to Step 2 (the conventional parameter recovery analysis) of 
# 2 The simulation protocol
# 2.1 Establishing the egoism learning task

rm(list=ls())

# Read task parameters remained after prescreen ----

setwd('D:/postdoc/work/UKW/projects/egoism/experiment/data/R/learning_phase')
# setwd('D:/data/R/learning_phase')
tp_cond_sub <- as.data.frame(read.table('lik_all_with_bias_2_thres_0.85_zero_fb_5_50_z_ll_m.txt', header = T))

# Perform model estimation for the learning process in Matlab ----


# Read model estimates from Matlab output ----

# Read missing data of model estimation
setwd('D:/postdoc/work/projects/egoism/experiment/data/matlab/learning_phase/est_data')
# setwd('D:/data/matlab/learning_phase/est_data')

missing_data <- as.data.frame(read.delim('missing_data.txt', header = F))[,1:3]
colnames(missing_data) <- c("lamda", 'dc', 'gn')


# Read learning parameters -----
setwd('D:/postdoc/work/UKW/projects/egoism/experiment/data/R/learning_phase/z_ll_m/gt_param')
# setwd('D:/data/R/learning_phase/z_ll_m/gt_param')

param_gt <- as.data.frame(read.table('param_gt.txt', header = F))
colnames(param_gt) <- c("alpha_gn", "beta_gn", "beta_0")

# Remove task parameter that has missing value in parameter estimation
library(dplyr)
tp_cond_sub <- anti_join(tp_cond, missing_data, by = c("lamda", "dc", "gn"))

path <- 'D:/postdoc/work/UKW/projects/egoism/experiment/data/matlab/z_ll_m/learning_phase/est_data'
# path <- 'D:/matlab/z_ll_m/learning_phase/est_data'



# Step 2: Parameter Recovery ----
library(R.matlab)
library(BayesFactor)

param_est <- as.data.frame(matrix(0, 30, 3))
colnames(param_est) <- c("alpha_gn", "beta_gn", "beta_0")

param_recov <- matrix(0, dim(tp_cond_sub)[1], 9)
colnames(param_recov) <- c("lamda", "dc", "gn", "cor_alpha_gn", "cor_beta_gn", "cor_beta_0",
                           "diff_alpha_gn", "diff_beta_gn", "diff_beta_0")


for (i in 1:dim(tp_cond_sub)[1]){
  
  for (sub in 1:30){
    
    fn <- paste0(path, '/subj_', sub, '_lamda_', tp_cond_sub[i,1], '_dc_', tp_cond_sub[i,2], '_gn_', tp_cond_sub[i,3], '.mat')
    data <- matrix(unlist(readMat(fn)), ncol = 4)
    colnames(data) <- c("ll", "alpha_gn", "beta_gn", "beta_0")
    
    param_est[sub, 1] <- data[1, 2]
    param_est[sub, 2] <- data[1, 3]
    param_est[sub, 3] <- data[1, 4]
    
  }
  
  param_recov[i, 1] <- tp_cond_sub[i, 1]
  param_recov[i, 2] <- tp_cond_sub[i, 2]
  param_recov[i, 3] <- tp_cond_sub[i, 3]
  
  param_recov[i, 4] <- cor(param_gt$alpha_gn, param_est$alpha_gn)
  param_recov[i, 5] <- cor(param_gt$beta_gn, param_est$beta_gn)
  param_recov[i, 6] <- cor(param_gt$beta_0, param_est$beta_0)
  
  param_recov[i, 7] <- extractBF(ttestBF(x = param_gt$alpha_gn, y = param_est$alpha_gn))$bf
  param_recov[i, 8] <- extractBF(ttestBF(x = param_gt$beta_gn, y = param_est$beta_gn))$bf
  param_recov[i, 9] <- extractBF(ttestBF(x = param_gt$beta_0, y = param_est$beta_0))$bf
  
}


# Save optimal task parameters -----
# Save optimal task parameters that give a BF10 smaller than 1 and 
# a correlation coefficient with the groundtruth parameter larger than 0.8
param_recov <- as.data.frame(param_recov)
tp_opt <- subset(param_recov, cor_alpha_gn >= 0.8 & cor_beta_gn >= 0.8 & diff_alpha_gn <= 1 & diff_beta_gn <= 1)

setwd('D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m')
# setwd('D:/data/R/decision_phase/z_ll_m')
fn <- 'optimal_task_parameter.txt'

write.table(tp_opt, fn, row.names = F, col.names = F)