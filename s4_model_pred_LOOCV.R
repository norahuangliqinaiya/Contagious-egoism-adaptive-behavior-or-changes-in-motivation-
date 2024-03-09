# The current script corresponds to 
# 2.4 Validating the two models

rm(list=ls())

library(R.matlab)

# Function to constrain behavior within range ----

est_constrain <- function(x){
  
  if (x > 90){
    u <- 90
  } else if (x < 0){
    u <- 0
  } else{
    u <- x
  }
  
  return(u)
}

allocate_constrain <- function(x, radius){
  
  if (x > radius){
    u <- radius
  } else if (x < 0){
    u <- 0
  } else{
    u <- x
  }
  
  return(u)
}

degree_to_radian <- function(d){
  
  r <- d * pi/180
  
  return(r)
  
}


# Simulate allocation behaviors ----

opt_task_param <- read.table('D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m/optimal_task_parameter.txt', header = F)
# opt_task_param <- read.table('D:/data/R/decision_phase/z_ll_m/optimal_task_parameter.txt', header = F)

path <- 'D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m'
# path <- 'D:/data/R/decision_phase/z_ll_m'

tp_idx <- 7
tn <- 45
subj_num <- 200
radius <- 50
sim_time <- 50

lamda <- opt_task_param[tp_idx,1]
dc <- opt_task_param[tp_idx,2]
gn <- opt_task_param[tp_idx,3]

bias_est <- 0
bias_allocate <- 0

param_gt <- matrix(0, subj_num, 7)
behav_est <- matrix(0, subj_num, tn)
behav_allocate_M1 <- matrix(0, subj_num, tn)
behav_allocate_M2 <- matrix(0, subj_num, tn)

for (st in 1:sim_time){
  
  for (subj in 1:subj_num){
    
    setwd('D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m/allocate_offer')
    # setwd('D:/data/R/decision_phase/z_ll_m/allocate_offer')
    
    allocate_offer <- read.table(paste0("lamda_", lamda, "_dc_", dc, "_gn_", gn, ".txt"), header = F)
    
    setwd('D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m/FB')
    # setwd('D:/data/R/decision_phase/z_ll_m/FB')
    
    FB <- read.table(paste0("lamda_", lamda, "_dc_", dc, "_gn_", gn, ".txt"), header = F)
    
    
    phi_gn <- matrix(0, 1, tn)
    phi_gn[1] <- 60
    
    phi_comb <- matrix(0, 1, tn)
    phi_self <- matrix(0, 1, tn+1)
    
    set.seed(st*1000 + subj)
    alpha_gn <- runif(1, min = 0, max = 1)
    beta_gn <- runif(1, min = 0, max = 1)
    beta_0 <- sample(seq(-bias_est,bias_est,1), 1, replace = T)
    
    w_self <- runif(1, min = 0, max = 1)
    beta_self <- runif(1, min = 0, max = 1)
    beta_1 <- sample(seq(-bias_allocate,bias_allocate,1), 1, replace = T)
    
    # define psb range
    phi_self_base <- sample(seq(0,90,1), 1, replace = T)
    
    param_gt[subj, 1] <- alpha_gn
    param_gt[subj, 2] <- beta_gn
    param_gt[subj, 3] <- beta_0
    
    param_gt[subj, 4] <- w_self
    param_gt[subj, 5] <- beta_self
    param_gt[subj, 6] <- beta_1
    param_gt[subj, 7] <- phi_self_base
    
    phi_self[1] <- phi_self_base
    
    for (t in 1:tn){
      
      phi_est <- matrix(0, 1, 4)
      phi_est[1] <- phi_gn[t]
      
      for (tau in 1:3){
        
        phi_est[tau+1] <- phi_est[tau] + alpha_gn * FB[t, tau] * (allocate_offer[t,tau] - phi_est[tau])
        
      }
      
      phi_gn[t+1] <- phi_est[4]
      behav_est[subj,t] <- round(est_constrain(beta_gn * phi_gn[t+1] + beta_0), 2)
      
      # Simulate behavior for M1
      phi_self[t+1] <- (1 - w_self) * phi_self[t] + w_self * phi_gn[t+1]
      behav_allocate_M1[subj,t] <- round(allocate_constrain(beta_self * radius * (1 - cos(degree_to_radian(phi_self[t+1]))) + beta_1, radius), 2)
      
      # Simulate behaviors for M2
      phi_comb[t] <- (1 - w_self) * phi_self_base + w_self * phi_gn[t+1]
      behav_allocate_M2[subj,t] <- round(allocate_constrain(beta_self * radius * (1 - cos(degree_to_radian(phi_comb[t]))) + beta_1, radius), 2)
      
    }
    
  }
  
  fn_param_gt <- paste0(path, '/behav_allocate/model_pred/param_gt_lamda_', lamda, '_dc_', dc, '_gn_', gn, '_st_',st,'.txt')
  write.table(param_gt, fn_param_gt, row.names = F, col.names = F)
  
  # psb full range: save behavior for M1
  fn_behav_allocate <- paste0(path, '/behav_allocate/model_pred/gt_M1/lamda_', lamda, '_dc_', dc, '_gn_', gn, '_st_',st,'.txt')
  write.table(behav_allocate_M1, fn_behav_allocate, row.names = F, col.names = F)
  
  # psb full range: save behavior for M2
  fn_behav_allocate <- paste0(path, '/behav_allocate/model_pred/gt_M2/lamda_', lamda, '_dc_', dc, '_gn_', gn, '_st_',st,'.txt')
  write.table(behav_allocate_M2, fn_behav_allocate, row.names = F, col.names = F)
  
}

hist(param_gt[,4])
hist(param_gt[,5])
hist(param_gt[,6])
hist(param_gt[,7])

colnames(param_gt) <- c('alpha_gn', 'beta_gn', 'beta_0',
                        'w_self', 'beta_self', 'beta_1', 
                        'phi_self_base')


# Perform model estimation for the decision process in Matlab ----


# LOO with cross-validation ----

# Function to plot normal distribution curves of two columns with peaks indicated
plot_normal_distribution_curves <- function(data, colnames, colors, p_title) {
  par(mfrow = c(1, 1))  # Set up the plot layout
  
  mean1 <- mean(data[, 1])
  sd1 <- sd(data[, 1])
  x_vals1 <- seq(min(data[, 1]), max(data[, 1]), length.out = 1000)
  y_vals1 <- dnorm(x_vals1, mean = mean1, sd = sd1)
  
  mean2 <- mean(data[, 2])
  sd2 <- sd(data[, 2])
  x_vals2 <- seq(min(data[, 2]), max(data[, 2]), length.out = 1000)
  y_vals2 <- dnorm(x_vals2, mean = mean2, sd = sd2)
  
  plot(x_vals1, y_vals1, col = colors[1], type = "l", main = p_title, 
       xlab = "Value", ylab = "Density", xlim = range(data), 
       ylim = c(0, max(y_vals1, y_vals2) * 1.1), lwd = 2)  # Plot normal distribution curve for column 1
  lines(x_vals2, y_vals2, col = colors[2], lwd = 2)  # Add normal distribution curve for column 2
  
  abline(v = mean1, col = colors[1], lwd = 2, lty = 3)  # Vertical line at mean of column 1
  abline(v = mean2, col = colors[2], lwd = 2, lty = 4)  # Vertical line at mean of column 2
  
  legend("topright", legend = colnames, col = colors, lty = c(1,1), lwd = 2)  # Add legend
}


# Read task stimulus and allocation behaviors ---

library(doParallel)

path <- 'D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m'
# path <- 'D:/data/R/decision_phase/z_ll_m'

opt_task_param <- read.table(paste0(path, '/optimal_task_parameter.txt'), header = F)

tp_idx <- 7
tn <- 45
radius <- 50
sim_time <- 50

bias_est <- 0
bias_allocate <- 0

lamda <- opt_task_param[tp_idx,1]
dc <- opt_task_param[tp_idx,2]
gn <- opt_task_param[tp_idx,3]

gt_model <- 1 # 1: change of motive; 2: change of behavior
sample_size <- seq(50,200,10)
es_mse <- matrix(0, sim_time, length(sample_size))
es_ll <- matrix(0, sim_time, length(sample_size))

numCores <- detectCores()
registerDoParallel(numCores)

# Perform LOOCV and generate MSEs ---

foreach (ss = 1:length(sample_size), .combine = rbind) %dopar% {
  
  library("boot")
  library("R.matlab")
  library("BayesFactor")
  
  set.seed(155)
  subj_pool <- sample(1:200, sample_size[ss], replace = F)
  
  for (st in 1:sim_time){
    
    # Read stimulus from gt_M1 (the change of motive model) or gt_M2 (the change of behavior model)
    setwd(paste0('D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m/behav_allocate/model_pred/gt_M', gt_model))
    # setwd(paste0('D:/data/R/decision_phase/z_ll_m/behav_allocate/model_pred/gt_M', gt_model))
    
    fn <- paste0('lamda_', lamda, '_dc_', dc, '_gn_', gn, '_st_',st,'.txt')
    behav_allocate <- as.matrix(read.table(fn, header = F))
    
    setwd('D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m/allocate_offer')
    # setwd('D:/data/R/decision_phase/z_ll_m/allocate_offer')
    
    allocate_offer <- read.table(paste0("lamda_", lamda, "_dc_", dc, "_gn_", gn, ".txt"), header = F)
    
    setwd('D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m/FB')
    # setwd('D:/data/R/decision_phase/z_ll_m/FB')
    
    FB <- read.table(paste0("lamda_", lamda, "_dc_", dc, "_gn_", gn, ".txt"), header = F)
    
    fn_param_gt <- paste0(path, '/behav_allocate/model_pred/param_gt_lamda_', lamda, '_dc_', dc, '_gn_', gn, '_st_',st,'.txt')
    param_gt <- read.table(fn_param_gt, header = F)
    
    cv.mse <- matrix(0, sample_size[ss], 2)
    
    for (idx in 1:sample_size[ss]){
      
      behav <- as.vector(behav_allocate[subj_pool[idx],])
      
      alpha_gn <- param_gt[subj_pool[idx], 1]
      beta_gn <- param_gt[subj_pool[idx], 2]
      beta_0 <- param_gt[subj_pool[idx], 3]
      phi_self_base <- param_gt[subj_pool[idx], 7]
      
      
      # read parameter estimates from M1
      setwd(paste0('D:/postdoc/work/UKW/projects/egoism/experiment/data/matlab/z_ll_m/decision_phase/est_data/gt_M',gt_model,'/M1/model_pred/beta_free_est'))
      # setwd(paste0('D:/data/matlab/z_ll_m/decision_phase/est_data/gt_M',gt_model,'/M1/model_pred/beta_free_est'))
      
      fn <- paste0('subj_', subj_pool[idx], '_lamda_', lamda, '_dc_', dc, '_gn_', gn, '_st_',st, '.mat')
      data <- matrix(unlist(readMat(fn)), ncol = 4)
      colnames(data) <- c("ll", "alpha_w", "beta_self", "beta_1")
      
      w_self_M1 <- data[2]
      beta_self_M1 <- data[3]
      beta_1_M1 <- data[4]
      
      # read parameter estimates from M2
      setwd(paste0('D:/postdoc/work/UKW/projects/egoism/experiment/data/matlab/z_ll_m/decision_phase/est_data/gt_M',gt_model,'/M2/model_pred/beta_free_est'))
      # setwd(paste0('D:/data/matlab/z_ll_m/decision_phase/est_data/gt_M',gt_model,'/M2/model_pred/beta_free_est'))
      
      fn <- paste0('subj_', subj_pool[idx], '_lamda_', lamda, '_dc_', dc, '_gn_', gn, '_st_',st, '.mat')
      data <- matrix(unlist(readMat(fn)), ncol = 4)
      colnames(data) <- c("ll", "alpha_w", "beta_self", "beta_1")
      
      w_self_M2 <- data[2]
      beta_self_M2 <- data[3]
      beta_1_M2 <- data[4]
      
      # calculate latent variables
      phi_gn <- matrix(0, 1, tn+1)
      phi_comb <- matrix(0, 1, tn)
      phi_self <- matrix(0, 1, tn+1)
      
      phi_gn[1] <- 60
      phi_self[1] <- phi_self_base
      
      for (t in 1:tn){
        
        phi_est <- matrix(0, 1, 4)
        phi_est[1] <- phi_gn[t]
        
        for (tau in 1:3){
          
          phi_est[tau+1] <- phi_est[tau] + alpha_gn * FB[t, tau] * (allocate_offer[t,tau] - phi_est[tau])
          
        }
        
        phi_gn[t+1] <- phi_est[4]
        
        # latent variable for M1
        phi_self[t+1] <- (1 - w_self_M1) * phi_self[t] + w_self_M1 * phi_gn[t+1]
        
        # latent variable for M2
        phi_comb[t] <- (1 - w_self_M2) * phi_self_base + w_self_M2 * phi_gn[t+1]
        
      }
      
      df_M <- data.frame(behav, phi_self[2:(tn+1)], phi_comb[1:tn])
      colnames(df_M) <- c("behav", "phi_self", "phi_comb")
      behav.loocv.M1 <- glm(behav ~ phi_self, data = df_M)
      cv.mse[idx,1] <- cv.glm(df_M, behav.loocv.M1)$delta[1]
      
      behav.loocv.M2 <- glm(behav ~ phi_comb, data = df_M)
      cv.mse[idx,2] <- cv.glm(df_M, behav.loocv.M2)$delta[1]
      
    }
    
    fn_cv_mse <- paste0(path, '/mse/sample_size_', ss, '_st_', st, '_M', gt_model,'.txt')
    write.table(cv.mse, fn_cv_mse, row.names = F, col.names = F)
    
  }
}


# Read negative log likelihood data

foreach (ss = 1:length(sample_size), .combine = rbind) %dopar% {
  
  set.seed(155)
  subj_pool <- sample(1:200, sample_size[ss], replace = F)
  
  for (st in 1:sim_time){
    
    neg_loglik <- matrix(0, sample_size[ss], 2)
    
    for (idx in 1:sample_size[ss]) {
      
      library(R.matlab)
      
      # read parameter estimates from M1
      setwd(paste0('D:/postdoc/work/UKW/projects/egoism/experiment/data/matlab/z_ll_m/decision_phase/est_data/gt_M',gt_model,'/M1/model_pred/beta_free_est'))
      # setwd(paste0('D:/data/matlab/z_ll_m/decision_phase/est_data/gt_M',gt_model,'/M1/model_pred/beta_free_est'))
      
      fn <- paste0('subj_', subj_pool[idx], '_lamda_', lamda, '_dc_', dc, '_gn_', gn, '_st_',st, '.mat')
      data <- matrix(unlist(readMat(fn)), ncol = 4)
      colnames(data) <- c("ll", "alpha_w", "beta_self", "beta_1")
      neg_loglik[idx,1] <- data[1]
      
      # read parameter estimates from M2
      setwd(paste0('D:/postdoc/work/UKW/projects/egoism/experiment/data/matlab/z_ll_m/decision_phase/est_data/gt_M',gt_model,'/M2/model_pred/beta_free_est'))
      # setwd(paste0('D:/data/matlab/z_ll_m/decision_phase/est_data/gt_M',gt_model,'/M2/model_pred/beta_free_est'))
      
      fn <- paste0('subj_', subj_pool[idx], '_lamda_', lamda, '_dc_', dc, '_gn_', gn, '_st_',st, '.mat')
      data <- matrix(unlist(readMat(fn)), ncol = 4)
      colnames(data) <- c("ll", "alpha_w", "beta_self", "beta_1")
      
      neg_loglik[idx,2] <- data[1]
      
      fn_neg_loglik <- paste0(path, '/neg_loglik/sample_size_', ss, '_st_', st, '_M', gt_model,'.txt')
      write.table(neg_loglik, fn_neg_loglik, row.names = F, col.names = F)
      
    }
  }
}


# Calculate Bayes Factor ----

# M1: the Change of Motive model
# M2: the Change of Behavior model
gt_model <- 1

path <- 'D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m'
# path <- 'D:/data/R/decision_phase/z_ll_m'

library(BayesFactor)
library(lsr)
library(ggplot2)

sample_size <- seq(50,200,10)
sim_time <- 50

es_mse_bf <- matrix(0, sim_time, length(sample_size))
es_mse_d <- matrix(0, sim_time, length(sample_size))

es_nll_bf <- matrix(0, sim_time, length(sample_size))
es_nll_d <- matrix(0, sim_time, length(sample_size))

for (ss in 1:length(sample_size)) {
  
  set.seed(155)
  subj_pool <- sample(1:200, sample_size[ss], replace = F)
  
  for (st in 1:sim_time){
    
    fn_cv_mse <- paste0(path, '/mse/sample_size_', ss, '_st_', st, '_M', gt_model,'.txt')
    cv.mse <- as.matrix(read.table(fn_cv_mse, header = F))
    es_mse_d[st, ss] <- cohensD(cv.mse[,2], cv.mse[,1], method = "paired")
    es_mse_bf[st, ss] <- extractBF(ttestBF(cv.mse[,1], cv.mse[,2]))$bf
    
    fn_neg_loglik <- paste0(path, '/neg_loglik/sample_size_', ss, '_st_', st, '_M', gt_model,'.txt')
    neg_loglik <- as.matrix(read.table(fn_neg_loglik, header = F))
    es_nll_d[st, ss] <- cohensD(neg_loglik[,2], neg_loglik[,1], method = "paired")
    es_nll_bf[st, ss] <- extractBF(ttestBF(neg_loglik[,1], neg_loglik[,2]))$bf
    
  }
}

# Density curves of MSE

if (gt_model == 1) {
  gt_model_fig = 2
} else {
  gt_model_fig = 1
}

column_names <- c("M2", "M1")
curve_colors <- c("blue", "red")

plot_normal_distribution_curves(log(cv.mse), column_names, curve_colors, paste0("Distribution of MSE (Groundtruth M", gt_model_fig, ")"))

# Density curves of negative logliklihood
column_names <- c("M2", "M1")
curve_colors <- c("blue", "red")

plot_normal_distribution_curves(log(neg_loglik), column_names, curve_colors , paste0("Distribution of Negative LL (Groundtruth M", gt_model_fig, ")"))


# Plot of Power analysis (BF ~ Sample size)

windowsFonts()

es_mse_bf_avg <- colMeans(es_mse_bf)
es_mse_bf_plot <- data.frame(x = sample_size, y = es_mse_bf_avg)

ggplot(es_mse_bf_plot[1:6,], aes(x, y)) + 
  geom_line() + 
  geom_text(aes(label = sprintf("%.2f", y)), 
            vjust = c(rep(-1,5), 0), 
            hjust = c(0.1,0.5,0.5,0.6,1,1),size = 4,
            family = "sans") + 
  labs(x = "Sample size",y = "BF10", title = "MSE") + 
  scale_x_continuous(breaks = seq(50, 100, 10),limits = c(50, 100)) + 
  theme_bw() + 
  theme(
    axis.text = element_text(family = "sans", size = 12),
    axis.title.y = element_text(family = "sans", size = 15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(family = "sans", size = 15, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill=NA, colour = "black", linewidth=1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  )



