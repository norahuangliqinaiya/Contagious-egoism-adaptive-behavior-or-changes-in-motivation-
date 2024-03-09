# The current script corresponds to 
# 2 The simulation protocol
# 2.2 Distinguishing the two decision models

rm(list=ls())

library(R.matlab)

# Function to constrain behavior within range ----

allocate_constrain <- function(x,radius){
  
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

# Simulate behaviors given all levels of decision parameters ----

setwd('D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m')
# setwd('D:/data/R/decision_phase/z_ll_m')
result <- read.table("optimal_task_parameter.txt",header=F)
tn <- 45
radius <- 50

path <- 'D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m'
# path <- 'D:/data/R/decision_phase/z_ll_m'

# alpha_gn; alpha_self; beta_self; phi_self_base

agn <- c(0.1, 0.3, 0.5, 0.9)
aw <- c(0.1, 0.3, 0.5, 0.9)
bs <- c(0.1, 0.3, 0.5, 0.9)
psb <- c(15, 30, 45, 60, 75, 90)


dp <- cbind(alpha_gn = rep(agn, each = length(aw) * length(bs) * length(psb)),
            alpha_w = rep(rep(aw, each = length(bs) * length(psb)), time = length(agn)),
            beta_self = rep(rep(bs, each = length(psb)), time = length(aw) * length(agn)),
            phi_self_base = rep(psb, time = length(bs) * length(aw) * length(agn)))


for (i in 1:dim(result)[1]){
  
  lamda <- result[i,1]
  dc <- result[i,2]
  gn <- result[i,3]
  
  # -------------------
  # Generate group norm
  # -------------------
  
  set.seed(1000+gn)
  phi_gt_1 <- sample(0:20, size = 1, replace = T)
  
  set.seed(2000+gn)
  phi_gt_2 <- sample(20:40, size = 1, replace = T)
  
  set.seed(3000+gn)
  phi_gt_3 <- sample(40:60, size = 1, replace = T)
  
  # -------------------------
  # Generate allocation plans
  # -------------------------
  
  set.seed(1000+gn)
  allocate_offer_1 <- as.matrix(phi_gt_1 + rpois(tn, lamda), tn, 1)
  
  set.seed(2000+gn)
  allocate_offer_2 <- as.matrix(phi_gt_2 + rpois(tn, lamda), tn, 1)
  
  set.seed(3000+gn)
  allocate_offer_3 <- as.matrix(phi_gt_3 + rpois(tn, lamda), tn, 1)
  
  allocate_offer <- cbind(allocate_offer_1, allocate_offer_2, allocate_offer_3)
  
  # -----------------
  # Simulate feedback
  # -----------------
  
  FB <- matrix(0, tn, 3)
  
  for (t in 1:tn){
    
    if ((allocate_offer[t,1] >= phi_gt_1 - dc) && (allocate_offer[t,1] <= phi_gt_1 + dc)){
      FB[t, 1] <- 1
    } else {
      FB[t, 1] <- 0
    }
    
    if ((allocate_offer[t,2] >= phi_gt_2 - dc) && (allocate_offer[t,2] <= phi_gt_2 + dc)){
      FB[t, 2] <- 1
    } else {
      FB[t, 2] <- 0
    }
    
    if ((allocate_offer[t,3] >= phi_gt_3 - dc) && (allocate_offer[t,3] <= phi_gt_3 + dc)){
      FB[t, 3] <- 1
    } else {
      FB[t, 3] <- 0
    }
  }
  
  allocate_offer_new <- matrix(0, tn, 3)
  FB_new <- matrix(0, tn, 3)
  
  set.seed(522)
  
  for (t in 1:tn){
    order <- sample(c(1,2,3), replace = F)
    allocate_offer_new[t,] <- allocate_offer[t, order]
    FB_new[t,] <- FB[t, order]
  }
  
  fn_fb <- paste0(path, '/FB/lamda_', lamda, '_dc_', dc, '_gn_', gn, '.txt')
  write.table(FB_new, fn_fb, row.names = F, col.names = F)
  
  fn_offer <- paste0(path, '/allocate_offer/lamda_', lamda, '_dc_', dc, '_gn_', gn, '.txt')
  write.table(allocate_offer_new, fn_offer, row.names = F, col.names = F)
  
  
  # ----------------------------
  # Simulate allocation decisions
  # ----------------------------
  
  phi_comb <- matrix(0, 1, tn)
  phi_self <- matrix(0, 1, tn+1)
  
  behav_allocate_M1 <- matrix(0, dim(dp)[1], tn+4)
  behav_allocate_M2 <- matrix(0, dim(dp)[1], tn+4)
  
  phi_gn <- matrix(0, 1, tn)
  phi_gn[1] <- 60
  
  for (subj in 1:dim(dp)[1]){
    
    set.seed(subj)
    self_bias <- runif(1, min = 0, max = 0)
    
    behav_allocate_M1[subj,1] <- dp[subj,1]
    behav_allocate_M1[subj,2] <- dp[subj,2]
    behav_allocate_M1[subj,3] <- dp[subj,3]
    behav_allocate_M1[subj,4] <- dp[subj,4]
    
    behav_allocate_M2[subj,1] <- dp[subj,1]
    behav_allocate_M2[subj,2] <- dp[subj,2]
    behav_allocate_M2[subj,3] <- dp[subj,3]
    behav_allocate_M2[subj,4] <- dp[subj,4]
    
    phi_self[1] <- dp[subj,4]
    
    for (t in 1:tn){
      
      phi_est <- matrix(0, 1, 4)
      phi_est[1] <- phi_gn[t]
      
      for (tau in 1:3){
        
        phi_est[tau+1] <- phi_est[tau] + dp[subj,1] * FB[t, tau] * (allocate_offer[t,tau] - phi_est[tau])
        
      }
      
      phi_gn[t+1] <- phi_est[4]
      
      # simulate behavior for M1
      phi_self[t+1] <- (1 - dp[subj,2]) * phi_self[t] + dp[subj,2] * phi_gn[t+1]
      behav_allocate_M1[subj,t+4] <- round(allocate_constrain(dp[subj,3] * radius * (1 - cos(degree_to_radian(phi_self[t+1]))) + self_bias, radius), 2)
      
      # simulate behavior for M2
      phi_comb[t] <- (1 - dp[subj,2]) * dp[subj,4] + dp[subj,2] * phi_gn[t+1]
      behav_allocate_M2[subj,t+4] <- round(allocate_constrain(dp[subj,3] * radius * (1 - cos(degree_to_radian(phi_comb[t]))) + self_bias, radius), 2)
      
    }
    
    fn_behav_allocate_M1 <- paste0(path, '/behav_allocate/model_dissociate/gt_M1/lamda_', lamda, '_dc_', dc, '_gn_', gn, '.txt')
    write.table(behav_allocate_M1, fn_behav_allocate_M1, row.names = F, col.names = F)
    
    fn_behav_allocate_M2 <- paste0(path, '/behav_allocate/model_dissociate/gt_M2/lamda_', lamda, '_dc_', dc, '_gn_', gn, '.txt')
    write.table(behav_allocate_M2, fn_behav_allocate_M2, row.names = F, col.names = F)
  }
}


# Perform model estimation for the decision process in Matlab ----


# To look for tp that distinguish the two models the most ----

library(R.matlab)

opt_task_param <- read.table('D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m/optimal_task_parameter.txt', header = F)
# opt_task_param <- read.table('D:/data/R/decision_phase/z_ll_m/optimal_task_parameter.txt', header = F)

path <- 'D:/postdoc/work/UKW/projects/egoism/experiment/data/matlab/z_ll_m/decision_phase/est_data'
# path <- 'D:/data/matlab/z_ll_m/decision_phase/est_data'


agn <- c(0.1, 0.3, 0.5, 0.9)
aw <- c(0.1, 0.3, 0.5, 0.9)
bs <- c(0.1, 0.3, 0.5, 0.9)
psb <- c(15, 30, 45, 60, 75, 90)

subj_num <- length(agn) * length(aw) * length(bs) * length(psb)

lik_diff <- matrix(0, dim(opt_task_param)[1],subj_num)
lik_diff_rank <- matrix(0, dim(opt_task_param)[1],subj_num)

gt_model <- 'gt_M1'

for (subj in 1:subj_num){
  
  for (i in (1:dim(opt_task_param)[1])){
    
    fn_M1 <- paste0(path, '/',gt_model, '/M1/model_dissociate/subj_', subj, '_lamda_', opt_task_param[i,1], '_dc_', opt_task_param[i,2], '_gn_', opt_task_param[i,3], '.mat')
    data_M1 <- matrix(unlist(readMat(fn_M1)), ncol = 4)
    colnames(data_M1) <- c("ll", "alpha_gn", "beta_gn", "beta_0")
    
    fn_M2 <- paste0(path, '/',gt_model, '/M2/model_dissociate/subj_', subj, '_lamda_', opt_task_param[i,1], '_dc_', opt_task_param[i,2], '_gn_', opt_task_param[i,3], '.mat')
    data_M2 <- matrix(unlist(readMat(fn_M2)), ncol = 4)
    colnames(data_M2) <- c("ll", "alpha_gn", "beta_gn", "beta_0")
    
    lik_diff[i, subj] <- data_M1[,1] - data_M2[,1]
    
  }
  
  lik_diff_rank[order(lik_diff[,subj], decreasing = T), subj] <- seq(1, dim(opt_task_param)[1], 1)
}


rank_sum <- rowSums(lik_diff_rank)
rank_sum_f <- cbind(tp_id = seq(1,length(rank_sum),1), rank_sum)
print(rank_sum_f[order(rank_sum_f[,2], decreasing = F),])


# Check parameter recovery for the selected TP ----

# -------------------
# Simulate behaviors
# ------------------

rm(list=ls())

opt_task_param <- read.table('D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m/optimal_task_parameter.txt', header = F)
# opt_task_param <- read.table('D:/data/R/decision_phase/z_ll_m/optimal_task_parameter.txt', header = F)


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

path <- 'D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m'
# path <- 'D:/data/R/decision_phase/z_ll_m'


tp_idx <- 7
tn <- 45
subj_num <- 30
radius <- 50

lamda <- opt_task_param[tp_idx,1]
dc <- opt_task_param[tp_idx,2]
gn <- opt_task_param[tp_idx,3]

bias_est <- 0
bias_allocate <- 0

param_gt <- matrix(0, subj_num, 7)
behav_est <- matrix(0, subj_num, tn)
behav_allocate_M1 <- matrix(0, subj_num, tn)  
behav_allocate_M2 <- matrix(0, subj_num, tn)  

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
  
  set.seed(subj)
  alpha_gn <- runif(1, min = 0, max = 1)
  beta_gn <- runif(1, min = 0, max = 1)
  beta_0 <- sample(seq(-bias_est,bias_est,1), 1, replace = T)
  
  w_self <- runif(1, min = 0, max = 1)
  beta_self <- runif(1, min = 0, max = 1)
  beta_1 <- sample(seq(-bias_allocate,bias_allocate,1), 1, replace = T)
  
  # Define psb range
  phi_self_base <- sample(seq(0,90,1), 1, replace = T)
  
  param_gt[subj, 1] <- alpha_gn
  param_gt[subj, 2] <- beta_gn
  param_gt[subj, 3] <- beta_0
  
  param_gt[subj, 4] <- w_self
  param_gt[subj, 5] <- beta_self
  param_gt[subj, 6] <- beta_1
  param_gt[subj, 7] <- phi_self_base
  

  fn_param_gt <- paste0(path, '/behav_allocate/param_recov/param_gt_lamda_', lamda, '_dc_', dc, '_gn_', gn, '.txt')
  write.table(param_gt, fn_param_gt, row.names = F, col.names = F)
  
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
  
  # Save behavior for M1
  fn_behav_allocate_M1 <- paste0(path, '/behav_allocate/param_recov/gt_M1/lamda_', lamda, '_dc_', dc, '_gn_', gn, '.txt')
  write.table(behav_allocate_M1, fn_behav_allocate_M1, row.names = F, col.names = F)
  
  # Save behavior for M2
  fn_behav_allocate_M2 <- paste0(path, '/behav_allocate/param_recov/gt_M2/lamda_', lamda, '_dc_', dc, '_gn_', gn, '.txt')
  write.table(behav_allocate_M2, fn_behav_allocate_M2, row.names = F, col.names = F)
}


# Perform model estimation for the decision process in Matlab ----


# Parameter recovery ----
tp_idx <- 7
tn <- 45
subj_num <- 30
radius <- 50

param_recov <- matrix(0, subj_num, 3)

lamda <- opt_task_param[tp_idx,1]
dc <- opt_task_param[tp_idx,2]
gn <- opt_task_param[tp_idx,3]

setwd('D:/postdoc/work/UKW/projects/egoism/experiment/data/R/decision_phase/z_ll_m/behav_allocate/param_recov')
# setwd('D:/data/R/decision_phase/z_ll_m/behav_allocate/param_recov')

param_gt <- read.table(paste0("param_gt_lamda_", lamda, "_dc_", dc, "_gn_", gn, ".txt"), header = F)


for (m in 1:2){
  
  for (subj in 1:subj_num){

    # read data
    setwd(paste0('D:/postdoc/work/UKW/projects/egoism/experiment/data/matlab/z_ll_m/decision_phase/est_data/gt_M', m, '/M', m, '/param_recov'))
    # setwd(paste0('D:/data/matlab/z_ll_m/decision_phase/est_data/gt_M', m, '/M', m, '/param_recov'))
    
    fn <- paste0('subj_', subj, '_lamda_', lamda, '_dc_', dc, '_gn_', gn, '.mat')
    data <- matrix(unlist(readMat(fn)), ncol = 4)
    colnames(data) <- c("ll", "alpha_w", "beta_self", "beta_1")
    
    param_recov[subj,1] <- data[2]
    param_recov[subj,2] <- data[3]
    param_recov[subj,3] <- data[4]
    
  }
  
  w_recov_df <- data.frame(x = param_gt[,4], y = param_recov[,1])
  p1 <- ggplot(w_recov_df, aes(x, y)) + 
    geom_point(shape = 21, stroke = 1, color = "black") + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(x = "", y = "") +  
    scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) + 
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) + 
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
 
  
  beta_recov_df <- data.frame(x = param_gt[,5], y = param_recov[,2])
  p2 <- ggplot(beta_recov_df, aes(x, y)) + 
    geom_point(shape = 21, stroke = 1, color = "black") + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(x = "", y = "") +  
    scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) + 
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) + 
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
  
  #setwd('D:/postdoc/work/UKW/projects/egoism/preregistration/maniscript/figs')
  #ggsave(paste0('param_recov_w_gt_m',m,'.jpg'), p1, width = 5, height = 5, units = c('in'))
  #ggsave(paste0('param_recov_beta_gt_m',m,'.jpg'), p2, width = 5, height = 5, units = c('in'))
  
  print(cor.test(param_gt[,4], param_recov[,1]))
  print(cor.test(param_gt[,5], param_recov[,2]))
}



plot(param_gt[,4], param_recov[,1])
text(param_gt[,4], param_recov[,1], labels = round(param_gt[,1],1), pos = 3, offset = 0.5, col = "red")
abline(a = 0, b = 1, col = "grey50")

cor.test(param_gt[,4], param_recov[,1])


plot(param_gt[,5], param_recov[,2])
text(param_gt[,5], param_recov[,2], labels = round(param_gt[,1],1), pos = 3, offset = 0.5, col = "red")
abline(a = 0, b = 1, col = "grey50")
cor.test(param_gt[,5], param_recov[,2])

cor.test(param_gt[,6], param_recov[,3])

