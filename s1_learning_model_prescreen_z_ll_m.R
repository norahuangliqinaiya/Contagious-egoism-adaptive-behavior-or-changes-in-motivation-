# The current script corresponds to Step 1 (the pre-screening process) of 
# 2 The simulation protocol
# 2.1 Establishing the egoism learning task

rm(list=ls())

library(maxLik)
library(foreach)
library(doParallel)

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


# Function to transfer degree into radian ----

degree_to_radian <- function(d){
  
  r <- d * pi/180
  
  return(r)
  
}

# Set path ----

path <- 'D:/postdoc/work/projects/egoism/experiment/data/R/learning_phase/'
# path <- 'D:/data/R/learning_phase/'


# Step 1: Pre-screen process ----

sim_gn_n <- 100
tn <- 45
subj_num <- 30
bias <- 2

sim_time <- 200

# Set up the task parameter space
lamda_rng <- seq(0,10,1)
dc_rng <- seq(0,10,0.5)
tp <- cbind(lamda = rep(lamda_rng, each = length(dc_rng)), dc = rep(dc_rng, time = length(lamda_rng)))

numCores <- detectCores()
registerDoParallel(numCores)

foreach (st = 1:sim_time, .combine = rbind) %dopar% {
  
  mn <- 1
  Loglik_all <- matrix(0, length(lamda_rng) * length(dc_rng) * sim_gn_n, 5)
  
  for (n in 1:dim(tp)[1]){
    
    for (i in 1:sim_gn_n){
      
      # -------------------------
      # Generate allocation offer
      #--------------------------
      
      phi_gt <- matrix(0, sim_gn_n, 3)
      
      set.seed(1000+i)
      phi_gt_1 <- sample(0:20, size = 1, replace = T)
      
      set.seed(2000+i)
      phi_gt_2 <- sample(20:40, size = 1, replace = T)
      
      set.seed(3000+i)
      phi_gt_3 <- sample(40:60, size = 1, replace = T)
      
      phi_gt[i, 1] <- phi_gt_1
      phi_gt[i, 2] <- phi_gt_2
      phi_gt[i, 3] <- phi_gt_3
      
      set.seed(1000+i)
      allocate_offer_1 <- as.matrix(phi_gt_1 + rpois(tn, lambda = tp[n,1]), tn, 1)
      
      set.seed(2000+i)
      allocate_offer_2 <- as.matrix(phi_gt_2 + rpois(tn, lambda = tp[n,1]), tn, 1)
      
      set.seed(3000+i)
      allocate_offer_3 <- as.matrix(phi_gt_3 + rpois(tn, lambda = tp[n,1]), tn, 1)
      
      allocate_offer <- cbind(allocate_offer_1, allocate_offer_2, allocate_offer_3)
      
      # -----------------
      # Simulate Feedback
      # -----------------
      
      FB <- matrix(0, tn, 3)
      
      for (t in 1:tn){
        
        if ((allocate_offer[t,1] >= phi_gt_1 - tp[n,2]) && (allocate_offer[t,1] <= phi_gt_1 + tp[n,2])){
          FB[t, 1] <- 1
        } else {
          FB[t, 1] <- 0
        }
        
        if ((allocate_offer[t,2] >= phi_gt_2 - tp[n,2]) && (allocate_offer[t,2] <= phi_gt_2 + tp[n,2])){
          FB[t, 2] <- 1
        } else {
          FB[t, 2] <- 0
        }
        
        if ((allocate_offer[t,3] >= phi_gt_3 - tp[n,2]) && (allocate_offer[t,3] <= phi_gt_3 + tp[n,2])){
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
      
      num_zeros <- sum(FB_new == 0)
      Loglik_all[mn,5] <- num_zeros
      
      # ------------------
      # Simulate behaviors
      # ------------------
      
      phi_gn <- matrix(0, 1, tn)
      phi_gn[1] <- 60
      
      param_est <- matrix(0, subj_num, 4)
      param_gt <- matrix(0, subj_num, 3)
      Loglik_sum <- matrix(0, 1, subj_num)
      
      for (subj in 1:subj_num){
        
        seed <- 100 * st + subj
        set.seed(seed)
        alpha_gn <- runif(1, min = 0, max = 1)
        beta_gn <- runif(1, min = 0, max = 1)
        beta_0 <- sample(seq(-bias,bias,1), 1, replace = T)
        
        param_gt[subj, 1] <- alpha_gn
        param_gt[subj, 2] <- beta_gn
        param_gt[subj, 3] <- beta_0
        
        behav_est <- matrix(0, 1, tn)
        
        for (t in 1:tn){
          
          phi_est <- matrix(0, 1, 4)
          phi_est[1] <- phi_gn[t]
          
          for (tau in 1:3){
            
            phi_est[tau+1] <- phi_est[tau] + alpha_gn * FB[t, tau] * (allocate_offer[t,tau] - phi_est[tau])
            
          }
          
          phi_gn[t+1] <- phi_est[4]
          behav_est[t] <- round(est_constrain(beta_gn * phi_gn[t+1] + beta_0), 2)
          
        }
        
        # --------------------------
        # likelihood of ground truth
        # --------------------------
        
        phi_gn_mf <- matrix(0, 1, tn+1)
        Loglik <- matrix(0, 1, tn)
        
        phi_gn_mf[1] <- 60
        
        for (j in 1:tn){
          
          phi_est_mf <- matrix(0, 1, 4)
          phi_est_mf[1] <- phi_gn_mf[j]
          
          for (tau in 1:3){
            
            phi_est_mf[tau+1] <- phi_est_mf[tau] + alpha_gn * FB[j,tau] * (allocate_offer[j,tau] - phi_est_mf[tau])
            
          }
          
          phi_gn_mf[j+1] <- phi_est_mf[4]
        }
        
        for (j in 1:tn){
          
          u <- behav_est[j] - beta_gn * phi_gn_mf[j+1] - beta_0
          Loglik[j] <- log10(dnorm(u))
        }
        
        Loglik_sum[subj] <- sum(Loglik)
        
      }
      
      Loglik_all[mn,1] <- tp[n,1]
      Loglik_all[mn,2] <- tp[n,2]
      Loglik_all[mn,3] <- i
      Loglik_all[mn,4] <- sum(Loglik_sum)
      
      mn <- mn + 1
      
    }
  }
  
  colnames(Loglik_all) <- c("lamda", "dc", "gn", "ll", 'zero_num')
  Loglik_all <- as.data.frame(Loglik_all)
  
  z_ll_st <- (Loglik_all$ll - mean(Loglik_all$ll))/sd(Loglik_all$ll)
  fn_z_ll_st <- paste0(path, 'lik_all_with_bias/z_ll_with_bias_',bias,'_st_',st,'_.txt')
  write.table(z_ll_st, fn_z_ll_st, row.names=F)

}



# Read data ----
rm(list=ls())

expand_matrix <- function(orig_matrix, n, mode){
  
  if (mode == 1){
    expand_matrix <- orig_matrix[rep(seq_len(nrow(orig_matrix)),each = n),]
  } else if (mode == 2){
    expand_matrix <- orig_matrix[rep(seq_len(nrow(orig_matrix)),n),]
  }
  
}

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


# Function to transfer degree into radian ----

degree_to_radian <- function(d){
  
  r <- d * pi/180
  
  return(r)
  
}


sim_gn_n <- 100
tn <- 45
subj_num <- 30
bias <- 2

sim_time <- 200

# Set up the task parameter space
lamda_rng <- seq(0,10,1)
dc_rng <- seq(0,10,0.5)

tp <- cbind(lamda = rep(lamda_rng, each = length(dc_rng)), dc = rep(dc_rng, time = length(lamda_rng)))
phi_gt <- matrix(0, sim_gn_n,4)
mn <- 1
fb_num_zeros <- matrix(0, length(lamda_rng) * length(dc_rng) * sim_gn_n, 1)

for (n in 1:dim(tp)[1]){

for (i in 1:sim_gn_n){
  
  set.seed(1000+i)
  phi_gt_1 <- sample(0:20, size = 1, replace = T)
  
  set.seed(2000+i)
  phi_gt_2 <- sample(20:40, size = 1, replace = T)
  
  set.seed(3000+i)
  phi_gt_3 <- sample(40:60, size = 1, replace = T)
  
  phi_gt[i,1] <- i
  phi_gt[i,2] <- phi_gt_1
  phi_gt[i,3] <- phi_gt_2
  phi_gt[i,4] <- phi_gt_3
  
  set.seed(1000+i)
  allocate_offer_1 <- as.matrix(phi_gt_1 + rpois(tn, lambda = tp[n,1]), tn, 1)
  
  set.seed(2000+i)
  allocate_offer_2 <- as.matrix(phi_gt_2 + rpois(tn, lambda = tp[n,1]), tn, 1)
  
  set.seed(3000+i)
  allocate_offer_3 <- as.matrix(phi_gt_3 + rpois(tn, lambda = tp[n,1]), tn, 1)
  
  allocate_offer <- cbind(allocate_offer_1, allocate_offer_2, allocate_offer_3)
  
  # -----------------
  # Simulate Feedback
  # -----------------
  
  FB <- matrix(0, tn, 3)
  
  for (t in 1:tn){
    
    if ((allocate_offer[t,1] >= phi_gt_1 - tp[n,2]) && (allocate_offer[t,1] <= phi_gt_1 + tp[n,2])){
      FB[t, 1] <- 1
    } else {
      FB[t, 1] <- 0
    }
    
    if ((allocate_offer[t,2] >= phi_gt_2 - tp[n,2]) && (allocate_offer[t,2] <= phi_gt_2 + tp[n,2])){
      FB[t, 2] <- 1
    } else {
      FB[t, 2] <- 0
    }
    
    if ((allocate_offer[t,3] >= phi_gt_3 - tp[n,2]) && (allocate_offer[t,3] <= phi_gt_3 + tp[n,2])){
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
  
  num_zeros <- sum(FB_new == 0)
  fb_num_zeros[mn,1] <- num_zeros
  mn <- mn + 1
  
  }
}

tp_exp <- expand_matrix(tp, n = sim_gn_n, mode = 1)
phi_gt_exp <- expand_matrix(phi_gt, n = dim(tp)[1], mode = 2)
tp_all <- cbind(tp_exp, phi_gt_exp)

colnames(tp_all) <- c("lamda", "dc", "gn", "phi_gt_1", "phi_gt_2", "phi_gt_3")

setwd('D:/postdoc/work/UKW/projects/egoism/experiment/data/R/learning_phase/lik_all_with_bias')
# setwd('D:/data/R/learning_phase/lik_all_with_bias')
z_ll <- matrix(0, length(lamda_rng) * length(dc_rng) * sim_gn_n, sim_time)

for (i in 1:sim_time){
  
  fn <- paste0('z_ll_with_bias_',bias,'_st_', i, '_', '.txt')
  z_ll[,i] <- as.matrix(read.table(fn, header = T))
  
}

z_ll <- as.data.frame(z_ll)
z_ll_m <- rowMeans(z_ll)

data_all <- cbind(tp_all, z_ll, fb_num_zeros, z_ll_m)
data <- data_all[order(data_all$z_ll_m, decreasing = T),]
data_new <- data[c("lamda", "dc", "gn", "phi_gt_1", "phi_gt_2", "phi_gt_3","fb_num_zeros", "z_ll_m")]



# Search for task parameter that meet the criteria of pre-screening ----

# Criteria 1: the likelihood z value rank higher than 0.85 of the distribution
thresh <- quantile(data_new$z_ll_m, probs = 0.85)
# Criteria 2:  the fictitious others disagree the suggested allocation plan for more than 5% and less than 50% of trials
result <-  subset(data_new, z_ll_m >= thresh & fb_num_zeros <= 45*3*0.5 & fb_num_zeros >= 45*3*0.05)

res_sub <- result[c("lamda", "dc", "phi_gt_1", "phi_gt_2", "phi_gt_3")]

result_u <- result[!duplicated(res_sub),]

path <- 'D:/postdoc/UKW/work/projects/egoism/experiment/data/R/learning_phase/'
# path <- 'D:/data/R/learning_phase/'
fn_lik_thres<- paste0(path, 'lik_all_with_bias_',bias,'_thres_0.85_zero_fb_5_50_z_ll_m.txt')
write.table(result_u, fn_lik_thres, row.names=F)



# Save stimulus for tp in result ----
path <- 'D:/postdoc/UKW/work/projects/egoism/experiment/data/R/learning_phase/z_ll_m/'
# path <- 'D:/data/R/learning_phase/z_ll_m/'

tn <- 45
subj_num <- 30
bias <- 0

for (i in 1:dim(result_u)[1]){
  
  lamda <- result_u[i,1]
  dc <- result_u[i,2]
  
  gn <- result_u[i,3]
  phi_gt_1 <- result_u[i,4]
  phi_gt_2 <- result_u[i,5]
  phi_gt_3 <- result_u[i,6]
  
  # ==========================
  # Generatve allocation offer
  # ==========================
  
  set.seed(1000+i)
  allocate_offer_1 <- as.matrix(phi_gt_1 + rpois(tn, lamda), tn, 1)
  
  set.seed(2000+i)
  allocate_offer_2 <- as.matrix(phi_gt_2 + rpois(tn, lamda), tn, 1)
  
  set.seed(3000+i)
  allocate_offer_3 <- as.matrix(phi_gt_3 + rpois(tn, lamda), tn, 1)
  
  allocate_offer <- cbind(allocate_offer_1, allocate_offer_2, allocate_offer_3)
  
  # ==================
  # Simulate Feedbacks
  # ==================
  
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
  
  
  # ==================
  # Simulate behaviors
  # ==================
  
  phi_gn <- matrix(0, 1, tn)
  phi_gn[1] <- 60
  
  param_est <- matrix(0, subj_num, 4)
  param_gt <- matrix(0, subj_num, 3)
  behav_est <- matrix(0, subj_num, tn)
  
  for (subj in 1:subj_num){
    
    set.seed(subj)
    alpha_gn <- runif(1, min = 0, max = 1)
    beta_gn <- runif(1, min = 0, max = 1)
    beta_0 <- sample(seq(-bias,bias,1), 1, replace = T)
    
    param_gt[subj, 1] <- alpha_gn
    param_gt[subj, 2] <- beta_gn
    param_gt[subj, 3] <- beta_0
    
    
    for (t in 1:tn){
      
      phi_est <- matrix(0, 1, 4)
      phi_est[1] <- phi_gn[t]
      
      for (tau in 1:3){
        
        phi_est[tau+1] <- phi_est[tau] + alpha_gn * FB[t, tau] * (allocate_offer[t,tau] - phi_est[tau])
        
      }
      
      phi_gn[t+1] <- phi_est[4]
      behav_est[subj,t] <- round(est_constrain(beta_gn * phi_gn[t+1] + beta_0), 2)
      
    }
  }
  
  
  fn_fb <- paste0(path, 'FB/lamda_', lamda, '_dc_', dc, '_gn_', gn, '.txt')
  write.table(FB_new, fn_fb, row.names = F, col.names = F)
  
  fn_offer <- paste0(path, 'allocate_offer/lamda_', lamda, '_dc_', dc, '_gn_', gn, '.txt')
  write.table(allocate_offer_new, fn_offer, row.names = F, col.names = F)
  
  fn_behav_est <- paste0(path, 'behav_est/lamda_', lamda, '_dc_', dc, '_gn_', gn, '.txt')
  write.table(behav_est, fn_behav_est, row.names = F, col.names = F)
  
}


fn_param_gt <- paste0(path,'gt_param/param_gt.txt')
write.table(param_gt, fn_param_gt, row.names = F, col.names = F)

