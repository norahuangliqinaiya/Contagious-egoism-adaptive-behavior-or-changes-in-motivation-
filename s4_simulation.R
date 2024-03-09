rm(list=ls())

library(maxLik)

# =====================
# Function of constrain
# =====================

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


allocate_constrain <- function(x){
  
  if (x > 5){
    u <- 5
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


# ===================================
# Model fitting of the learning phase
# ===================================

function_est <- function(param){
  
  alpha_gn_mf <- param[1]
  beta_gn_mf <- param[2]
  beta_0_mf <- param[3]
  
  phi_gn_mf <- matrix(0, 1, tn+1)
  Loglik <- matrix(0, 1, tn)
  
  phi_gn_mf[1] <- 60
  
  for (j in 1:tn){
    
    phi_est_mf <- matrix(0, 1, 4)
    phi_est_mf[1] <- phi_gn_mf[j]
    
    for (tau in 1:3){
      
      phi_est_mf[tau+1] <- phi_est_mf[tau] + alpha_gn_mf * FB_mf[(j-1)*3 + tau] * (allocate_offer[j,tau] - phi_est_mf[tau])
      
    }
    
    phi_gn_mf[j+1] <- phi_est_mf[4]
  }
  
  for (j in 1:tn){
    
    u <- behav_est_mf[j] - beta_gn_mf * phi_gn_mf[j+1] - beta_0_mf
    Loglik[j] <- pnorm(u, log = T)
  }
  
  return(sum(Loglik))
}


# ================
# Model fit of M1
# ================

function_allocate_M1 <- function(param){
  
  alpha_w_mf <- param[1]
  beta_self_mf <- param[2]
  beta_1_mf <- param[3]
  
  phi_self_mf <- matrix(0, 1, tn+1)
  phi_self_mf[1] <- phi_self_base_mf
  
  Loglik <- matrix(0, 1, tn)
  
  for (j in 1:tn){
    
    phi_self_mf[j+1] <- phi_self_mf[j] + alpha_w_mf * (phi_gn_mf_M1[j+1] - phi_self_mf[j])
    
  }
  
  for (j in 1:tn){

    u <- behav_allocate_mf[j]/(radius - behav_allocate_mf[j]) - beta_self_mf * (1 - cos(degree_to_radian(phi_self_mf[j+1])))/cos(degree_to_radian(phi_self_mf[j+1])) - beta_1_mf
    Loglik[j] <- pnorm(u, log = T)
  }
  
  return(sum(Loglik))
  
}

# ================
# Model fit of M2
# ===============

function_allocate_M2 <- function(param){
  
  alpha_w_mf <- param[1]
  beta_self_mf <- param[2]
  beta_1_mf <- param[3]
  
  Loglik <- matrix(0, 1, tn)
  phi_comb_mf <- matrix(0, 1, tn)
  
  
  for (j in 1:tn){
    
    phi_comb_mf[j] <- (1 - alpha_w_mf) * phi_self_base_mf + alpha_w_mf * phi_gn_mf_M2[j+1]
    
  }
  
  for (j in 1:tn){
    
    u <- behav_allocate_mf[j]/(radius - behav_allocate_mf[j]) - beta_self_mf * (1 - cos(degree_to_radian(phi_comb_mf[j])))/cos(degree_to_radian(phi_comb_mf[j])) - beta_1_mf
    Loglik[j] <- pnorm(u, log = T)
    
  }
  return(sum(Loglik))
}


# ====================
# array multiplication
# ====================

element_wise_matrix_multiply <- function(array) {
  
  num_matrices <- dim(array)[3]
  result <- matrix(1, dim(array)[1], dim(array)[2])
  
  for (i in 1:num_matrices){
    result <- result * array[, , i]
  }
  
  return(result)
}

# =================
# Set path
# =================

path <- 'F:/Research/Postdoc/work/projects/egoism/experiment/data/simulation/'

# ==================
# Main program
# ==================

subj_id <- 1
 
radius <- 5
sim_gn_n <- 2
tn <- 45

agn <- c(0.1, 0.3, 0.5, 0.9)
aw <- c(0.1, 0.3, 0.5, 0.9)
bgn <- c(0.1, 0.3, 0.5, 0.9)
bs <- c(0.1, 0.3, 0.5, 0.9)
psb <- c(15, 30, 45, 60, 75, 90)


lamda_rng <- c() # inputs from 'param_recov_learning_model'
dc_rng <- c() # inputs from 'param_recov_learning_model'
seed_i <- c() # inputs from 'param_recov_learning_model'
tp <- data.frame(lamda = lamda_rng, dc = dc_rng, i = seed_i)


lik_diff_maps <- array(0, dim = c(length(dc), sim_gn_n, length(agn) * length(aw) * length(bgn) * length(bs) * length(psb)))

mn <- 1

dp <- cbind(alpha_gn = rep(agn, each = length(aw) * length(bgn) * length(bs) * length(psb)),
            alpha_w = rep(rep(aw, each = length(bgn) * length(bs) * length(psb)), time = length(agn)),
            beta_gn = rep(rep(bgn, each = length(bs) * length(psb)), time = length(aw) * length(agn)),
            beta_self = rep(rep(bs, each = length(psb)), time = length(bgn) * length(aw) * length(agn)),
            phi_self_base = rep(psb, time = length(bs) * length(bgn) * length(aw) * length(agn)))


for (n in 1:dim(dp)[1]){
  
  cc <- 1
  likelihood <- matrix(0, length(dc), sim_gn_n)
  
  for (r in 1:dim(tp)[1]){
        
        # ==================
        # Simulate Feedbacks
        # ==================
        
        phi_gt <- matrix(0, sim_gn_n, 3)
        
        behav_est <- matrix(0, sim_gn_n, tn)
        behav_allocate <- matrix(0, sim_gn_n, tn)
        param_est <- matrix(0, sim_gn_n, 4)
        
        set.seed(1000+tp[r,3])
        phi_gt_1 <- sample(0:60, size = 1, replace = T)
        
        set.seed(2000+tp[r,3])
        phi_gt_2 <- sample(0:60, size = 1, replace = T)
        
        set.seed(3000+tp[r,3])
        phi_gt_3 <- sample(0:60, size = 1, replace = T)
            
        phi_gt[i, 1] <- phi_gt_1
        phi_gt[i, 2] <- phi_gt_2
        phi_gt[i, 3] <- phi_gt_3
              
        
        # =========================
        # Generate allocation offer
        # =========================
              
        set.seed(1000+tp[r,3])
        allocate_offer_1 <- as.matrix(phi_gt_1 + rpois(tn, lambda = tp[r,1]), tn, 1)
        
        set.seed(2000+tp[r,3])
        allocate_offer_2 <- as.matrix(phi_gt_2 + rpois(tn, lambda = tp[r,1]), tn, 1)
              
        set.seed(3000+tp[r,3])
        allocate_offer_3 <- as.matrix(phi_gt_3 + rpois(tn, lambda = tp[r,1]), tn, 1)
              
        allocate_offer <- cbind(allocate_offer_1, allocate_offer_2, allocate_offer_3)
              
        phi_gn <- matrix(0, 1, tn)
        phi_gn[1] <- 60
              
        phi_self <- matrix(0, 1, tn)
        phi_self[1] <- phi_self_base
              
              
        FB <- matrix(0, tn, 3)
        phi_comb <- matrix(0, 1, tn)
              
              
        for (t in 1:tn){
          
          if ((allocate_offer[t,1] >= phi_gt_1 - tp[r,2]) && (allocate_offer[t,1] <= phi_gt_1 + tp[r,2])){
            FB[t, 1] <- 1
            } else {
              FB[t, 1] <- 0
              }
                
          if ((allocate_offer[t,2] >= phi_gt_2 - tp[r,2]) && (allocate_offer[t,2] <= phi_gt_2 + tp[r,2])){
            FB[t, 2] <- 1
            } else {
              FB[t, 2] <- 0
              }
                
          if ((allocate_offer[t,3] >= phi_gt_3 - tp[r,2]) && (allocate_offer[t,3] <= phi_gt_3 + tp[r,2])){
            FB[t, 3] <- 1
            } else {
              FB[t, 3] <- 0
              }
                
          # ==============================
          # Simulate behaviors based on M2
          # ==============================
                
          phi_est <- matrix(0, 1, 4)
          phi_est[1] <- phi_gn[t]
                
            for (tau in 1:3){
                  
                  phi_est[tau+1] <- phi_est[tau] + alpha_gn * FB[t, tau] * (allocate_offer[t,tau] - phi_est[tau])
                  
                }
                
                phi_gn[t+1] <- phi_est[4]
                
                set.seed(subj_id)
                behav_est[i, t] <- round(est_constrain(beta_gn * phi_gn[t+1] + runif(1, min = 0, max = 0)), 2)
                
                
                # phi_self[t+1] <- phi_self[t] + alpha_w * (phi_gn[t+1] - phi_self[t])
                # behav_allocate[i, t] <- round(allocate_constrain(beta_self * radius * (1 - cos(degree_to_radian(phi_self[t+1]))) + runif(1, min = 0, max = 5)), 2)
                
                phi_comb[t] <- (1 - alpha_w) * phi_self_base + alpha_w * phi_gn[t+1]
                behav_allocate[r, t] <- round(allocate_constrain(beta_self * radius * (1 - cos(degree_to_radian(phi_comb[t]))) + runif(1, min = 0, max = 0)), 2)
                
                subj_id <- subj_id + 1
                
              }
              
              # ===================================
              # Model fit of the learning phase
              # ===================================
              
              # Input data
              behav_est_mf <- behav_est[i, 1:tn]
              FB_mf <- FB
              init_num <- 500
              
              result <- matrix(0,init_num,5)
              
              for (m in 1:init_num){
                
                sucess <- FALSE
                
                while (!success){
                  
                  tryCatch({
                    
                    alpha_gn_mf <- as.numeric(sample(1:9.9, size = 1, replace = T)/10)
                    beta_gn_mf <- as.numeric(sample(1:9.9, size = 1, replace = T)/10)
                    beta_0_mf <- as.numeric(sample(-99.9:99.9, size = 1, replace = T)/10)
                    
                    param <- c(alpha_gn_mf, beta_gn_mf, beta_0_mf)
                    
                    A <- matrix(c(1,0,0,
                                  -1,0,0,
                                  0,1,0,
                                  0,-1,0,
                                  0,0,1,
                                  0,0,-1), 6, 3, byrow = T)
                    
                    B <- matrix(c(0, 1, 0, 1, 10, 10), 6, 1, byrow = T)
                    
                    model.fit <- maxLik(function_est, start = param, constraint = list(ineqA=A, ineqB=B), method = 'NM')
                    model.summary <- summary(model.fit)
                    
                    data <- cbind(m, as.numeric(coef(model.summary)[1,1]), as.numeric(coef(model.summary)[2,1]),
                                  as.numeric(coef(model.summary)[3,1]), as.numeric(summary(model.fit)[5]))
                    
                    result[m,] <- data
                    
                    success <- TRUE
                  
                  }, error = function(e) {
                  
                    cat("Error occurred in optimization with initial value", m, "\n")
                  })
                }
              }
              
              alpha_gn_mf <- result[(which.max(result[,5])),2]
              beta_gn_mf <- result[(which.max(result[,5])),3]
              beta_0_mf <- result[(which.max(result[,5])),4]
              lh <- result[(which.max(result[,5])),5]
              
              param_est[i,] <- cbind(alpha_gn_mf, beta_gn_mf, beta_0_mf, lh)
              
              
              # Calculate the phi_gn
              phi_gn_mf <- matrix(0, 1, tn+1)
              phi_gn_mf[1] <- 60
              
              for (j in 1:tn){
                
                phi_est_mf <- matrix(0, 1, 4)
                phi_est_mf[1] <- phi_gn_mf[j]
                
                for (tau in 1:3){
                  
                  phi_est_mf[tau+1] <- phi_est_mf[tau] + alpha_gn_mf * FB_mf[(j-1)*3 + tau] * (allocate_offer[j, tau] - phi_est_mf[tau])
                  
                }
                
                phi_gn_mf[j+1] <- phi_est_mf[4]
                
              }
              
              # ================
              # Model fit of M1
              # ================
              
              # Input data
              phi_self_base_mf <- phi_self_base
              phi_gn_mf_M1 <- phi_gn_mf
              behav_allocate_mf <- behav_allocate[i, 1:tn]
              
              result <- matrix(0,init_num,5)
              likelihood_M1 <- matrix(0, length(dc), sim_gn_n)
              
              for (m in 1:init_num){
                
                sucess <- FALSE
                
                while(!sucess){
                  
                  tryCatch({
                    
                    alpha_w_mf <- as.numeric(sample(1:9.9, size = 1, replace = T)/10)
                    beta_self_mf <- as.numeric(sample(1:9.9, size = 1, replace = T)/10)
                    beta_1_mf <- as.numeric(sample(1:49.9, size = 1, replace = T)/10)
                    
                    param <- c(alpha_w_mf, beta_self_mf, beta_1_mf)
                    
                    A <- matrix(c(1,0,0,
                                  -1,0,0,
                                  0,1,0,
                                  0,-1,0,
                                  0,0,1,
                                  0,0,-1), 6, 3, byrow = T)
                    
                    B <- matrix(c(0, 1, 0, 1, 0, 5), 6, 1, byrow = T)
                    
                    model.fit <- maxLik(function_allocate_M1, start = param, constraint = list(ineqA=A, ineqB=B), method = 'NM')
                    
                    success <- TRUE
                    
                  }, error = function(e) {
                    
                    cat("Error occurred in optimization with initial value", m, "\n")
                    
                  })
                }
                
                likelihood_M1[cc,i] <- as.numeric(summary(model.fit)[5])
              }
              
              # ================
              # Model fit of M2
              # ================
              
              # Input data
              phi_gn_mf_M2 <- phi_gn_mf
              behav_allocate_mf <- behav_allocate[i, 1:tn]
              
              result <- matrix(0,init_num,5)
              likelihood_M2 <- matrix(0, length(dc), sim_gn_n)
              
              for (m in 1:init_num){
                
                success <- FALSE
                
                while(!success){
                  
                  tryCatch({
                    
                    alpha_w_mf <- as.numeric(sample(1:9.9, size = 1, replace = T)/10)
                    beta_self_mf <- as.numeric(sample(1:9.9, size = 1, replace = T)/10)
                    beta_1_mf <- as.numeric(sample(1:49.9, size = 1, replace = T)/10)
                    
                    param <- c(alpha_w_mf, beta_self_mf, beta_1_mf)
                    
                    A <- matrix(c(1,0,0,
                                  -1,0,0,
                                  0,1,0,
                                  0,-1,0,
                                  0,0,1,
                                  0,0,-1), 6, 3, byrow = T)
                    
                    B <- matrix(c(0, 1, 0, 1, 0, 5), 6, 1, byrow = T)
                    
                    model.fit <- maxLik(function_allocate_M2, start = param, constraint = list(ineqA=A, ineqB=B), method = 'NM')
                    
                  }, error = function(e) {
                    
                    cat("Error occurred in optimization with initial value", m, "\n")
                    
                  })
                }
                likelihood_M2[cc,i] <- as.numeric(summary(model.fit)[5])
              }
            }
            
            cc <- cc + 1
            
            # ==========================
            # Output simulated behaviors
            # ==========================
            
            # fn_est <- paste0(path, 'behav_est_', alpha_gn, '_', alpha_w, '_', beta_gn, '_', beta_self, '_', c, '.csv')
            # write.csv(behav_est, fn_est, row.names=F)
            # 
            # fn_allocate <- paste0(path, 'behav_allocate_', alpha_gn, '_', alpha_w, '_', beta_gn, '_', beta_self, '_', c, '.csv')
            # write.csv(behav_allocate, fn_allocate, row.names=F)
            # 
            # fn_gn_gt <- paste0(path, 'Group_norm_ground_truth.csv')
            # write.csv(phi_gt, fn_gn_gt, row.names=F)
            # 
            # fn_param_est <- paste0(path, 'param_est_', alpha_gn, '_', beta_gn,  '_', c, '.csv')
            # colnames(param_est) <- c('alpha_gn', 'beta_gn', 'beta_0', 'lh')
            # write.csv(param_est, fn_param_est, row.names=F)

        likelihood_diff <- abs(likelihood_M1 - likelihood_M2)
        
        # fn_lik <- paste0(path, 'Loglik_', alpha_gn, '_', alpha_w, '_', beta_gn, '_', beta_self, '.csv')
        # write.csv(likelihood, fn_lik, row.names = F)

        likelihood_diff[likelihood_diff <= quantile(likelihood_diff, probs = 0.25)] <- 0
        
        lik_diff_maps[,,mn] <- likelihood_diff
        mn <- mn + 1
      }


lik_map_final <- element_wise_matrix_multiply(lik_diff_maps)
