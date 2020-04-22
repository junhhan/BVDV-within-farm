rm(list= ls())

setwd("D:/Temp/Within")

library(parallel)
library(pscl)
library(MASS)
library(lmtest)
library(ggplot2)
library(grid)
library(gridExtra)

# Setting the parallel computing
detectCores(logical= T)
n_core <- 10 #5
cl <- makeCluster(n_core)

###################################################
# Loss due to BVDV (4% Trojan + 2% Neighbouring PIs)
###################################################
# Setting the variables required for the model
lims <- 1000/n_core
n_year <- 10

list_ma <- c(50, 150, 300)
list_beta <- c(0.11, 0.05, 0.4, 0.4)
list_rho <- c(0, 0.04)
list_nb <- c(0, 0.02)
list_ctrl <- c(0, 0, 0)

result_ti <- NULL
result_tti <- NULL
result_pi <- NULL
result_tpi <- NULL
result_wgt <- NULL
result_twgt <- NULL
result_sup <- NULL
result_tsup <- NULL
result_pur <- NULL
result_tpur <- NULL
result_gfr <- NULL
result_tgfr <- NULL
result_ctest <- NULL
result_tctest <- NULL
result_cvac <- NULL
result_tcvac <- NULL
result_cdf <- NULL
result_tcdf <- NULL
result_cost <- NULL
result_tcost <- NULL

clusterExport(cl, c("lims", "n_year", "list_ma", "list_beta", "list_rho", "list_nb", "list_ctrl",
                    "result_ti", "result_pi", "result_wgt", "result_sup", "result_pur", "result_gfr", 
                    "result_ctest", "result_cvac", "result_cdf", "result_cost",
                    "result_tti", "result_tpi", "result_twgt", "result_tsup", "result_tpur", 
                    "result_tgfr", "result_tctest", "result_tcvac", "result_tcdf", "result_tcost"))

result_ma <- NULL; result_rho <- NULL; result_nb <- NULL
result_wp <- NULL; result_wt <- NULL; result_bh <- NULL; result_bf <- NULL;
clusterExport(cl, c("result_ma", "result_rho", "result_nb", "result_wp", "result_wt", "result_bh", "result_bf"))

n_ti <- n_pi <- c_test <- c_vac <- c_df <- rep(0, n_year)
wgt <- sup <- pur <- rep(0, n_year) 

clusterExport(cl, c("n_ti", "n_pi", "c_test", "c_vac", "c_df", "wgt", "sup", "pur"))

# Setting the BVDV simulation fuction
bvdsim <- function(seed, n_ma, betas, rho, SCR, VAC, DF) {
  
  .C("outbreak_beef", as.integer(seed), as.integer(n_ma), as.double(betas), as.double(rho), as.double(p_neighbour),
     as.integer(SCR), as.integer(VAC), as.integer(DF),
     as.integer(n_ti), as.integer(n_pi), as.double(wgt), as.double(sup), as.double(pur),
     as.double(c_test), as.double(c_vac), as.double(c_df))
}

clusterExport(cl, c('bvdsim'))

# Set the different seeds for the core
clusterApply(cl, seq_along(cl), function(i) seed <<- (i-1)*lims)

# Test whether BVDV is endemic and estimate the economic loss
t_init <- proc.time()
result <- 
  clusterEvalQ(cl, {
    
    # Loading the compiled C code
    dyn.load("D:/Temp/Within/outbreak_beef.dll")
    
    for (i in 2:2) {
      for (j in 1:2) {
        for (k in 1:2) {
          n_ma <- list_ma[i]
          betas <- c(list_beta[1], list_beta[2], list_beta[3], list_beta[4])
          rho <- list_rho[j]
          p_neighbour <- list_nb[k]
          
          iter <- 0
          while (iter < lims) {
            result <- bvdsim(seed, n_ma, betas, rho, list_ctrl[1], list_ctrl[2], list_ctrl[3])
            
            q <- 9
            n_ti <- result[[q]]
            n_pi <- result[[q+1]]
            wgt <- result[[q+2]]
            sup <- result[[q+3]]
            pur <- result[[q+4]]
            gfr <- wgt + sup - pur
            
            test <- result[[q+5]]
            vac <- result[[q+6]]
            df <- result[[q+7]]
            cost <- test + vac + df
            
            tti <- 0; tpi <- 0;
            twgt <- 0; tsup <-0; tpur <-0; tgfr <- 0; ttest <- 0; tvac <- 0; tdf <- 0; tcost <- 0;
            period <- 4
            for (q in (n_year-period):n_year) {
              tti <- tti + n_ti[q];
              tpi <- tpi + n_pi[q];
              
              twgt <- twgt + wgt[q]/((1.019)^(q-(n_year-period)))
              tsup <- tsup + sup[q]/((1.019)^(q-(n_year-period)))
              tpur <- tpur + pur[q]/((1.019)^(q-(n_year-period)))
              tgfr <- tgfr + (wgt[q] + sup[q] - pur[q])/((1.019)^(q-(n_year-period)))
              ttest <- ttest + test[q]/((1.019)^(q-(n_year-period)))
              tvac <- tvac + vac[q]/((1.019)^(q-(n_year-period)))
              tdf <- tdf + df[q]/((1.019)^(q-(n_year-period)))
              tcost <- tcost + (test[q] + vac[q]+ df[q])/((1.019)^(q-(n_year-period)))
            }
            tti <- tti / n_ma
            tpi <- tpi / n_ma
            
            iter <- iter + 1
            seed <- seed + 1
            
            result_ti <- c(result_ti, n_ti[c((n_year-period):n_year)])
            result_tti <- c(result_tti, tti)
            result_pi <- c(result_pi, n_pi[c((n_year-period):n_year)])
            result_tpi <- c(result_tpi, tpi)
            
            result_ma <- c(result_ma, n_ma)
            result_rho <- c(result_rho, rho)
            result_nb <- c(result_nb, p_neighbour)

            result_wgt <- c(result_wgt, wgt[c((n_year-period):n_year)])
            result_twgt <- c(result_twgt, twgt)
            result_sup <- c(result_sup, sup[c((n_year-period):n_year)])
            result_tsup <- c(result_tsup, tsup)
            result_pur <- c(result_pur, pur[c((n_year-period):n_year)])
            result_tpur <- c(result_tpur, tpur)
            result_gfr <- c(result_gfr, gfr[c((n_year-period):n_year)]); 
            result_tgfr <- c(result_tgfr, tgfr); 
            result_ctest <- c(result_ctest, test[c((n_year-period):n_year)])
            result_tctest <- c(result_tctest, ttest)
            result_cvac <- c(result_cvac, vac[c((n_year-period):n_year)]); 
            result_tcvac <- c(result_tcvac, tvac); 
            result_cdf <- c(result_cdf, df[c((n_year-period):n_year)])
            result_tcdf <- c(result_tcdf, tdf)
            result_cost <- c(result_cost, cost[c((n_year-period):n_year)]);
            result_tcost <- c(result_tcost, tcost);
          }
        }
      }
    }
    list(
      result_ma, result_rho, result_nb, 
      result_ti, result_pi, result_wgt, result_sup, result_pur, result_gfr, result_ctest, result_cvac, result_cdf, result_cost,
      result_tti, result_tpi, result_twgt, result_tsup, result_tpur, result_tgfr, result_tctest, result_tcvac, result_tcdf, result_tcost)
  })
t_end <- proc.time(); t_end - t_init
stopCluster(cl)

for (i in 1:n_core) {
  result_ma <- c(result_ma, result[[i]][[1]])
  result_rho <- c(result_rho, result[[i]][[2]])
  result_nb <- c(result_nb, result[[i]][[3]])
  result_ti <- c(result_ti, result[[i]][[4]])
  result_pi <- c(result_pi, result[[i]][[5]])
  result_wgt <- c(result_wgt, result[[i]][[6]])
  result_sup <- c(result_sup, result[[i]][[7]])
  result_pur <- c(result_pur, result[[i]][[8]])
  result_gfr <- c(result_gfr, result[[i]][[9]])
  result_ctest <- c(result_ctest, result[[i]][[10]])
  result_cvac <- c(result_cvac, result[[i]][[11]])
  result_cdf <- c(result_cdf, result[[i]][[12]])
  result_cost <- c(result_cost, result[[i]][[13]])
  result_tti <- c(result_tti, result[[i]][[14]])
  result_tpi <- c(result_tpi, result[[i]][[15]])
  result_twgt <- c(result_twgt, result[[i]][[16]])
  result_tsup <- c(result_tsup, result[[i]][[17]])
  result_tpur <- c(result_tpur, result[[i]][[18]])
  result_tgfr <- c(result_tgfr, result[[i]][[19]])
  result_tctest <- c(result_tctest, result[[i]][[20]])
  result_tcvac <- c(result_tcvac, result[[i]][[21]])
  result_tcdf <- c(result_tcdf, result[[i]][[22]])
  result_tcost <- c(result_tcost, result[[i]][[23]])
}

###################################################
# Control options
###################################################
rm(list= ls())

setwd("D:/Temp/Within")

detectCores(logical= T)
n_core <- 10 #5
cl <- makeCluster(n_core)

# Setting the variables required for the model
lims <- 1000/n_core
n_year <- 10

list_ma <- c(50, 150, 300)
list_beta <- c(0.11, 0.05, 0.4, 0.4)
list_rho <- c(0, 0.04)
list_nb <- c(0, 0.02)

list_scr <- c(0, 1, 2)
list_vac <- c(0, 1, 2, 3, 4)
list_df <- c(0, 1)

result_tti <- NULL
result_tpi <- NULL
result_twgt <- NULL
result_tsup <- NULL
result_tpur <- NULL
result_tgfr <- NULL
result_tctest <- NULL
result_tcvac <- NULL
result_tcdf <- NULL
result_tcost <- NULL

clusterExport(cl, c("lims", "n_year", "list_ma", "list_beta", "list_rho", "list_nb", "list_scr", "list_vac", "list_df", 
                    "result_tti", "result_tpi", "result_twgt", "result_tsup", "result_tpur", "result_tgfr", 
                    "result_tctest", "result_tcvac", "result_tcdf", "result_tcost"))

result_ma <- NULL; result_rho <- NULL; result_scr <- NULL; result_vac <- NULL; result_df <- NULL
clusterExport(cl, c("result_ma", "result_rho", "result_scr", "result_vac", "result_df"))

n_ti <- n_pi <- c_test <- c_vac <- c_df <- rep(0, n_year)
wgt <- sup <- pur <- rep(0, n_year) 

clusterExport(cl, c("n_ti", "n_pi", "c_test", "c_vac", "c_df", "wgt", "sup", "pur"))

# Setting the BVDV simulation fuction
bvdsim <- function(seed, n_ma, betas, rho, SCR, VAC, DF) {
  
  .C("outbreak_beef", as.integer(seed), as.integer(n_ma), as.double(betas), as.double(rho), as.double(p_neighbour),
     as.integer(SCR), as.integer(VAC), as.integer(DF),
     as.integer(n_ti), as.integer(n_pi), as.double(wgt), as.double(sup), as.double(pur),  
     as.double(c_test), as.double(c_vac), as.double(c_df))
}

clusterExport(cl, c('bvdsim'))

# Set the different seeds for the core
clusterApply(cl, seq_along(cl), function(i) seed <<- (i-1)*lims)

# Test whether BVDV is endemic and estimate the economic loss
t_init <- proc.time()
result <- 
  clusterEvalQ(cl, {
    
    # Loading the compiled C code
    dyn.load("D:/Temp/Within/outbreak_beef.dll")
    
    for (i in 2:2) {
      for (j in 2:2) {
        for (k in 2:2) {
          for (l in 1:3) {
            for (m in 1:5) {
              for (n in 1:2) {
                n_ma <- list_ma[i]
                betas <- c(list_beta[1], list_beta[2], list_beta[3], list_beta[4])
                rho <- list_rho[j]
                p_neighbour <- list_nb[k]
                
                iter <- 0
                while (iter < lims) {
                  result <- bvdsim(seed, n_ma, betas, rho, list_scr[l], list_vac[m], list_df[n])
                  
                  q <- 9
                  n_ti <- result[[q]]
                  n_pi <- result[[q+1]]
                  wgt <- result[[q+2]]
                  sup <- result[[q+3]]
                  pur <- result[[q+4]]

                  test <- result[[q+5]]
                  vac <- result[[q+6]]
                  df <- result[[q+7]]

                  tti <- 0; tpi <- 0;
                  twgt <- 0; tsup <-0; tpur <-0; tgfr <- 0; ttest <- 0; tvac <- 0; tdf <- 0; tcost <- 0;
                  period <- 4
                  for (q in (n_year-period):n_year) {
                    tti <- tti + n_ti[q];
                    tpi <- tpi + n_pi[q];
                    
                    twgt <- twgt + wgt[q]/((1.019)^(q-(n_year-period)))
                    tsup <- tsup + sup[q]/((1.019)^(q-(n_year-period)))
                    tpur <- tpur + pur[q]/((1.019)^(q-(n_year-period)))
                    tgfr <- tgfr + (wgt[q] + sup[q] - pur[q])/((1.019)^(q-(n_year-period)))
                    ttest <- ttest + test[q]/((1.019)^(q-(n_year-period)))
                    tvac <- tvac + vac[q]/((1.019)^(q-(n_year-period)))
                    tdf <- tdf + df[q]/((1.019)^(q-(n_year-period)))
                    tcost <- tcost + (test[q] + vac[q]+ df[q])/((1.019)^(q-(n_year-period)))
                  }
                  tti <- tti / n_ma
                  tpi <- tpi / n_ma
                  
                  iter <- iter + 1
                  seed <- seed + 1
                  
                  result_ma <- c(result_ma, list_ma[i])
                  result_rho <- c(result_rho, list_rho[j])
                  result_scr <- c(result_scr, list_scr[l])
                  result_vac <- c(result_vac, list_vac[m])
                  result_df <- c(result_df, list_df[n])

                  result_tti <- c(result_tti, tti)
                  result_tpi <- c(result_tpi, tpi)
                  
                  result_twgt <- c(result_twgt, twgt); result_tsup <- c(result_tsup, tsup); result_tpur <- c(result_tpur, tpur); 
                  result_tgfr <- c(result_tgfr, tgfr); result_tctest <- c(result_tctest, ttest); result_tcvac <- c(result_tcvac, tvac); 
                  result_tcdf <- c(result_tcdf, tdf); result_tcost <- c(result_tcost, tcost);
                }
              }
            }      
          }
        }
      }
    }
    list(
      result_ma, result_rho, result_scr, result_vac, result_df,
      result_tti, result_tpi, result_twgt, result_tsup, result_tpur, result_tgfr, result_tctest, result_tcvac, result_tcdf, result_tcost)
  })
t_end <- proc.time(); t_end - t_init
stopCluster(cl)

for (i in 1:n_core) {
  result_ma <- c(result_ma, result[[i]][[1]])
  result_rho <- c(result_rho, result[[i]][[2]])
  result_scr <- c(result_scr, result[[i]][[3]])
  result_vac <- c(result_vac, result[[i]][[4]])
  result_df <- c(result_df, result[[i]][[5]])
  result_tti <- c(result_tti, result[[i]][[6]])
  result_tpi <- c(result_tpi, result[[i]][[7]])
  result_twgt <- c(result_twgt, result[[i]][[8]])
  result_tsup <- c(result_tsup, result[[i]][[9]])
  result_tpur <- c(result_tpur, result[[i]][[10]])
  result_tgfr <- c(result_tgfr, result[[i]][[11]])
  result_tctest <- c(result_tctest, result[[i]][[12]])
  result_tcvac <- c(result_tcvac, result[[i]][[13]])
  result_tcdf <- c(result_tcdf, result[[i]][[14]])
  result_tcost <- c(result_tcost, result[[i]][[15]])
}

dat <- data.frame("n_ma"= result_ma, "rho"= result_rho, "scr"= result_scr, "vac"= result_vac, "df"= result_df, 
                  "ti"= result_tti, "pi"= result_tpi, "wgt"= result_twgt, "sup"= result_tsup, "pur"= result_tpur, "gfr"= result_tgfr, 
                  "ctest"= result_tctest, "cvac"= result_tcvac, "cdf"= result_tcdf, "cost"= result_tcost)
dat$bal <- dat$gfr - dat$cost
