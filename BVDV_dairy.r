rm(list= ls())

setwd("D:/Temp/Within")

library(parallel)
library(pscl)
library(MASS)
library(lmtest)
library(ggplot2)
library(grid)
library(ggrepel)
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

list_ma <- c(100, 250, 400)
list_beta <- c(0.5, 0.05, 0.4, 0.4)
list_rho <- c(0, 0.04)
list_nb <- c(0, 0.02)
list_ctrl <- c(0, 0, 0)

result_ti <- NULL
result_pi <- NULL
result_ms <- NULL
result_wgt <- NULL
result_sup <- NULL
result_pur <- NULL
result_bob <- NULL
result_ai <- NULL
result_tx <- NULL
result_gfr <- NULL
result_ctest <- NULL
result_cvac <- NULL
result_cdf <- NULL
result_cost <- NULL
result_tti <- NULL
result_tpi <- NULL
result_tms <- NULL
result_twgt <- NULL
result_tsup <- NULL
result_tpur <- NULL
result_tbob <- NULL
result_tai <- NULL
result_ttx <- NULL
result_tgfr <- NULL
result_tctest <- NULL
result_tcvac <- NULL
result_tcdf <- NULL
result_tcost <- NULL

clusterExport(cl, c("lims", "n_year", "list_ma", "list_beta", "list_rho", "list_nb", "list_ctrl", 
                    "result_ti", "result_pi", "result_ms", "result_wgt", "result_sup", "result_pur", "result_bob", 
                    "result_ai", "result_tx", "result_gfr", "result_ctest", "result_cvac", "result_cdf", "result_cost",
                    "result_tti", "result_tpi", "result_tms", "result_twgt", "result_tsup", "result_tpur", "result_tbob", 
                    "result_tai", "result_ttx", "result_tgfr", "result_tctest", "result_tcvac", "result_tcdf", "result_tcost"))

result_ma <- NULL; result_rho <- NULL; result_nb <- NULL
clusterExport(cl, c("result_ma", "result_rho", "result_nb"))

n_ti <- n_pi <- c_test <- c_vac <- c_df <- rep(0, n_year)
ms <- wgt <- sup <- pur <- bob <- ai <- tx <- rep(0, n_year) 

clusterExport(cl, c("n_ti", "n_pi", "c_test", "c_vac", "c_df", "ms", "wgt", "sup", "pur", "bob", "ai", "tx"))

# Setting the BVDV simulation fuction
bvdsim <- function(seed, n_ma, betas, rho, SCR, VAC, DF) {
  
  .C("outbreak_dairy", as.integer(seed), as.integer(n_ma), as.double(betas), as.double(rho), as.double(p_neighbour),
     as.integer(SCR), as.integer(VAC), as.integer(DF),
     as.integer(n_ti), as.integer(n_pi), as.double(ms), as.double(wgt), as.double(sup), as.double(pur), as.double(bob), 
     as.double(ai), as.double(tx), as.double(c_test), as.double(c_vac), as.double(c_df))
}

clusterExport(cl, c('bvdsim'))

# Set the different seeds for the core
clusterApply(cl, seq_along(cl), function(i) seed <<- (i-1)*lims)

# Test whether BVDV is endemic and estimate the economic loss
t_init <- proc.time()
result <- 
  clusterEvalQ(cl, {
    
    # Loading the compiled C code
    dyn.load("D:/Temp/Within/outbreak_dairy.dll")
    
    for (i in 3:3) {
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
            ms <- result[[q+2]]
            wgt <- result[[q+3]]
            sup <- result[[q+4]]
            pur <- result[[q+5]]
            bob <- result[[q+6]]
            ai <- result[[q+7]]
            tx <- result[[q+8]]
            gfr <- ms + wgt + sup + bob - pur - ai - tx
            
            test <- result[[q+9]]
            vac <- result[[q+10]]
            df <- result[[q+11]]
            cost <- test + vac + df
            
            tti <- 0; tpi <- 0;
            tms <- 0; twgt <- 0; tsup <- 0; tpur <- 0; tbob <- 0; 
            tai <- 0; ttx <- 0; tgfr <- 0; ttest <- 0; tvac <- 0; tdf <- 0; tcost <- 0;
            period <- 4
            for (q in (n_year-period):n_year) {
              tti <- tti + n_ti[q];
              tpi <- tpi + n_pi[q];
              
              tms <- tms + ms[q]/((1.019)^(q-(n_year-period)))
              twgt <- twgt + wgt[q]/((1.019)^(q-(n_year-period)))
              tsup <- tsup + sup[q]/((1.019)^(q-(n_year-period)))
              tpur <- tpur + pur[q]/((1.019)^(q-(n_year-period)))
              tbob <- tbob + bob[q]/((1.019)^(q-(n_year-period)))
              tai <- tai + ai[q]/((1.019)^(q-(n_year-period)))
              ttx <- ttx + tx[q]/((1.019)^(q-(n_year-period)))
              tgfr <- tgfr + (ms[q] + wgt[q] + sup[q] + bob[q] - pur[q] - ai[q] - tx[q])/((1.019)^(q-(n_year-period)))
              ttest <- ttest + test[q]/((1.019)^(q-(n_year-period)))
              tvac <- tvac + vac[q]/((1.019)^(q-(n_year-period)))
              tdf <- tdf + df[q]/((1.019)^(q-(n_year-period)))
              tcost <- tcost + (test[q] + vac[q] + df[q])/((1.019)^(q-(n_year-period)))
            }
            tti <- tti
            tpi <- tpi
            
            iter <- iter + 1
            seed <- seed + 1
            
            result_ti <- c(result_ti, n_ti[c((n_year-period):n_year)])
            result_tti <- c(result_tti, tti)
            result_pi <- c(result_pi, n_pi[c((n_year-period):n_year)])
            result_tpi <- c(result_tpi, tpi)
            
            result_ma <- c(result_ma, n_ma)
            result_rho <- c(result_rho, rho)
            result_nb <- c(result_nb, p_neighbour)
            
            result_ms <- c(result_ms, ms[c((n_year-period):n_year)])
            result_tms <- c(result_tms, tms)
            result_wgt <- c(result_wgt, wgt[c((n_year-period):n_year)])
            result_twgt <- c(result_twgt, twgt)
            result_sup <- c(result_sup, sup[c((n_year-period):n_year)])
            result_tsup <- c(result_tsup, tsup)
            result_pur <- c(result_pur, pur[c((n_year-period):n_year)])
            result_tpur <- c(result_tpur, tpur)
            result_bob <- c(result_bob, bob[c((n_year-period):n_year)])
            result_tbob <- c(result_tbob, tbob)
            result_ai <- c(result_ai, ai[c((n_year-period):n_year)])
            result_tai <- c(result_tai, tai)
            result_tx <- c(result_tx, tx[c((n_year-period):n_year)])
            result_ttx <- c(result_ttx, ttx)
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
      result_ti, result_pi, result_ms, result_wgt, result_sup, result_pur, result_bob, result_ai, result_tx, result_gfr, 
      result_ctest, result_cvac, result_cdf, result_cost,
      result_tti, result_tpi, result_tms, result_twgt, result_tsup, result_tpur, result_tbob, result_tai, result_ttx, result_tgfr, 
      result_tctest, result_tcvac, result_tcdf, result_tcost)
  })
t_end <- proc.time(); t_end - t_init
stopCluster(cl)

for (i in 1:n_core) {
  result_ma <- c(result_ma, result[[i]][[1]])
  result_rho <- c(result_rho, result[[i]][[2]])
  result_nb <- c(result_nb, result[[i]][[3]])
  result_ti <- c(result_ti, result[[i]][[4]])
  result_pi <- c(result_pi, result[[i]][[5]])
  result_ms <- c(result_ms, result[[i]][[6]])
  result_wgt <- c(result_wgt, result[[i]][[7]])
  result_sup <- c(result_sup, result[[i]][[8]])
  result_pur <- c(result_pur, result[[i]][[9]])
  result_bob <- c(result_bob, result[[i]][[10]])
  result_ai <- c(result_ai, result[[i]][[11]])
  result_tx <- c(result_tx, result[[i]][[12]])
  result_gfr <- c(result_gfr, result[[i]][[13]])
  result_ctest <- c(result_ctest, result[[i]][[14]])
  result_cvac <- c(result_cvac, result[[i]][[15]])
  result_cdf <- c(result_cdf, result[[i]][[16]])
  result_cost <- c(result_cost, result[[i]][[17]])
  result_tti <- c(result_tti, result[[i]][[18]])
  result_tpi <- c(result_tpi, result[[i]][[19]])
  result_tms <- c(result_tms, result[[i]][[20]])
  result_twgt <- c(result_twgt, result[[i]][[21]])
  result_tsup <- c(result_tsup, result[[i]][[22]])
  result_tpur <- c(result_tpur, result[[i]][[23]])
  result_tbob <- c(result_tbob, result[[i]][[24]])
  result_tai <- c(result_tai, result[[i]][[25]])
  result_ttx <- c(result_ttx, result[[i]][[26]])
  result_tgfr <- c(result_tgfr, result[[i]][[27]])
  result_tctest <- c(result_tctest, result[[i]][[28]])
  result_tcvac <- c(result_tcvac, result[[i]][[29]])
  result_tcdf <- c(result_tcdf, result[[i]][[30]])
  result_tcost <- c(result_tcost, result[[i]][[31]])
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

list_ma <- c(100, 250, 400)
list_beta <- c(0.5, 0.05, 0.4, 0.4)
list_rho <- c(0, 0.04)
list_nb <- c(0, 0.02)

list_scr <- c(0, 1, 2)
list_vac <- c(0, 1, 2, 3, 4)
list_df <- c(0, 1)

result_tti <- NULL
result_tpi <- NULL
result_tms <- NULL
result_twgt <- NULL
result_tsup <- NULL
result_tpur <- NULL
result_tbob <- NULL
result_tai <- NULL
result_ttx <- NULL
result_tgfr <- NULL
result_tctest <- NULL
result_tcvac <- NULL
result_tcdf <- NULL
result_tcost <- NULL

clusterExport(cl, c("lims", "n_year", "list_ma", "list_beta", "list_rho", "list_nb", "list_scr", "list_vac", "list_df", 
                    "result_tti", "result_tpi", "result_tms", "result_twgt", "result_tsup", "result_tpur", "result_tbob", 
                    "result_tai", "result_ttx", "result_tgfr", 
                    "result_tctest", "result_tcvac", "result_tcdf", "result_tcost"))

result_ma <- NULL; result_rho <- NULL; result_scr <- NULL; result_vac <- NULL; result_df <- NULL
clusterExport(cl, c("result_ma", "result_rho", "result_scr", "result_vac", "result_df"))

n_ti <- n_pi <- c_test <- c_vac <- c_df <- rep(0, n_year)
ms <- wgt <- sup <- pur <- bob <- ai <- tx <- rep(0, n_year) 

clusterExport(cl, c("n_ti", "n_pi", "c_test", "c_vac", "c_df", "ms", "wgt", "sup", "pur", "bob", "ai", "tx"))

# Setting the BVDV simulation fuction
bvdsim <- function(seed, n_ma, betas, rho, SCR, VAC, DF) {
  
  .C("outbreak_dairy", as.integer(seed), as.integer(n_ma), as.double(betas), as.double(rho), as.double(p_neighbour),
     as.integer(SCR), as.integer(VAC), as.integer(DF),
     as.integer(n_ti), as.integer(n_pi), as.double(ms), as.double(wgt), as.double(sup), as.double(pur), 
     as.double(bob), as.double(ai), as.double(tx), as.double(c_test), as.double(c_vac), as.double(c_df))
}

clusterExport(cl, c('bvdsim'))

# Set the different seeds for the core
clusterApply(cl, seq_along(cl), function(i) seed <<- (i-1)*lims)

# Test whether BVDV is endemic and estimate the economic loss
t_init <- proc.time()
result <- 
  clusterEvalQ(cl, {
    
    # Loading the compiled C code
    dyn.load("D:/Temp/Within/outbreak_dairy.dll")
    
    for (i in 3:3) {
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
                  ms <- result[[q+2]]
                  wgt <- result[[q+3]]
                  sup <- result[[q+4]]
                  pur <- result[[q+5]]
                  bob <- result[[q+6]]
                  ai <- result[[q+7]]
                  tx <- result[[q+8]]

                  test <- result[[q+9]]
                  vac <- result[[q+10]]
                  df <- result[[q+11]]

                  tti <- 0; tpi <- 0;
                  tms <- 0; twgt <- 0; tsup <- 0; tpur <- 0; tbob <- 0; 
                  tai <- 0; ttx <- 0; tgfr <- 0; ttest <- 0; tvac <- 0; tdf <- 0; tcost <- 0;
                  period <- 4
                  for (q in (n_year-period):n_year) {
                    tti <- tti + n_ti[q];
                    tpi <- tpi + n_pi[q];
                    
                    tms <- tms + ms[q]/((1.019)^(q-(n_year-period)))
                    twgt <- twgt + wgt[q]/((1.019)^(q-(n_year-period)))
                    tsup <- tsup + sup[q]/((1.019)^(q-(n_year-period)))
                    tpur <- tpur + pur[q]/((1.019)^(q-(n_year-period)))
                    tbob <- tbob + bob[q]/((1.019)^(q-(n_year-period)))
                    tai <- tai + ai[q]/((1.019)^(q-(n_year-period)))
                    ttx <- ttx + tx[q]/((1.019)^(q-(n_year-period)))
                    tgfr <- tgfr + (ms[q] + wgt[q] + sup[q] + bob[q] - pur[q] - ai[q] - tx[q])/((1.019)^(q-(n_year-period)))
                    ttest <- ttest + test[q]/((1.019)^(q-(n_year-period)))
                    tvac <- tvac + vac[q]/((1.019)^(q-(n_year-period)))
                    tdf <- tdf + df[q]/((1.019)^(q-(n_year-period)))
                    tcost <- tcost + (test[q] + vac[q] + df[q])/((1.019)^(q-(n_year-period)))
                  }
                  
                  tti <- tti
                  tpi <- tpi
                  
                  iter <- iter + 1
                  seed <- seed + 1
                  
                  result_ma <- c(result_ma, list_ma[i])
                  result_rho <- c(result_rho, list_rho[j])
                  result_scr <- c(result_scr, list_scr[l])
                  result_vac <- c(result_vac, list_vac[m])
                  result_df <- c(result_df, list_df[n])
                  
                  result_tti <- c(result_tti, tti)
                  result_tpi <- c(result_tpi, tpi)
                  
                  result_tms <- c(result_tms, tms); result_twgt <- c(result_twgt, twgt); result_tsup <- c(result_tsup, tsup)
                  result_tpur <- c(result_tpur, tpur); result_tbob <- c(result_tbob, tbob) 
                  result_tai <- c(result_tai, tai); result_ttx <- c(result_ttx, ttx); result_tgfr <- c(result_tgfr, tgfr); 
                  result_tctest <- c(result_tctest, ttest); result_tcvac <- c(result_tcvac, tvac); 
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
      result_tti, result_tpi, result_tms, result_twgt, result_tsup, result_tpur, result_tbob, result_tai, result_ttx, result_tgfr, 
      result_tctest, result_tcvac, result_tcdf, result_tcost)
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
  result_tms <- c(result_tms, result[[i]][[8]])
  result_twgt <- c(result_twgt, result[[i]][[9]])
  result_tsup <- c(result_tsup, result[[i]][[10]])
  result_tpur <- c(result_tpur, result[[i]][[11]])
  result_tbob <- c(result_tbob, result[[i]][[12]])
  result_tai <- c(result_tai, result[[i]][[13]])
  result_ttx <- c(result_ttx, result[[i]][[14]])
  result_tgfr <- c(result_tgfr, result[[i]][[15]])
  result_tctest <- c(result_tctest, result[[i]][[16]])
  result_tcvac <- c(result_tcvac, result[[i]][[17]])
  result_tcdf <- c(result_tcdf, result[[i]][[18]])
  result_tcost <- c(result_tcost, result[[i]][[19]])
}

dat <- data.frame("n_ma"= result_ma, "rho"= result_rho, "scr"= result_scr, "vac"= result_vac, "df"= result_df, 
                  "ti"= result_tti, "pi"= result_tpi, "ms"= result_tms, "wgt"= result_twgt, "sup"= result_tsup, "pur"= result_tpur,
                  "bob"= result_tbob, "ai"= result_tai, "tx"= result_ttx, "gfr"= result_tgfr, "ctest"= result_tctest, 
                  "cvac"= result_tcvac, "cdf"= result_tcdf, "cost"= result_tcost)
dat$bal <- dat$gfr - dat$cost
