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
#clusterSetRNGStream(cl, iseed= 1) # Setup seed for RNG

#stopCluster(cl)

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

# Individual year
yr_sim <- 5
dat <- data.frame("year"= c(1:yr_sim), "n_ma"= rep(result_ma, each= yr_sim), "rho"= rep(result_rho, each= yr_sim), 
                  "nb"= rep(result_nb, each= yr_sim), 
                  "ti"= result_ti, "pi"= result_pi, "ms"= result_ms, "wgt"= result_wgt, "sup"= result_sup, "pur"= result_pur, "bob"= result_bob, 
                  "ai"= result_ai, "tx"= result_tx, "gfr"= result_gfr,
                  "test"= result_ctest, "vac"= result_cvac, "df"= result_cdf, "cost"= result_cost)
dat$type <- paste(dat$rho, dat$nb, sep= "-")
dat <- dat[dat$type == "0-0" | dat$type == "0.04-0.02", ]
#dat <- dat[dat$type == "0-0" | dat$type == "0.04-0.02" | dat$type == "0-0.02", ]
#dat_1 <- dat[dat$n_ma == 100, ]
#dat_2 <- dat[dat$n_ma == 250, ]
dat_3 <- dat[dat$n_ma == 400, ]

gd <- ggplot(data= dat_3, aes(x= as.factor(year), y= gfr, fill= type)) + 
  geom_boxplot(outlier.shape= NA, width= 0.8) +
  theme_bw() +
  labs(x= "\nProduction year", y= "Gross farm revenue (NZ$)\n", colour= "black") +
  scale_fill_brewer(palette= "Set2") +
  scale_x_discrete(expand= c(0.12, 0)) +
  scale_y_continuous(expand= c(0,0), breaks= c(1100000, 1150000, 1200000, 1250000), limits= c(1080000, 1250000), 
                     labels= c("1.10M", "1.15M", "1.20M", "1.25M")) +
  theme(
    legend.position= "none",
    axis.line= element_line(colour= "black"),
    #axis.line.y= element_blank(),
    #axis.line.x= element_blank(),
    panel.grid.minor= element_blank(),
    panel.border= element_blank(),
    panel.background= element_blank(),
    axis.title.x= element_text(size= 11),
    axis.title.y= element_text(size= 11),
    axis.text.x= element_text(size= 10, colour= "black"),
    axis.text.y= element_text(size= 10, colour= "black"),
    aspect.ratio= 3/1)
gd

gd <- ggplotGrob(gd); gb <- ggplotGrob(gb)
g <- cbind(gd, gb, size= "first")
grid.arrange(g, ncol= 1)
g <- arrangeGrob(g, ncol= 1)
ggsave("D:/Temp/Within/_Impact.tiff", g, device= "tiff", dpi= 500)

yr_sim <- 5
dat <- data.frame("year"= c(1:yr_sim), "n_ma"= rep(result_ma, each= yr_sim), "rho"= rep(result_rho, each= yr_sim), 
                  "nb"= rep(result_nb, each= yr_sim), "n"= result_ti, "typ"= "ti")
dat <- rbind(dat, data.frame("year"= c(1:yr_sim), "n_ma"= rep(result_ma, each= yr_sim), "rho"= rep(result_rho, each= yr_sim), 
                             "nb"= rep(result_nb, each= yr_sim), "n"= result_pi, "typ"= "pi"))
dat$type <- paste(dat$rho, dat$nb, sep= "-")
dat <- dat[dat$type == "0.04-0.02", ]
#dat_1 <- dat[dat$n_ma == 100, ]
#dat_2 <- dat[dat$n_ma == 250, ]
dat_3 <- dat[dat$n_ma == 400, ]
lab <- c("Transiently infected", "Persistently infected")
names(lab) <- c("ti", "pi")
gd2 <- ggplot(data= dat_3, aes(x= as.factor(year), y= n)) + 
  geom_boxplot(outlier.shape= NA, width= 0.8) +
  facet_wrap( ~ typ, scales= "free_y", labeller= labeller(typ= lab)) +
  theme_bw() +
  labs(x= "\nProduction year", y= "Number of cattle (Dairy)\n", colour= "black") +
  theme(
    aspect.ratio= 3/1)
gd2
  
gd2 <- ggplotGrob(gd2); gb2 <- ggplotGrob(gb2)
g <- cbind(gd2, gb2, size= "first")
grid.arrange(g, ncol= 1)
g <- arrangeGrob(g, ncol= 1)
ggsave("D:/Temp/Within/_Dynamics.tiff", g, device= "tiff", dpi= 500)

# As a total sum
dat <- data.frame("n_ma"= result_ma, "rho"= result_rho, "nb"= result_nb, "ti"= result_tti, "pi"= result_tpi,
                  "ms"= result_tms, "wgt"= result_twgt, "sup"= result_tsup, "pur"= result_tpur, "bob"= result_tbob, 
                  "ai"= result_tai, "tx"= result_ttx, "gfr"= result_tgfr, 
                  "test"= result_tctest, "vac"= result_tcvac, "df"= result_tcdf, "cost"= result_tcost)
dat$type <- paste(dat$rho, dat$nb, sep= "-")
dat <- dat[dat$type == "0-0" | dat$type == "0.04-0.02", ]
#dat_1 <- dat[dat$n_ma == 100, ]
#dat_2 <- dat[dat$n_ma == 250, ]
dat_3 <- dat[dat$n_ma == 400, ]
options(scipen= 999)
print(summary(lm(gfr ~ as.factor(type), data= dat_3)), digits= 10)
42197.46/(5*400); (42197.46 - 1551.01*1.96)/(5*400); (42197.46 + 1551.01*1.96)/(5*400)

#save.image("D:/Temp/Within/bvd_impact_dairy.RData")
rm(list= ls())
load("D:/Temp/Within/bvd_impact_dairy.RData")
load("D:/Temp/Within/bvd_impact_beef.RData")

#############################
dat <- data.frame("n_ma"= result_ma, "rho"= result_rho, "nb"= result_nb, "ti"= result_tti, "pi"= result_tpi,
                  "ms"= result_tms, "wgt"= result_twgt, "sup"= result_tsup, "pur"= result_tpur, "bob"= result_tbob, 
                  "ai"= result_tai, "tx"= result_ttx, "gfr"= result_tgfr, 
                  "test"= result_tctest, "vac"= result_tcvac, "df"= result_tcdf, "cost"= result_tcost)
dat$type <- paste(dat$rho, dat$nb, sep= "-")
dat <- dat[dat$type == "0-0" | dat$type == "0.04-0.02", ]
summary(lm(dat$ms ~ dat$type))
summary(lm(dat$wgt ~ dat$type))
summary(lm(dat$sup ~ dat$type))
summary(lm(dat$pur ~ dat$type))
summary(lm(dat$bob ~ dat$type))
summary(lm(dat$ai ~ dat$type))
summary(lm(dat$tx ~ dat$type))

#############################

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
dat$control <- paste(dat$scr, "-", dat$vac, "-", dat$df, sep= "")

#dat_1 <- dat[dat$n_ma == 100, ]
#dat_2 <- dat[dat$n_ma == 250, ]
dat_3 <- dat[dat$n_ma == 400, ]
options(scipen= 999)
#summary(lm(cost ~ as.factor(control), data= dat_1))
#summary(lm(cost ~ as.factor(control), data= dat_2))
print(summary(lm(cost ~ as.factor(control), data= dat_3)), digits= 10)

#summary(lm(gfr ~ as.factor(control), data= dat_1))
#summary(lm(gfr ~ as.factor(control), data= dat_2))
print(summary(lm(gfr ~ as.factor(control), data= dat_3)), digits= 10)
#vcov(lm(gfr ~ as.factor(control), data= dat_2))

# Benefit-cost ratio
#as.numeric(lm(gfr ~ as.factor(control), data= dat_1)$coefficient[2:30] / lm(cost ~ as.factor(control), data= dat_1)$coefficient[2:30])
#as.numeric(lm(gfr ~ as.factor(control), data= dat_2)$coefficient[2:30] / lm(cost ~ as.factor(control), data= dat_2)$coefficient[2:30])
as.numeric(lm(gfr ~ as.factor(control), data= dat_3)$coefficient[2:30] / lm(cost ~ as.factor(control), data= dat_3)$coefficient[2:30])

# Large
gfr <- glm(gfr ~ as.factor(control), data= dat_3)
benefit <- as.double(gfr$coefficient[2:30])
pval <- coef(summary(gfr))[2:30, 4]

benefit <- ifelse(benefit < 0, 0, benefit)
benefit <- benefit / (5*400)
benefit[pval > 0.01] <- 0

cost <- as.numeric(glm(cost ~ as.factor(control), data= dat_3)$coefficient[2:30])
cost <- cost / (5*400)
lab <- NULL; one <- c("", "CT:", "SCR:"); two <- c("", "V1:", "V1M:", "V2:", "V2M:"); three <- c("", "DF")
for (i in 1:3) {
  for (j in 1:5) {
    for (k in 1:2) {
      lab <- c(lab, paste(one[i], two[j], three[k], sep= ""))
    }
  }
}
lab
lab[3] <- "V1"; lab[5] <- "V1M"; lab[7] <- "V2"; lab[9] <- "V2M"; lab[11] <- "CT"; lab[13] <- "CT:V1"; lab[15] <- "CT:V1M";
lab[17] <- "CT:V2"; lab[19] <- "CT:V2M"; lab[21] <- "SCR"; lab[23] <- "SCR:V1"; lab[25] <- "SCR:V1M"; lab[27] <- "SCR:V2"; 
lab[29] <- "SCR:V2M"
lab <- lab[2:30]

bpc <- benefit/cost
#bpc <- bpc[benefit != 0]
#cost <- cost[benefit != 0]
#lab <- lab[benefit != 0]
#benefit <- benefit[benefit != 0]
benefit <- benefit + 20
#benefit <- benefit[bpc > 0.5]
#lab <- lab[bpc > 0.5]
#cost <- cost[bpc > 0.5]
#bpc <- bpc[bpc > 0.5]
ord <- order(benefit)

g_dat <- data.frame("type"= 1, "value"= -cost[rev(ord)], "lab"= lab[rev(ord)], "bpc"= "")
g_dat <- rbind(g_dat, data.frame("type"= 2, "value"= benefit[rev(ord)], 
                                 "lab"= "", "bpc"= as.character(round(bpc[rev(ord)], 2))))
g_dat$control <- rep(c(length(bpc):1), 2)

gd <- ggplot(data= g_dat, aes(x= control, y= value)) +
  geom_bar(stat= "identity", aes(fill= as.factor(type)), width= 0.85, alpha= 0.75) +
  #geom_text(aes(label= lab),  size= 2.5, fontface= "bold", position= position_stack(1.0), hjust= 1.1, colour= "black") +
  geom_text(aes(label= bpc),  size= 2.5, fontface= "bold", hjust= -0.25, 
            colour= c(rep("black", length(bpc)), 
                      rep("red", 5), "black", rep("red", 4), "black", rep("red", 3), rep("black", 2), rep("red", 3), "black", rep("red", 2), "black", rep("white", 6))) +
            #colour= c(rep("black", length(bpc)), rep("black", length(bpc)))) +
  scale_fill_manual("legend", values= c("1"= "firebrick1", "2"= "limegreen")) +
  scale_x_continuous(expand= c(0, 0.25)) +
  scale_y_continuous(expand= c(0, 0), breaks= c(-20, -10, 0, 20, 20+10, 20+20), limits= c(-20-2, 2+20+20), 
                     labels= c("$20", "$10", "Cost", "Benefit", "$10", "$20")) +
  theme(legend.position= "none",
        axis.line= element_line(colour= "black"),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        panel.border= element_blank(),
        axis.line.y= element_blank(),
        panel.background= element_blank(),
        axis.title.x= element_blank(),
        axis.title.y= element_blank(),
        axis.text.x = element_text(size= 7.5, colour= "black"),
        axis.text.y = element_blank(),
        aspect.ratio= 3/2) +
  coord_flip()
gd
#ggsave("D:/Temp/Within/_Cb_dairy_ref.tiff", gd, device= "tiff", dpi= 100)
#save.image("D:/Temp/Within/bvd_control_dairy.RData")
rm(list= ls())
load("D:/Temp/Within/bvd_control_dairy.RData")
load("D:/Temp/Within/bvd_control_beef.RData")

gd <- ggplotGrob(gd); gb <- ggplotGrob(gb)
g <- cbind(gd, gb, size= "first")
grid.arrange(g, ncol= 1)
g <- arrangeGrob(g, ncol= 1)
ggsave("D:/Temp/Within/_CostBenefit.tiff", g, device= "tiff", dpi= 500)


###############################################################

# Cost effectiveness
table(dat_1$pi)
pois <- hurdle(pi ~ as.factor(control) + offset(log(rep(list_ma[1], nrow(dat_1)))), data= dat_1, dist= "poisson", zero.dist= "binomial")
nb <- hurdle(pi ~ as.factor(control) + offset(log(rep(list_ma[1], nrow(dat_1)))), data= dat_1, dist= "negbin", zero.dist= "binomial")
lrtest(pois, nb)
summary(nb)
(1-exp(nb$coefficient[[1]][2:30]))[coef(summary(nb))[[1]][2:30, 4] < 0.05]
(1 - exp(nb$coefficient[[2]][1] + nb$coefficient[[2]][2:30])/(1+exp(nb$coefficient[[2]][1] + nb$coefficient[[2]][2:30])) / 
    exp(nb$coefficient[[2]][1])/(1+exp(nb$coefficient[[2]][1])))[coef(summary(nb))[[2]][2:30, 4] < 0.05]

table(dat_2$pi)
pois <- hurdle(pi ~ as.factor(control) + offset(log(rep(list_ma[2], nrow(dat_2)))), data= dat_2, dist= "poisson", zero.dist= "binomial")
nb <- hurdle(pi ~ as.factor(control) + offset(log(rep(list_ma[2], nrow(dat_2)))), data= dat_2, dist= "negbin", zero.dist= "binomial")
lrtest(pois, nb)
summary(nb)
(1-exp(nb$coefficient[[1]][2:30]))[coef(summary(nb))[[1]][2:30, 4] < 0.05]
(1 - exp(nb$coefficient[[2]][1] + nb$coefficient[[2]][2:30])/(1+exp(nb$coefficient[[2]][1] + nb$coefficient[[2]][2:30])) / 
    exp(nb$coefficient[[2]][1])/(1+exp(nb$coefficient[[2]][1])))[coef(summary(nb))[[2]][2:30, 4] < 0.05]

table(dat_3$pi)
pois <- hurdle(pi ~ as.factor(control) + offset(log(rep(list_ma[3], nrow(dat_3)))), data= dat_2, dist= "poisson", zero.dist= "binomial")
nb <- hurdle(pi ~ as.factor(control) + offset(log(rep(list_ma[3], nrow(dat_3)))), data= dat_2, dist= "negbin", zero.dist= "binomial")
lrtest(pois, nb)
summary(nb)
(1-exp(nb$coefficient[[1]][2:30]))[coef(summary(nb))[[1]][2:30, 4] < 0.05]
(1 - exp(nb$coefficient[[2]][1] + nb$coefficient[[2]][2:30])/(1+exp(nb$coefficient[[2]][1] + nb$coefficient[[2]][2:30])) / 
  exp(nb$coefficient[[2]][1])/(1+exp(nb$coefficient[[2]][1])))[coef(summary(nb))[[2]][2:30, 4] < 0.05]

summary(lm(gfr ~ ti + pi, data= dat_1))
summary(lm(gfr ~ ti + pi, data= dat_2))
summary(lm(gfr ~ ti + pi, data= dat_3))

#####
gfr <- glm(gfr ~ as.factor(control), data= dat_1)
benefit <- as.double(gfr$coefficient)
pval <- coef(summary(gfr))[, 4]
benefit <- ifelse(benefit < 0, 0, benefit)
benefit <- benefit / (5*100)
benefit[pval > 0.01] <- 0
cost <- as.numeric(glm(cost ~ as.factor(control), data= dat_1)$coefficient)
cost <- cost / (5*100)
tmp <- data.frame("gfr"= benefit, "cost"= cost, "bcr"= benefit / cost)
one <- c("No", "CT", "SCR"); two <- c("No", "V1", "V1M", "V2", "V2M"); three <- c("No", "DF"); l <- 1
for (i in 1:3) {
  for (j in 1:5) {
    for (k in 1:2) {
      tmp$Test[l] <- paste(one[i], "-", three[k], sep= "")
      tmp$Vaccine[l] <- two[j]
      tmp$bcrcol[l] <- ifelse(tmp$bcr[l] > 0, "black", "white")
      l <- l + 1
    }
  }
}
tmp <- tmp[2:30, ]

ggplot(data= tmp, aes(x= cost, y= gfr, size= bcr)) +
  geom_point(aes(colour= Vaccine, fill= Vaccine, shape= as.factor(Test)), alpha= 0.5, stroke= 1.5) +
  scale_size(range = c(0.1, 15)) +
  scale_color_brewer(palette= "Set1", name= "Vaccination", breaks= c("No", "V1", "V2", "V1M", "V2M")) +
  scale_shape_manual(values= c(16, 1, 17, 2, 15, 0), name= "Annual test", breaks= c("No-DF", "CT-DF", "SCR-DF"),
                     labels= c("No test", "Calf test", "Screening")) +
  scale_x_continuous(expand= c(0, 0), limits= c(0, 30), breaks= c(0, 10, 20, 30)) +
  scale_y_continuous(expand= c(0, 0), limits= c(0, 30), breaks= c(0, 10, 20, 30)) +
  xlab("\nControl cost (NZ$/cow)") + ylab("Benefit (NZ$/cow)\n") +
  geom_text_repel(aes(label= round(bcr, digits= 2)), size= 3, 
                  colour= tmp$bcrcol) +
  guides(size= F, fill= F) +
  theme(#legend.position= "bottom",
    axis.line= element_line(colour= "black"),
    panel.grid.major= element_line(colour= "grey"),
    panel.grid.minor= element_blank(),
    panel.border= element_blank(),
    panel.background= element_blank(),
    axis.title.x= element_text(size= 12.5, colour= "black"),
    axis.title.y= element_text(size= 12.5, colour= "black"),
    axis.text.x= element_text(size= 10, colour= "black"),
    axis.text.y= element_text(size= 10, colour= "black"),
    aspect.ratio= 1/1)
#####