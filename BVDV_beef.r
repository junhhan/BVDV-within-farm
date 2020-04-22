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
#clusterSetRNGStream(cl, iseed= 1) # Setup seed for RNG

#stopCluster(cl)

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

# Individual year
yr_sim <- 5
dat <- data.frame("year"= c(1:yr_sim), "n_ma"= rep(result_ma, each= yr_sim), "rho"= rep(result_rho, each= yr_sim), 
                  "nb"= rep(result_nb, each= yr_sim), 
                  "ti"= result_ti, "pi"= result_pi, "wgt"= result_wgt, "sup"= result_sup, "pur"= result_pur, "gfr"= result_gfr, 
                  "test"= result_ctest, "vac"= result_cvac, "df"= result_cdf, "cost"= result_cost)
dat$type <- paste(dat$rho, dat$nb, sep= "-")
dat <- dat[dat$type == "0-0" | dat$type == "0.04-0.02", ]
#dat <- dat[dat$type == "0-0" | dat$type == "0.04-0.02" | dat$type == "0-0.02", ]
#dat_1 <- dat[dat$n_ma == 50, ]
dat_2 <- dat[dat$n_ma == 150, ]
#dat_3 <- dat[dat$n_ma == 300, ]

gb <- ggplot(data= dat_2, aes(x= as.factor(year), y= gfr, fill= type)) + 
  geom_boxplot(outlier.shape= NA, width= 0.8) +
  theme_bw() +
  labs(x= "\nProduction year", y= "\n\n", colour= "black") +
  scale_fill_brewer(palette= "Set2") +
  scale_x_discrete(expand= c(0.12, 0)) +
  scale_y_continuous(expand= c(0,0), breaks= c(175000, 225000, 275000, 325000), limits= c(160000, 325000), 
                     labels= c("175K", "225K", "275K", "325K")) +
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
gb

yr_sim <- 5
dat <- data.frame("year"= c(1:yr_sim), "n_ma"= rep(result_ma, each= yr_sim), "rho"= rep(result_rho, each= yr_sim), 
                  "nb"= rep(result_nb, each= yr_sim), "n"= result_ti, "typ"= "ti")
dat <- rbind(dat, data.frame("year"= c(1:yr_sim), "n_ma"= rep(result_ma, each= yr_sim), "rho"= rep(result_rho, each= yr_sim), 
                             "nb"= rep(result_nb, each= yr_sim), "n"= result_pi, "typ"= "pi"))
dat$type <- paste(dat$rho, dat$nb, sep= "-")
dat <- dat[dat$type == "0.04-0.02", ]
#dat_1 <- dat[dat$n_ma == 50, ]
dat_2 <- dat[dat$n_ma == 150, ]
#dat_3 <- dat[dat$n_ma == 300, ]
lab <- c("Transiently infected", "Persistently infected")
names(lab) <- c("ti", "pi")
gb2 <- ggplot(data= dat_2, aes(x= as.factor(year), y= n)) + 
  geom_boxplot(outlier.shape= NA, width= 0.8) +
  facet_wrap( ~ typ, scales= "free_y", labeller= labeller(typ= lab)) +
  theme_bw() +
  labs(x= "\nProduction year", y= "Number of cattle (Beef)\n", colour= "black") +
  theme(
    aspect.ratio= 3/1)
gb2

# As a total sum
dat <- data.frame("n_ma"= result_ma, "rho"= result_rho, "nb"= result_nb, "ti"= result_tti, "pi"= result_tpi,
                  "wgt"= result_twgt, "sup"= result_tsup, "pur"= result_tpur, "gfr"= result_tgfr, 
                  "test"= result_tctest, "vac"= result_tcvac, "df"= result_tcdf, "cost"= result_tcost)
dat$type <- paste(dat$rho, dat$nb, sep= "-")
dat <- dat[dat$type == "0-0" | dat$type == "0.04-0.02", ]
#dat_1 <- dat[dat$n_ma == 50, ]
dat_2 <- dat[dat$n_ma == 150, ]
#dat_3 <- dat[dat$n_ma == 300, ]
options(scipen= 999)
print(summary(lm(gfr ~ as.factor(type), data= dat_2)), digits= 10)
30422.98/(5*150); (30422.98 - 906.85*1.96)/(5*150); (30422.98 + 906.85*1.96)/(5*150)

#save.image("D:/Temp/Within/bvd_impact_beef.RData")
rm(list= ls())
load("D:/Temp/Within/bvd_impact_beef.RData")


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
dat$control <- paste(dat$scr, "-", dat$vac, "-", dat$df, sep= "")

#dat_1 <- dat[dat$n_ma == 50, ]
dat_2 <- dat[dat$n_ma == 150, ]
#dat_3 <- dat[dat$n_ma == 300, ]
options(scipen= 999)
#summary(lm(cost ~ as.factor(control), data= dat_1))
print(summary(lm(cost ~ as.factor(control), data= dat_2)), digits= 10)
#summary(lm(cost ~ as.factor(control), data= dat_3))

#summary(lm(gfr ~ as.factor(control), data= dat_1))
print(summary(lm(gfr ~ as.factor(control), data= dat_2)), digits= 10)
#summary(lm(gfr ~ as.factor(control), data= dat_3))

# Benefit-cost ratio
#as.numeric(lm(gfr ~ as.factor(control), data= dat_1)$coefficient[2:30] / lm(cost ~ as.factor(control), data= dat_1)$coefficient[2:30])
as.numeric(lm(gfr ~ as.factor(control), data= dat_2)$coefficient[2:30] / lm(cost ~ as.factor(control), data= dat_2)$coefficient[2:30])
#as.numeric(lm(gfr ~ as.factor(control), data= dat_3)$coefficient[2:30] / lm(cost ~ as.factor(control), data= dat_3)$coefficient[2:30])

# Moderate
gfr <- glm(gfr ~ as.factor(control), data= dat_2)
benefit <- as.double(gfr$coefficient[2:30])
pval <- coef(summary(gfr))[2:30, 4]

benefit <- ifelse(benefit < 0, 0, benefit)
benefit <- benefit / (5*150)
benefit[pval > 0.01] <- 0

cost <- as.numeric(glm(cost ~ as.factor(control), data= dat_2)$coefficient[2:30])
cost <- cost / (5*150)
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
benefit <- benefit + 50
#benefit <- benefit[bpc > 0.5]
#lab <- lab[bpc > 0.5]
#cost <- cost[bpc > 0.5]
#bpc <- bpc[bpc > 0.5]
ord <- order(benefit)

g_dat <- data.frame("type"= 1, "value"= -cost[rev(ord)], "lab"= lab[rev(ord)], "bpc"= "")
g_dat <- rbind(g_dat, data.frame("type"= 2, "value"= benefit[rev(ord)], 
                                 "lab"= "", "bpc"= as.character(round(bpc[rev(ord)], 2))))
g_dat$control <- rep(c(length(bpc):1), 2)

gb <- ggplot(data= g_dat, aes(x= control, y= value)) +
  geom_bar(stat= "identity", aes(fill= as.factor(type)), width= 0.85, alpha= 0.75) +
  #geom_text(aes(label= lab),  size= 2.5, fontface= "bold", position= position_stack(1.0), hjust= 1.1, colour= "black") +
  geom_text(aes(label= bpc),  size= 2.5, fontface= "bold", hjust= -0.25, 
            colour= c(rep("black", length(bpc)), 
                      rep("red", 7), "black", rep("red", 6), rep(c("black", "red"), 2), rep("black", 2), rep("red", 2), rep("black", 4), rep("white", 3))) +
            #colour= c(rep("black", length(bpc)), rep("black", length(bpc)))) +
  scale_fill_manual("legend", values= c("1"= "firebrick1", "2"= "limegreen")) +
  scale_x_continuous(expand= c(0, 0.25)) +
  scale_y_continuous(expand= c(0, 0), breaks= c(-50, -25, 0, 50, 50+25, 50+50), limits= c(-50-5, 5+50+50), 
                     labels= c("$50", "$25", "Cost", "Benefit", "$25", "$50")) +
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
gb
#ggsave("D:/Temp/Within/_Cb_beef_ref.tiff", gb, device= "tiff", dpi= 100)
save.image("D:/Temp/Within/bvd_control_beef.RData")
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

