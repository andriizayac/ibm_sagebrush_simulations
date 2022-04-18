# This script simulates the outcomes of MajorsFlat CG growth with manipulated plant interspaces

# === load packages
pkgs <- c("car", "som.nn", "tidyverse", "ggplot2", "viridisLite")
sapply(pkgs, require, character.only = T)

# === helper functions
source("helper_fns.R")

# === read in the model and data
# rstan object from Zaiats et al. 2021 (10.1002/ecy.3502)
modf <- readRDS("data/fit_mf3_segment_c.rds") 
# list of parameters from modf
pars <- readRDS("data/parameters.rds") 
# data from the common garden from Zaiats et al. 2021 (10.1002/ecy.3502)
dat <- readRDS("data/list_mf3.rds") 

# === load field data from the common garden and select the 2011-2012 censuses
dfcomp <- read.csv("data/MajorFlats_comp_df1.csv")
dfcomp2 <- dfcomp[dfcomp$census == 2, ]
dfcomp3 <- dfcomp[dfcomp$census == 3, ]

# === extract plant size measurements (exclude dead plants)
dead <- which(dfcomp3$surv == 0)

growth <- dfcomp3$growth_std[-dead]
size1 <- dfcomp2$size_t[-dead]
size2 <- dfcomp3$size_t[-dead]
subspp <- dfcomp3$subsppcyt[-dead]

# total number of live plants
N <- length(size1)

set.seed(123)
# perturbation results in qualitatively similar results, and are not presented in the paper
pert <- sample(1:N, N) # perturb the spatial arrangement
# comment in `[pert]` to reproduce perturbed arrangement
growth_std <- growth#[pert]
sizet <- size1#[pert]  
sizemat <- matrix( rep(sizet, N), N)
subspp <- subspp#[pert]

type <- as.numeric(as.factor(subspp))

# === simulate distance matrices (toroid wrapper not used but created)
nr <- 18 # number of rows in the common garden
nc <- 26 # number of cols in the common garden

# --- create a list of distance matrices for each simulation treatment
ndistmat_list <- list()
tdistmat_list <- list()
sizemat_list <- list()
# spacing array for plant interspaces from 0.5 to 4m
spacing <- seq(.5, 4, by = .25) 

for(l in 1:length(spacing)){
  x <- rep(seq(0, l = nr, by = spacing[l]), times = nc) 
  y <- sort(rep(seq(0, l = nc, by = spacing[l]), times = nr))
  
  set.seed(123)
  outs <- sample(1:468, 468-N)
  
  ndistmat=matrix(NA, 468, 468)[-outs,-outs]
  tdistmat=matrix(NA, 468, 468)[-outs,-outs]
  for(i in 1:N){
    for(j in 1:N){
      ndistmat[i,j] = dist.fx(x[i], x[j], y[i], y[j])
      tdistmat[i,j] = toroid.dist(x[i], y[i], x[j], y[j], max(x) + spacing[l], max(y) + spacing[l])
    }
  }
  ndistmat <- ifelse(ndistmat > 4 | ndistmat == 0,  NA, ndistmat)
  ndistmat_list[[l]] <- ndistmat
  tdistmat <- ifelse(tdistmat > 4 | tdistmat == 0, NA, tdistmat)
  tdistmat_list[[l]] <- tdistmat
  sizemat1 <- ifelse(is.na(ndistmat), NA, sizemat)
  sizemat_list[[l]] <- sizemat1
}

# === create an input list with demographic variables and simulations
# --- a template 
post <- pars
set.seed(123)
ind <- 1:2000 # sample(1:2000, 150) # for shorter simulations
alist <- list(N = N, k = 6, 
              sizet = sizet, growth_std = growth_std,
              type = as.numeric(as.factor(subspp)), sizemat = sizemat_list[[1]], 
              distmat = ndistmat_list[[1]], tdistmat = tdistmat_list[[1]],
              # parameter estimates
              Itype = with(post, Itype)[ind,], 
              a0 = with(post, a0)[ind], 
              sigma_Itype = with(post, sigma_Itype)[ind],
              c = with(post, c)[ind,], 
              c0 = with(post, c0)[ind], 
              sigma_c = with(post, sigma_c)[ind],
              a3 = with(post, a3)[ind,], 
              a2 = with(post, a2)[ind], 
              a1 = with(post, a1)[ind,],
              a03 = with(post, a03)[ind], 
              a01 = with(post, a01)[ind],
              sigma = with(post, sigma)[ind], 
              sigma_a3 = with(post, sigma_a3)[ind], 
              nn = length(post$sigma_Itype[ind])) 

# using the template, create data inputs for all simulations (each spacing treatment)
aa_list <- list()
for(m in 1:length(spacing)){
  aa_list[[m]] <- list(N = N, k = 6, 
                       sizet = sizet, growth_std = growth_std,
                       type = as.numeric(as.factor(subspp)), sizemat = sizemat_list[[m]], 
                       distmat = ndistmat_list[[m]], tdistmat = tdistmat_list[[m]],
                       # parameter estimates
                       Itype = with(post, Itype)[ind,], 
                       a0 = with(post, a0)[ind], 
                       sigma_Itype = with(post, sigma_Itype)[ind],
                       c = with(post, c)[ind,], 
                       c0 = with(post, c0)[ind], 
                       sigma_c = with(post, sigma_c)[ind],
                       a3 = with(post, a3)[ind,], 
                       a2 = with(post, a2)[ind], 
                       a1 = with(post,a1)[ind,],
                       a03 = with(post, a03)[ind], 
                       a01 = with(post, a01)[ind],
                       sigma = with(post, sigma)[ind], 
                       sigma_a3 = with(post, sigma_a3)[ind], 
                       nn = length(post$sigma_Itype[ind])) 
}


# --- create a crowding variable based on the data and parameters
# - used for re-scaling
cf_dat <- cf.fn(pars, dat)
cf_datmu = apply(cf_dat, 2, mean) # posterior averages
# ---

N = aa_list[[1]]$N # sample size
nn = aa_list[[1]]$nn # number of posterior samples

# empty arrays for storage
cf <- array(dim = c(length(spacing), nn, N)) 
mu_dd <- array(dim = c(length(spacing), nn, N))
mu_base <- array(dim = c(length(spacing), nn, N))
pred_dd <- array(dim = c(length(spacing), nn, N))
pred_base <- array(dim = c(length(spacing), nn, N))

# === R simulation loop

# - the following loop may take awhile. 
# - Consider loading the pre-saved simulation output: 
# res <- readRDS("predictions_res.rds")
# - decode the list after the loop
# cf <- res[[1]]; mu_dd <- res[[2]]; mu_base <- res[[3]]; pred_dd <- res[[2]]; pred_base <- res[[5]];

for(ll in 1:length(spacing)){
  print(ll)
  # --- input parameters
  N = aa_list[[ll]]$N
  nn = aa_list[[ll]]$nn
  type = aa_list[[ll]]$type
  sizet = aa_list[[ll]]$sizet
  sizet_std = (sizet - mean(sizet)) / (2*sd(sizet))
  Itype = aa_list[[ll]]$Itype
  c = aa_list[[ll]]$c
  a1 = aa_list[[ll]]$a1
  a01 = aa_list[[ll]]$a01
  a2 = aa_list[[ll]]$a2 
  a3 = aa_list[[ll]]$a3 
  sigma = aa_list[[ll]]$sigma
  smat = aa_list[[ll]]$sizemat
  dmat = aa_list[[ll]]$distmat

  # --- main loop
  for(n in 1:nn){
    for(i in 1:N){
      cf[ll, n, i] = sum( (smat[i, ]*a01[n]) / exp(dmat[i, ]^2*a2[n]), na.rm = TRUE); 
    }
    for(l in 1:N){
      
      mu_dd[ll, n, l] = Itype[n, type[l]] + 
        # back transform parameters to predict on the original scale
        -(c[n,type[l]] / (2*sd(sizet))) * mean(sizet) +  
        (c[n,type[l]] / (2*sd(sizet))) * sizet[l] +  
        (a3[n,type[l]] / (2*sd(cf_datmu))) * cf[ll, n, l];  
      
      mu_base[ll, n, l] = Itype[n, type[l]] + 
        # back transform parameters to predict on the original scale
        -(c[n,type[l]] / (2*sd(sizet))) * mean(sizet) + #  
        (c[n,type[l]]  / (2*sd(sizet))) * sizet[l];
    }
    pred_dd[ll, n, ] = rnorm(N, mu_dd[ll, n, ], sigma[n])
    pred_base[ll, n, ] = rnorm(N, mu_base[ll, n, ], sigma[n])
  }
}
# --- save results
#saveRDS(list(cf, mu_dd, mu_base, pred_dd, pred_base), "predictions_res.rds")
res <- readRDS("predictions_res.rds")
cf <- res[[1]]; mu_dd <- res[[2]]; mu_base <- res[[3]]; pred_dd <- res[[2]]; pred_base <- res[[5]];

# Calculate F-statisic's
fstatf <- matrix(NA, nn, length(spacing)) # full model
fstatb <- matrix(NA, nn, length(spacing)) # base model
corcoef <- matrix(NA, nn, length(spacing)) # correlation between full and base models
fstatdiff <- matrix(NA, nn, length(spacing)) # difference between full and base
fstatdiff1 <- matrix(NA, nn, length(spacing)) # difference between full and base models standardized by the median of the base model
resdiff <- array(dim = c(nn, length(spacing), N))

# --- The following loop may take awhile to run
# consider loading the pre-saved results:
outf <- readRDS("predictions_fstats.rds")
# - decode the results
fstatf <- outf$fstatf; fstatb <- outf$fstatb; corcoef <- outf$corcoef;
fstatdiff <- outf$fstatdiff; fstatdiff1 <- outf$fstatdiff1;

for(i in 1:length(spacing)) {
  for(j in 1:nn) {
  fstatf[j,i] <- Anova(lm(pred_dd[i,j,] ~ type))[[3]][1]
  fstatb[j,i] <- Anova(lm(pred_base[i,j,] ~ type))[[3]][1]
  corcoef[j,i] <- cor(mu_base[i,j,], mu_dd[i,j,])
  fstatdiff[j,i] <- Anova(lm( (pred_dd[i,j,]-pred_base[i,j,]) ~ type))[[3]][1]
  resdiff[j,i,] <- pred_base[i,j,] - pred_dd[i,j,]
  } 
  fstatdiff1[,i] <- fstatdiff[,i]/median(fstatb[,i])
}
# save the results of the f-statistic analysis
# saveRDS(list(fstatf = fstatf, fstatb = fstatb, 
#              corcoef = corcoef, fstatdiff = fstatdiff, 
#              fstatdiff1 = fstatdiff1), "predictions_fstats.rds")


res <- apply(resdiff, c(2,3), mean)

# === basic figures and summaries
cat("Range of F-statistic (full model): ", range(fstatf)) 
cat("Range of F-statistic (base model): ", range(fstatb)) 
cat("Range of Spearman correlations (base and full models): ", range(corcoef)) 
cat("Difference b/w average base and full F-statistics: ", range(fstatdiff))
plot(apply(fstatdiff,2, median), type = "l", xlab = "Spacing, [m]", ylab = "F-statistic difference")

  

# === baseline ANOVA results
pd <- apply(pred_base, c(2,3), mean)
fbase <- apply(pd, 1, function(x) { Anova(lm(x ~ type))[[3]][1] })


