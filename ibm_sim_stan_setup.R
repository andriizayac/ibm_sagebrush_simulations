library(rstan)
library(som.nn)
library(dplyr)
library(ggplot2)
library(reshape2)
library(rethinking)
library(viridisLite)
library(ggridges)
library(forcats)
virpalette <- viridis(5)
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE) 

dist.fx <- function(x1, x2, y1, y2) {sqrt((x2-x1)^2 + (y2-y1)^2)}
toroid.dist <- function(x1,y1,x2,y2,xmax,ymax){
  xcutoff=xmax/2
  ycutoff=ymax/2
  dx=abs(x2-x1)
  dy=abs(y2-y1)
  if (dx > xcutoff){
    dx=xmax-dx
  } 
  if(dy > ycutoff){
    dy=ymax-dy
  }
  return(sqrt(dx^2+dy^2))
}
cf.fn <- function(model=NULL,data=NULL){
  n = data$N
  nssp = data$k
  post = extract(model)
  iter = dim(post[[1]])[1]
  # extract parameters
  a01 = with(post,a01)
  a2 = with(post,a2)
  # spatial data
  sobs = data$size_observations
  dobs = data$dist_observations^2
  n_nb = data$n_nb
  pos = data$pos
  type = data$type
  cf_mod = matrix(NA,nrow=iter,ncol=n)
  for(i in 1:iter){
    for(j in 1:n){
      smat = sobs[pos[j] : (pos[j] + n_nb[j] - 1)]
      dmat = dobs[pos[j] : (pos[j] + n_nb[j] - 1)]
      cf_mod[i,j]=sum(smat^a01[i] / exp(dmat * a2[i]))
    }
  }
  return(cf_mod)
}

# ======================================== Simulations with mf3 model 

# === export figure ####
# library(export)
# path="~/Desktop/Orchard/Writing_IBM/ibm_figures/"
# graph2svg(file=paste0(path,"fig2.svg"),width=4.5,height=3)

# === read in the models
modf = readRDS("fit_mf3_segment_c.rds")
#modb = readRDS("~/Desktop/Orchard/OrchardProjectR/CaughlinLab_sagebrushABM/mf3base.rds")
dat = readRDS("list_mf3.rds")
dfcomp = read.csv("MajorFlats_comp_df1.csv")
dfcomp2 = dfcomp[dfcomp$census == 2, ]
dfcomp3 = dfcomp[dfcomp$census == 3, ]
pars=extract(modf)

# === re-run the models on the untransformed growth scale 
# --- can wait for now

# === extract size measurements
dead <- which(dfcomp3$surv == 0)
growth <- dfcomp3$growth_std[-dead]
size1 <- dfcomp2$size_t[-dead]
size2 <- dfcomp3$size_t[-dead]
set.seed(123)
growth_std <- growth # sample(growth)
set.seed(123)
sizet <- size1 # sample(size1)
set.seed(123)
sizetmat <- size2 # sample(size2)
set.seed(123)
subspp <- dfcomp3$subsppcyt[-dead] # sample(dfcomp3$subsppcyt[-dead])

N <- length(size1)

sizemat <- matrix(NA, N, N)
for(i in 1:N){ sizemat[, i] = sizetmat[i] }

# === simulate dist matrices with a toroid wrapper
#outs = which(dfcomp3$surv == 0) # create index to take out dead plants

nr <- 18 # number of rows in the common garden
nc <- 26 # number of cols in the common garden

# = create a list of distmatrixes
ndistmat_list <- list()
tdistmat_list <- list()
sizemat_list <- list()
spacing <- seq(.5, 4, by = .5) # spacing
for(l in 1:length(spacing)){
  #dist_increment=seq(1,l=n,by=l)
  #start_point=dist_increment[2]-dist_increment[1]
  x <- rep(seq(0, l = nr, by = spacing[l]), times = nc) 
  y <- sort(rep(seq(0, l = nc, by = spacing[l]), times = nr))

  set.seed(123)
  outs <- sample(1:468,468-N)

  ndistmat=matrix(NA, 468, 468)[-outs,-outs]
  tdistmat=matrix(NA, 468, 468)[-outs,-outs]
  for(i in 1:N){
    for(j in 1:N){
      ndistmat[i,j] = dist.fx(x[i], x[j], y[i], y[j])
      tdistmat[i,j] = toroid.dist(x[i], y[i], x[j], y[j], max(x)+spacing[l], max(y)+spacing[l])
    }
  }
  ndistmat <- ifelse(ndistmat > 4 | ndistmat == 0,  NA, ndistmat)
  ndistmat_list[[l]] <- ndistmat
  tdistmat <- ifelse(tdistmat > 4 | tdistmat == 0, NA, tdistmat)
  tdistmat_list[[l]] <- tdistmat
  sizemat1 <- ifelse(is.na(ndistmat), NA, sizemat)
  sizemat_list[[l]] <- sizemat1
}

# === run stan model
post <- extract(modf)
set.seed(123)
ind <- sample(1:2000, 150) # 1:2000 # 
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
          a1 = with(post,a1)[ind,],
          a03 = with(post, a03)[ind], 
          a01 = with(post, a01)[ind],
          sigma = with(post, sigma)[ind], 
          sigma_a3 = with(post, sigma_a3)[ind], 
          # sigma_a1 = with(post, sigma_a1)[ind],
          nn = length(post$sigma_Itype[ind])) 

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
             #sigma_a1 = with(post, sigma_a1)[ind],
             nn = length(post$sigma_Itype[ind])) 
}


# === STAN loop
# simmod <- stan_model("ibm_simulation_stan.stan")
# gc()
# predf_list = list()
# predb_list = list()
# for(ll in 1:length(spacing)){
#   print(ll)
#   gc()
#   m <- sampling(simmod, data = aa_list[[ll]], iter=2, warmup=1, chains=1,
#                 save_warmup = FALSE, seed = 123, algorithm = "Fixed_param",
#            pars = c("pred_dd","pred_base"), include = TRUE)
# 
#   mpost <- extract(m)
#   rm(m)
#   predf_list[[ll]] <- mpost$pred_dd[1, , ]
#   predb_list[[ll]] <- mpost$pred_base[1, , ]
# }

# === R loop
cf_dat <- cf.fn(modf, dat)
cf_datmu = apply(cf_dat, 2, mean)

N = aa_list[[1]]$N
nn = aa_list[[1]]$nn
cf <- array(dim = c(length(spacing), nn, N))
mu_dd <- array(dim = c(length(spacing), nn, N))
mu_base <- array(dim = c(length(spacing), nn, N))
pred_dd <- array(dim = c(length(spacing), nn, N))
pred_base <- array(dim = c(length(spacing), nn, N))

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
  #cf_std <- matrix(NA, nrow = nn, ncol = N)

  # --- main loop
  for(n in 1:nn){
    for(i in 1:N){
      cf[ll, n, i] = sum( (smat[i, ]*a01[n]) / exp(dmat[i, ]^2*a2[n]), na.rm = TRUE); #
    }
    #cf_std[n, ] = (cf[n, ] - mean(cf[n, ])) / (2 * sd(cf[n, ]))
    for(l in 1:N){
      
      mu_dd[ll, n, l] = Itype[n, type[l]] + #c[n,type[l]]*sizet_std[l] + a3[n,type[l]]*cf_std[n, l]
         # back transform patameters to predict on the original scale
         -(c[n,type[l]] / (2*sd(sizet))) * mean(sizet) +  
         (c[n,type[l]] / (2*sd(sizet))) * sizet[l] +  
         #-(a3[n,type[l]] / (2*sd(cf_dat[l, ]))) * mean(cf_dat[l, ]) +  #
         (a3[n,type[l]] / (2*sd(cf_datmu))) * cf[ll, n, l];  #
      mu_base[ll, n, l] = Itype[n, type[l]] + #c[n,type[l]]*sizet_std[l]
        # back transform patameters to predict on the original scale
        -(c[n,type[l]] / (2*sd(sizet))) * mean(sizet) + #  
        (c[n,type[l]]  / (2*sd(sizet))) * sizet[l];
        #-(a3[n,type[l]] / (2*sd(cf_dat[l, ]))) * mean(cf_dat[l, ]);
    }
    pred_dd[ll, n, ] = rnorm(N, mu_dd[ll, n, ], sigma[n])
    pred_base[ll, n, ] = rnorm(N, mu_base[ll, n, ], sigma[n])
  }
}
# with increased spacing the crowding gets progressively smaller
crowding <- apply(cf, c(1,2), mean)
matplot(as.matrix(crowding), type = "l") # confirms crowding decreases with spacing
# ---
# cf_std <- array(dim = c(length(spacing), nn, N))
# for(ll in 1:length(spacing)){
#   for(n in 1:nn){
#     cf_std[ll, n, ] = (cf[ll, n, ] - mean(cf_datmu)) / (2*sd(cf_datmu))
#   }
# }
# crowding_std <- apply(cf_std, c(1,2), mean)
# matplot(as.matrix(crowding_std), type = "l") # confirms crowding decreases with spacing

ix <- 3 # spacing option
boxplot(dat$growth_std~dat$subsppcyt)
boxplot(apply(pred_base[ix, , ], 2, mean) ~ type, pch = 19, col = rgb(.5,0,.5,.5)) # full model
plot(apply(pred_dd[ix, , ], 2, mean) ~ jitter(type), pch = 16, col = rgb(.5,1,.5,.5)) # base model


# explore correlations
matplot(corcoef, type = "l")

# F-st on the difference
matplot(fstatdiff, type = "l")

# === figures and summaries
print("Range of F-statistic (full model): "); range(fstatf)
print("Range of F-statistic (base model): "); range(fstatb)
print("Range of Spearman correlations (base and full models): "); range(corcoef)
print("Difference b/w average base and full F-statistics: "); plot(apply(fstatdiff,1,median), type = "l")

# full model F-statistic
rownames(fstatf) <- spacing
fstatf %>% t() %>% as.data.frame() %>% melt() %>% 
  ggplot(aes(y = value, x = variable)) +
  geom_boxplot(fill = virpalette[1]) +
  geom_hline(yintercept = median(fstatb), colour = virpalette[3], size = 1) +
  geom_hline(yintercept = quantile(fstatb, probs = c(.025, .975)), colour = virpalette[3], size = 1, linetype = "dashed") +
  labs(x = "Distance (m)", y = "F-statistic", title = "Full model") +
  theme_bw() + 
  theme(axis.title = element_text(size = 16), title = element_text(size = 18))

#  full model F-statistic
par(mgp = c(2,.75,0), mar = c(3,3,0,0)+0.5)
mu.mean <- apply(fstatf, 1, median)
mu.HPDI <- apply(fstatf, 1, HPDI, prob=.95)
plot(mu.mean ~ spacing, type = "l", lty = "dashed", bty = "n", ylim = c(0,300), 
     lwd = 5.5, col = virpalette[1],
     xlab="Distance [m]", ylab="F-statistic", main = "")
abline(h = median(fstatb), lwd = 2, col = virpalette[3], lty = "solid")
abline(h = quantile(fstatb, probs = c(.025, .975)), lwd = 1.5, col = virpalette[3], lty = "dotted")
shade(mu.HPDI, spacing, col = col.alpha(virpalette[1], .15))

# difference betwenn full and base models in F-statistic
mu.mean <- apply(fstatdiff, 1, median)
mu.HPDI <- apply(fstatdiff, 1, PI, prob=.95)
plot(mu.mean ~ spacing, type = "l", lty = "dashed", bty = "n",
     lwd = 5.5, col = virpalette[2],
     xlab="Distance [m]", ylab="Difference [F-statistic]", main = "Full vs Base model")
shade(mu.HPDI, spacing, col = col.alpha(virpalette[2], .15))


# correlation betwenn full and base models
pal <- viridis(2000)
par(mgp = c(2,.75,0), mar = c(3,3,0,0)+0.5)

mu.mean <- apply(corcoef, 1, median)
mu.HPDI <- apply(corcoef, 1, PI, prob=.95)
plot(1,type = "n", ylim = c(-.5, 1), xlim = range(spacing), bty = "n",
     xlab = "Distance [m]", ylab = "Correlation")
# for(i in 1:nn) { 
#   lines(spacing, corcoef[, i], col = col.alpha(pal[i], .05))
# }
# lines(spacing, mu.mean, lty = "dashed", lwd = 6, col = "black", 
#      xlab = "Distance [m]", ylab = "Correlation")
lines(spacing, mu.mean, lty = "dashed", col = virpalette[3], lwd = 5)
abline(v = c(1, 1.5), lty = "dotted", col = virpalette[1], lwd = 2)
shade(mu.HPDI,spacing,col=col.alpha(virpalette[3],.15))

# Fstatistic full model
pal <- viridis(2000)
mu.mean <- apply(fstatf, 1, median)
mu.HPDI <- apply(fstatf, 1, PI, prob=.95)
plot(1,type = "n", ylim = c(0, 325), xlim = range(spacing), bty = "n",
     xlab="Distance [m]", ylab = "F-statistic")
for(i in 1:nn) { 
  lines(spacing, fstatf[, i], col = col.alpha(pal[i], .05))
}
lines(spacing, mu.mean, lty = "dashed", lwd = 6, col = "black", 
      xlab = "Distance [m]", ylab = "Correlation")
abline(v = c(1, 1.5), lty = "dotted", col = virpalette[1], lwd = 2)
real <- summary(lm(dfcomp3$growth_std ~ dfcomp3$subsppcyt))$fstatistic[1]
abline(h = real, lty = "solid", col = virpalette[5], lwd = 2)
# shade(mu.HPDI,spacing,col=col.alpha(virpalette[3],.15))


# === correlation by species between mean predictions (no stochasticity)
cormu <-matrix(NA, ncol=nn, nrow = length(spacing))
corsubspp <- array(dim = c(length(spacing), nn, 5))
t2x <- which(subspp == "T2n")
t4x <- which(subspp == "T4n")
v2x <- which(subspp == "V2n")
v4x <- which(subspp == "V4n")
w4x <- which(subspp == "W4n")
subspplist <- list(t2x, t4x, v2x, v4x, w4x)
for(i in 1:length(spacing)){
  for(j in 1:nn){
    for(h in 1:5){
      corsubspp[i,j, h] = cor(mu_dd[i, j, subspplist[[h]]], mu_base[i, j, subspplist[[h]]])
    }
    cormu[i,j] = cor(mu_dd[i, j, ], mu_base[i, j, ])
  }
}

mu.mean <- apply(cormu, 1, mean)
mu.HPDI <- apply(cormu, 1, HPDI, prob=.95)
mu.mean1 <- apply(corsubspp[, ,1], 1, mean)
mu.mean2 <- apply(corsubspp[, ,2], 1, mean)
mu.mean3 <- apply(corsubspp[, ,3], 1, mean)
mu.mean4 <- apply(corsubspp[, ,4], 1, mean)
mu.mean5 <- apply(corsubspp[, ,5], 1, mean)

plot(1, type="n", xlim = range(spacing), ylim = c(0, 1), xlab="", ylab="")
shade(mu.HPDI,spacing,col=col.alpha('black',.15))
lines(spacing, mu.mean1, type = "l", lty = "dashed", lwd = 5, col = virpalette[1], xlab="", ylab="")
lines(spacing, mu.mean2, type = "l", lty = "dashed", lwd = 5, col = virpalette[2], xlab="", ylab="")
lines(spacing, mu.mean3, type = "l", lty = "dashed", lwd = 5, col = virpalette[3], xlab="", ylab="")
lines(spacing, mu.mean4, type = "l", lty = "dashed", lwd = 5, col = virpalette[4], xlab="", ylab="")
lines(spacing, mu.mean5, type = "l", lty = "dashed", lwd = 5, col = virpalette[5], xlab="", ylab="")

# --- raw common garden
pal <- viridis(6)
dfcomp3 %>% 
  #filter(!Subspecies=="AR") %>% droplevels() %>%
  mutate(subsppcyt = fct_recode(subsppcyt, "A.arbuscula" = "AR",
                                "A.tridentata:2x" = "T2n",
                                "A.tridentata:4x" = "T4n", "A.vaseyana:2x" = "V2n", 
                                "A.vaseyana:4x" = "V4n", "A.wyomingensis:4x" = "W4n")) %>%
  ggplot(aes(x = x, y = y, size = size_t, colour = subsppcyt)) + 
  geom_point() +
  scale_size(range=c(0,2))+
  labs(x = "x coordinate [m]", y = "y coordinate [m]") + #, title = "Common garden map") +
  scale_colour_manual(name = 'Suspecies:cytotype', values = pal) + #, labels = c('T2n','T4n','V2n','V4n','W4n'))
  guides(size = FALSE, colour = guide_legend("Subspecies:cytotype")) + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.minor=element_line(size=.2,colour="gray",linetype="dashed"),
        panel.grid.major=element_line(size=.2,colour="gray",linetype="dashed"),
        #panel.background =element_rect(fill=col.alpha("black",0.05)),
        text=element_text(size=16),
        #axis.ticks = element_blank(),
        axis.line  = element_line(color="black", size = .1),
        #axis.text = element_blank(),
        axis.title=element_text(size=12))


