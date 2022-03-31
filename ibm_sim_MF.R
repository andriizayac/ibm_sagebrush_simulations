# This script simulates the outcomes of MajorsFlat CG growth with manipulated plant interspaces

pkgs <- c("som.nn", "dplyr", "tidyverse", "rstan", "ggplot2", "viridisLite")
sapply(pkgs, require, character.only = T)

# === helper functions
# euclidean distance function
dist.fx <- function(x1, x2, y1, y2) {sqrt((x2-x1)^2 + (y2-y1)^2)}
# eucledean distance on the toroid wrap
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
# competition kernel function based reproducing the relationship from the statistical model
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


# === read in the models and data
modf <- readRDS("fit_mf3_segment_c.rds")
dat <- readRDS("list_mf3.rds")
# actual data from the common garden
dfcomp <- read.csv("MajorFlats_comp_df1.csv")
dfcomp2 <- dfcomp[dfcomp$census == 2, ]
dfcomp3 <- dfcomp[dfcomp$census == 3, ]
pars <- rstan::extract(modf)

# === observed group differences
anova(lm(dat$growth_std~dat$subsppcyt))
fig2 <- dfcomp3 %>% 
  mutate(Type = as.factor(subsppcyt)) %>%
  mutate(Type = fct_recode(Type, "A.tridentata:2x"="T2n",
                                      "A.tridentata:4x"="T4n","A.vaseyana:2x"="V2n",
                                      "A.vaseyana:4x"="V4n","A.wyomingensis:4x"="W4n",
                                    "A.arbuscula"="AR")) %>%
  ggplot(aes(y = growth_std, x = Type)) +
  geom_boxplot() + 
  labs(subtitle="", 
       y=expression(paste("Growth [",m^{3}~month^{-1},"]")), 
       x = "Plant type") + theme_bw() +
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# ggsave("fig2.pdf",fig2,
#        width=160,height=120,units=c("mm"),dpi=300,path = "~/Downloads/")

# === extract size measurements
dead <- which(dfcomp3$surv == 0)
growth <- dfcomp3$growth_std[-dead]
size1 <- dfcomp2$size_t[-dead]
size2 <- dfcomp3$size_t[-dead]
subspp <- dfcomp3$subsppcyt[-dead]
N <- length(size1)

set.seed(123)
# perturbation results in qualitatively similar results, and are not presented in the paper
pert <- sample(1:N, N) # perturb the spatial arrangement
growth_std <- growth#[pert]
sizet <- size1#[pert]  
sizemat <- matrix( rep(sizet, N), N)
subspp <- subspp#[pert]

type <- as.numeric(as.factor(subspp))

# === simulate dist matrices with a toroid wrapper
#outs = which(dfcomp3$surv == 0) # create index to take out dead plants

nr <- 18 # number of rows in the common garden
nc <- 26 # number of cols in the common garden

# = create a list of distmatrixes
ndistmat_list <- list()
tdistmat_list <- list()
sizemat_list <- list()
spacing <- seq(.5, 4, by = .25) # spacing
for(l in 1:length(spacing)){
  #dist_increment=seq(1,l=n,by=l)
  #start_point=dist_increment[2]-dist_increment[1]
  x <- rep(seq(0, l = nr, by = spacing[l]), times = nc) 
  y <- sort(rep(seq(0, l = nc, by = spacing[l]), times = nr))
  
  set.seed(123)
  outs <- sample(1:468, 468-N)
  
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
              a1 = with(post,a1)[ind,],
              a03 = with(post, a03)[ind], 
              a01 = with(post, a01)[ind],
              sigma = with(post, sigma)[ind], 
              sigma_a3 = with(post, sigma_a3)[ind], 
              # sigma_a1 = with(post, sigma_a1)[ind],
              nn = length(post$sigma_Itype[ind])) 

# data inputs for the simulation
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

# === R simulation loop
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
# save results
#saveRDS(list(cf, mu_dd, mu_base, pred_dd, pred_base), "predictions_res.rds")
res <- readRDS("predictions_res.rds")
cf <- res[[1]]; mu_dd <- res[[2]]; mu_base <- res[[3]]; pred_dd <- res[[2]]; pred_base <- res[[5]];
# --- with increased spacing the crowding gets progressively smaller
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

ix <- 1 # spacing option
boxplot(dat$growth_std~dat$subsppcyt, ylim = c(-.1, .1))
points(apply(pred_base[ix, , ], 2, mean) ~ type, pch = 19, col = rgb(.5,0,.5,.5)) # full model
points(apply(pred_dd[ix, , ], 2, mean) ~ jitter(type), pch = 16, col = rgb(.5,1,.5,.5)) # base model


# Calculate F-statisic

fstatf <- matrix(NA, nn, length(spacing)) # full model
fstatb <- matrix(NA, nn, length(spacing)) # base model
corcoef <- matrix(NA, nn, length(spacing)) # correlation btw full and base models
fstatdiff <- matrix(NA, nn, length(spacing)) # difference btw full and base
fstatdiff1 <- matrix(NA, nn, length(spacing))

resdiff <- array(dim = c(nn, length(spacing), N))
for(i in 1:length(spacing)) {
  for(j in 1:nn) {
  fstatf[j,i] <- summary(lm(pred_dd[i,j,] ~ type))$fstatistic[[1]]
  fstatb[j,i] <- summary(lm(pred_base[i,j,] ~ type))$fstatistic[[1]]
  corcoef[j,i] <- cor(mu_base[i,j,], mu_dd[i,j,])
  fstatdiff[j,i] <- summary(lm( (pred_dd[i,j,]-pred_base[i,j,]) ~ type))$fstatistic[[1]]
  resdiff[j,i,] <- pred_base[i,j,] - pred_dd[i,j,]
  } 
  fstatdiff1[,i] <- fstatdiff[,i]/median(fstatb[,i])
}

res <- apply(resdiff, c(2,3), mean)
boxplot(t(res)[,15] ~ dat$subsppcyt)

boxplot(apply(pred_base, c(1,3), mean)[1,] ~ dat$subsppcyt)

# === figures and summaries
cat("Range of F-statistic (full model): ", range(fstatf)) 
cat("Range of F-statistic (base model): ", range(fstatb)) 
cat("Range of Spearman correlations (base and full models): ", range(corcoef)) 
cat("Difference b/w average base and full F-statistics: ", range(fstatdiff))
plot(apply(fstatdiff,2, median), type = "l")

# === Figure 1: F-statisti::difference 
pal <- viridis(nn)
mu.mean <- apply(fstatdiff1, 2, meadian)
mu.HPDI <- apply(fstatdiff1, 2, rethinking::PI, prob=.95)
mu.HPDI.sd <- apply(fstatdiff1, 2, rethinking::PI, prob=.68)

# export figure
# pdf("~/Desktop/Orchard/Writing_IBM/ibm_figures/fig1.pdf", width=6, height=4.5)
plot(1,type = "n", ylim = c(0, 1.5), #ylim = c(0, 150),#
     xlim = range(spacing), bty = "n",
     xlab="Distance [m]", 
     ylab = "Relative F-statistic")
# for(i in 1:nn) { 
#   lines(spacing, fstatdiff[i, ], col = rethinking::col.alpha(pal[i], .05))
# }
lines(spacing, mu.mean, lty = "dashed", lwd = 6, col = "black", 
      xlab = "Distance [m]", ylab = "Correlation")
# lines(spacing, mu.HPDI[1,], lty = "dashed", lwd = 1, 
#       xlab = "Distance [m]", ylab = "Correlation", 
#       col = pal[1])
# lines(spacing, mu.HPDI[2,], lty = "dashed", lwd = 1,
#       xlab = "Distance [m]", ylab = "Correlation", 
#       col = pal[1])
rethinking::shade(mu.HPDI, spacing, lty = "dotted", 
                  col = rethinking::col.alpha(pal[1], 0.25))
rethinking::shade(mu.HPDI.sd, spacing, lty = "dotted", 
                  col = rethinking::col.alpha(pal[1], 0.35))
# dev.off()


# === Figure 2: correlation betwenn full and base models
pal <- viridis(2000)
par(mgp = c(2,.75,0), mar = c(3,3,0,0)+0.5)

#pdf("~/Desktop/Orchard/Writing_IBM/ibm_figures/fig2.pdf", width=6, height=4.5)
mu.mean <- apply(corcoef, 2, median)
mu.HPDI <- apply(corcoef, 2, rethinking::PI, prob=.95)
mu.HPDI.sd <- apply(corcoef, 2, rethinking::PI, prob=.68)
plot(1,type = "n", ylim = c(-1, 1), xlim = range(spacing), bty = "n",
     xlab = "Distance [m]", ylab = "Correlation")
# for(i in 1:nn) { 
#   lines(spacing, corcoef[, i], col = col.alpha(pal[i], .05))
# }
# lines(spacing, mu.mean, lty = "dashed", lwd = 6, col = "black", 
#        xlab = "Distance [m]", ylab = "Correlation")
# abline(v = c(1, 1.5), lty = "dotted", col = pal[1], lwd = 2)
rethinking::shade(mu.HPDI, spacing, col=rethinking::col.alpha(pal[1000],.25))
rethinking::shade(mu.HPDI.sd, spacing, col=rethinking::col.alpha(pal[1000],.35))
lines(spacing, mu.mean, lty = "dashed", col = "black", lwd = 6)

dev.off()

# === Figure 3: 
dd_summ <- mu_dd[1,,] %>% t() %>% 
  bind_cols(type = as.factor(type)) %>% 
  pivot_longer(cols = 1:2000) %>% 
  select(-name) %>%
  group_by(type) %>% 
  mutate(mean = mean(value), 
         upper = mean(value) + sd(value), 
         lower = mean(value) - sd(value)) %>%
  select(-value) %>%
  distinct(.keep_all = TRUE) %>% 
  mutate(trtm = 0) 

mu_dd[15,,] %>% t() %>%
  bind_cols(type = as.factor(type)) %>% 
  pivot_longer(cols = 1:2000) %>% 
  select(-name) %>%
  group_by(type) %>% 
  mutate(mean = mean(value), 
         upper = mean(value) + sd(value), 
         lower = mean(value) - sd(value)) %>%
  select(-value) %>%
  distinct(.keep_all = TRUE) %>% 
  mutate(trtm = 1) %>% 
  bind_rows(dd_summ) %>% 
  mutate(type = fct_recode(type, "A.arbuscula" = "1",
                                "A.tridentata:2x" = "2",
                                "A.tridentata:4x" = "3", "A.vaseyana:2x" = "4", 
                                "A.vaseyana:4x" = "5", "A.wyomingensis:4x" = "6")) %>%
  ggplot(aes(x = trtm, y = mean, colour = type)) + 
  geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.1) + 
  geom_point(size = 3) +
  geom_line() + 
  scale_color_manual(values = viridis(alist$k)) +
  scale_x_discrete(limits = c(0, 1), labels = c("0.5", "4"), expand=c(0,.1))+
  labs(y = expression(paste("Average growth, [", m^3/month^{-1}, "]")), 
       x = "Distance between plants [m]",
       colour = "Subspecies:cytotype") + 
  theme_classic(base_size = 14) 
  #theme(legend.position = "none")

ggsave("fig4.pdf",
       width=160,height=90,units=c("mm"),dpi=300,
       path = "~/Desktop/Orchard/Writing_IBM/ibm_figures/")

# === Figure 5: difference between full and base models
fig5 <- fstatf %>% 
  as.data.frame() %>% 
  pivot_longer(cols = 1:length(spacing)) %>% 
  ggplot(aes(x = name, y = value, color = name, fill = name)) +
  see::geom_violinhalf() +
  labs(x = "Spacing [m]", y = "F-statistic") +
  scale_x_discrete(labels = as.character(spacing)) +
  scale_colour_manual(values = rep(pal[4], length(spacing))) +
  scale_fill_manual(values = rep(pal[4], length(spacing))) +
  geom_hline(yintercept = mean(fstatb), colour = pal[2], size = 1, linetype = "solid") +
  geom_hline(yintercept = quantile(fstatb, probs = c(0.025, 0.975)), colour = pal[2], size = .75, linetype = "dashed") +  
  annotate("rect", xmin = -Inf, xmax = Inf, 
            ymin = quantile(fstatb, 0.025), 
            ymax = quantile(fstatb, 0.975), fill = pal[2], alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=14), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave("fig5.pdf",fig5,
              width=140,height=80,units=c("mm"),dpi=300,path = "~/Downloads/")
       
  

# === baseline ANOVA results
pd <- apply(pred_base, c(2,3), mean)
fbase <- apply(pd, 1, function(x) { summary(lm(x ~ dat.df$Type))$fstatistic[[1]] })
# === densest ANOVA results
f05 <- apply(pred_dd[1,,], 1, function(x) { summary(lm(x ~ dat.df$Type))$fstatistic[[1]] })
mean(f05)
quantile(f05, probs = c(0.025, 0.975))

# === raw common garden plot
pal <- viridisLite::viridis(6)
fig1 <- dfcomp3 %>% 
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
# ggsave("fig1.pdf",fig1,
#       width=140,height=80,units=c("mm"),dpi=300,path = "~/Downloads/")
