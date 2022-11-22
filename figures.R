# === load packages
pkgs <- c("car", "som.nn", "tidyverse", "ggplot2", "viridisLite")
sapply(pkgs, require, character.only = T)

# === helper functions
source("helper_fns.R")



# === read in the model and data
modf <- readRDS("data/fit_mf3_segment_c.rds")
pars <- readRDS("data/parameters.rds")
dat <- readRDS("data/list_mf3.rds")

# === read simulation results and decode the entries
res <- readRDS("predictions_res.rds")
# - decode
cf <- res[[1]]; mu_dd <- res[[2]]; mu_base <- res[[3]]; pred_dd <- res[[2]]; pred_base <- res[[5]];
# ---

fres <- readRDS("predictions_fstats.rds")
# - decode
fstatf <- fres$fstatf; fstatb <- fres$fstatb;  
corcoef <- fres$corcoef; fstatdiff <- fres$fstatdiff; 
fstatdiff1 <- fres$fstatdiff1
# ---
# === load field data from the common garden and select the 2011-2012 censuses
dfcomp <- read.csv("data/MajorFlats_comp_df1.csv")
dfcomp2 <- dfcomp[dfcomp$census == 2, ]
dfcomp3 <- dfcomp[dfcomp$census == 3, ]

# === extract size measurements
dead <- which(dfcomp3$surv == 0)
growth <- dfcomp3$growth_std[-dead]
size1 <- dfcomp2$size_t[-dead]
size2 <- dfcomp3$size_t[-dead]
subspp <- dfcomp3$subsppcyt[-dead]

# 
spacing <- seq(.5, 4, by = .25) 
N <- dat$N
nn <- nrow(pars[[1]])
type <- as.numeric(as.factor(subspp))


# === Figure 1: F-statistic::difference 
pal <- viridis(nn)
mu.mean <- apply(fstatdiff1, 2, mean)
mu.HPDI <- apply(fstatdiff1, 2, rethinking::PI, prob=.95)
mu.HPDI.sd <- apply(fstatdiff1, 2, rethinking::PI, prob=.68)

df <- data.frame(spacing = spacing, mu = mu.mean, 
                 mu.l = mu.HPDI[1,], mu.u = mu.HPDI[2,],
                 musd.l = mu.HPDI.sd[1,], musd.u = mu.HPDI.sd[2,])
df %>% 
  ggplot() +
  geom_ribbon(aes(x = spacing, ymin = mu.l, ymax = mu.u), fill = pal[1], alpha = 0.25) +
  geom_ribbon(aes(x = spacing, ymin = musd.l, ymax = musd.u), fill = pal[1], alpha = 0.35) +
  geom_line(aes(x = spacing, y = mu), size = 2, linetype = "dashed") + 
  geom_hline(yintercept = 1, size = .5, linetype = "dashed", colour = "gray") + 
  labs(x = "Spacing between plants, [m]", y = "Relative F-statistic") +
  theme_bw(base_size = 14) + 
  theme(panel.border = element_blank(), panel.grid = element_blank(), 
        axis.line = element_line())

# export figure
# ggsave("Figures/fig1.svg", width = 130, height = 100, units = "mm", dpi = 300)

# === Figure 2: correlation betwenn full and base models
mu.mean <- apply(corcoef, 2, mean)
mu.HPDI <- apply(corcoef, 2, rethinking::PI, prob=.95)
mu.HPDI.sd <- apply(corcoef, 2, rethinking::PI, prob=.68)

df <- data.frame(spacing = spacing, mu = mu.mean, 
                 mu.l = mu.HPDI[1,], mu.u = mu.HPDI[2,],
                 musd.l = mu.HPDI.sd[1,], musd.u = mu.HPDI.sd[2,])

df %>% 
  ggplot() +
  geom_ribbon(aes(x = spacing, ymin = mu.l, ymax = mu.u), fill = pal[1000], alpha = 0.25) +
  geom_ribbon(aes(x = spacing, ymin = musd.l, ymax = musd.u), fill = pal[1000], alpha = 0.35) +
  geom_line(aes(x = spacing, y = mu), size = 2, linetype = "dashed") + 
  geom_vline(xintercept = 2.25, size = .5, linetype = "dashed", colour = "gray") + 
  labs(x = "Spacing between plants, [m]", y = "Pearson correlation") +
  theme_bw(base_size = 14) + 
  theme(panel.border = element_blank(), panel.grid = element_blank(), 
        axis.line = element_line())

# ggsave("Figures/fig2.png", width = 130, height = 100, units = "mm", dpi = 300)

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

dd_summ125 <- mu_dd[4,,] %>% t() %>% 
  bind_cols(type = as.factor(type)) %>% 
  pivot_longer(cols = 1:2000) %>% 
  select(-name) %>%
  group_by(type) %>% 
  mutate(mean = mean(value), 
         upper = mean(value) + sd(value), 
         lower = mean(value) - sd(value)) %>%
  select(-value) %>%
  distinct(.keep_all = TRUE) %>% 
  mutate(trtm = 1)

dd_summ225 <- mu_dd[8,,] %>% t() %>% 
  bind_cols(type = as.factor(type)) %>% 
  pivot_longer(cols = 1:2000) %>% 
  select(-name) %>%
  group_by(type) %>% 
  mutate(mean = mean(value), 
         upper = mean(value) + sd(value), 
         lower = mean(value) - sd(value)) %>%
  select(-value) %>%
  distinct(.keep_all = TRUE) %>% 
  mutate(trtm = 2)

dd_summ400 <- mu_dd[15,,] %>% t() %>%
  bind_cols(type = as.factor(type)) %>% 
  pivot_longer(cols = 1:2000) %>% 
  select(-name) %>%
  group_by(type) %>% 
  mutate(mean = mean(value), 
         upper = mean(value) + sd(value), 
         lower = mean(value) - sd(value)) %>%
  select(-value) %>%
  distinct(.keep_all = TRUE) %>% 
  mutate(trtm = 3) %>% 
  bind_rows(dd_summ, dd_summ125, dd_summ225) %>% 
  mutate(type = fct_recode(type, "A.arbuscula" = "1",
                           "A.tridentata-2x" = "2",
                           "A.tridentata-4x" = "3", "A.vaseyana-2x" = "4", 
                           "A.vaseyana-4x" = "5", "A.wyomingensis-4x" = "6")) 

dd_summ400 %>% 
  ggplot(aes(x = trtm, y = mean, colour = type)) + 
  geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.1) + 
  geom_point(size = 3) +
  geom_line() + 
  scale_color_manual(values = viridis(ncol(pars[[1]]))) +
  scale_x_discrete(limits = c(0, 1, 2, 3), labels = c("0.5","1.25","2.25", "4"), expand=c(0,.1))+
  labs(y = expression(paste("Average growth, [", m^3, " ", month^{-1}, "]")), 
       x = "Distance between plants [m]",
       colour = "Plant type") + 
  theme_classic(base_size = 14) 

ggsave("Figures/fig3.png", width = 160, height = 90, units = "mm", dpi = 300)



# === [not used] difference between full and base models
fstatf %>% 
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

# ggsave("fig5.pdf",fig5,
#        width=140,height=80,units=c("mm"),dpi=300,path = "~/Downloads/")



# === Figure S2: observed group differences
Anova(lm(dat$growth_std~dat$subsppcyt))
dfcomp3 %>% 
  filter(!is.na(growth_std)) %>% 
  mutate(Type = as.factor(subsppcyt)) %>%
  mutate(Type = fct_recode(Type, "A.tridentata-2x"="T2n",
                           "A.tridentata-4x"="T4n","A.vaseyana-2x"="V2n",
                           "A.vaseyana-4x"="V4n","A.wyomingensis-4x"="W4n",
                           "A.arbuscula"="AR")) %>%
  ggplot(aes(y = growth_std, x = Type)) +
  geom_boxplot() + 
  labs(subtitle="", 
       y=expression(paste("Growth [",m^{3}~month^{-1},"]")), 
       x = "Plant type") + theme_bw() +
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle = 35, vjust = 1, hjust=1))

# ggsave("Figures/figS2.png", width = 160, height = 120, units = "mm", dpi = 300)

