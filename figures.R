# === load packages
pkgs <- c("car", "som.nn", "tidyverse", "ggplot2", "viridisLite")
sapply(pkgs, require, character.only = T)

# === helper functions
source("helper_fns.R")

# === read in the model and data
modf <- readRDS("data/fit_mf3_segment_c.rds")
pars <- readRDS("data/parameters.rds")
dat <- readRDS("data/list_mf3.rds")

# === load field data from the common garden and select the 2011-2012 censuses
dfcomp <- read.csv("data/MajorFlats_comp_df1.csv")
dfcomp2 <- dfcomp[dfcomp$census == 2, ]
dfcomp3 <- dfcomp[dfcomp$census == 3, ]

# === observed group differences
Anova(lm(dat$growth_std~dat$subsppcyt))
fig2 <- dfcomp3 %>% 
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

# ggsave("Figures/figS2.pdf",
#        width = 160, height = 120, units = "mm", dpi = 300)
