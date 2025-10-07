library(tidyverse)
library(lme4)
library(lmerTest)

exp_pxg <- read_csv("mme_exp_main_linkage_inv_gts.csv") %>% filter(!(pop=="RGen0"))

sub_exp_pxg <- exp_pxg %>% dplyr::select(ind_id,pop,STD_Length,inv11_gt, inv18.1_gt,inv24.1_gt) %>% 
  pivot_longer(cols=inv11_gt:inv24.1_gt, names_to="inversion",values_to="ktype") %>% 
  mutate(inversion=str_remove(inversion,"_gt")) %>% 
  mutate(inversion = str_remove(inversion, "\\.1")) %>%
  mutate(trait="std. length (mm)") %>%
  mutate(pop2 = case_when(
    grepl("^R", pop) ~ "random",
    grepl("^U", pop) ~ "up", 
    grepl("^D", pop) ~ "down")) %>% 
  mutate(ktype.int = case_when(ktype=="SS"~"0",
                               ktype=="NS"~"1",
                               ktype=="NN"~"2")) %>% 
  mutate(ktype.int=as.numeric(ktype.int)) %>% 
  mutate(ktype.f = factor(ktype, levels=c("NS","SS","NN")))%>% 
  mutate(ktype = factor(ktype, levels=c("SS","NS","NN")))


### linear mixed-models
for(inv in c("inv11", "inv18", "inv24")) {
  cat("\n=== Linear mixed model for", inv, "===\n")
  dat <- subset(sub_exp_pxg, inversion == inv)
  
  lmm <- lmer(STD_Length ~ ktype.f + (1|pop), data = dat, REML = FALSE)
  print(summary(lmm))
}

### linear regression with selection treatment as an ordered predictor, 
## testing whether N allele count increases linearly from down→random→up

for(inv in c("inv11", "inv18", "inv24")) {
  cat("\n=== Linear trend for", inv, "===\n")
  dat <- subset(sub_exp_pxg, inversion == inv)
  dat$regime_ordered <- as.numeric(factor(dat$pop2, levels = c("down", "random", "up")))
  
  m <- lm(ktype.int ~ regime_ordered, data = dat)
  print(summary(m))
}


### plot experiment data

ggplot()+
  geom_violin(data=sub_exp_pxg,aes(x=ktype,y=STD_Length),fill="red",alpha=0.4,linewidth=0.41)+
  geom_jitter(data=sub_exp_pxg,aes(x=ktype,y=STD_Length),shape=1,size=1.5,stroke=0.51,width = 0.1,height=0)+
  stat_summary(data=sub_exp_pxg,
               aes(x=ktype,y=STD_Length),fun = "mean",
               geom = "crossbar", 
               width = 0.7)+
  #geom_text(data=annotation_data_exp, aes(x=1.5, y=-Inf, label=stats),vjust=-1, size=3) +
  #geom_text(data=annotation_data_exp, aes(x=1.5, y=-Inf, label=PVE), vjust=-3, size=3) +
  facet_wrap(~inversion,scales="free",nrow=1)+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position="none",
        panel.spacing.x = unit(0.6, "cm"),
        strip.text=element_text(color="black",size=10),
        axis.text = element_text(color="black")) +
  xlab("Genotype")+ylab("standard length (mm)")


