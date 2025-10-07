library(tidyverse)
library(qtl2)

setwd("~/Desktop/Science Revision/QTL_code/")

#load data with 221 individuals (unknown sex assumed female - DROP BEFORE ANALYSIS)
mme_qtl <- read_cross2("mme_qtl.yaml")
summary(mme_qtl)

#remove individuals with no sex data
mme_qtl <- mme_qtl[c("-JP3I269","-JP3I272","-JP3I275","-JP3I276","-JP3I278","-JP3J281")]


#set covariate
sex.co <- setNames( (mme_qtl$covar$sex == "m")*1, rownames(mme_qtl$covar) )

#insert positions between markers (“pseudomarkers”) into genetic map 
map <- insert_pseudomarkers(mme_qtl$gmap, step=1)

#calculate conditional genotype probabilities
pr <- calc_genoprob(mme_qtl, map, error_prob=0.002, cores=6)

# qtl scan for all traits
out <- scan1(pr, mme_qtl$pheno, addcovar=sex.co)

outpeak <- find_peaks(out, map, threshold=4.2,drop=1.5)

out.covar <- scan1(pr, mme_qtl$pheno, addcovar=sex.co, intcovar=sex.co)

outdiff <- abs(out-out.covar) 

diffpeaks <- find_peaks(outdiff, map, threshold=2.5,drop=0.5) 

comout <- bind_rows("sex.add"=outpeak,"sex.co"=diffpeaks, .id="groups")%>% arrange(lodindex,chr)

newmap <- interp_map(map,mme_qtl$gmap,mme_qtl$pmap)

outpeakMB <- find_peaks(out, newmap, threshold=4.2,drop=1.5)

diffpeaksMB <- find_peaks(outdiff, newmap, threshold=2.5,drop=0.5) 

comoutMB <- bind_rows("sex.add"=outpeakMB,"sex.co"=diffpeaksMB, .id="groups") %>% arrange(lodindex,chr)

main_qtl_table <- left_join(comout,comoutMB,by=c("groups", "lodindex", "lodcolumn", "chr", "lod"), suffix = c(".cm",".mb"))  %>% 
  arrange(groups,chr) %>% 
  mutate(trait = recode(
      lodcolumn,
      vert_no = "Vertebral Count",
      grow_rate = "Growth Rate",
      logMO2 = "O2 Consumption",
      swim_speed = "Swim Speed",
      mineral = "Mineral Content",
      sqrt_lipid = "Lipid Content",
      morph.PC1 = "Shape PC1",
      res.morph.PC2 = "Shape PC2",
      res.lean = "Lean Content")) %>% dplyr::select(-lodcolumn,-lodindex) %>%
  dplyr::select(groups,trait,chr,pos.cm,ci_lo.cm,ci_hi.cm,pos.mb,ci_lo.mb,ci_hi.mb )


Chr11_87cM<- maxmarg(pr, map, chr=11, pos=86.92, return_char=TRUE) %>% 
  bind_cols(mme_qtl$covar$sex,mme_qtl$pheno[,"grow_rate"]) %>% dplyr::rename(Genotype=...1,Sex=...2,'Growth Rate'=...3) %>% 
  drop_na() %>% mutate(Genotype=as.factor(Genotype)) %>% mutate(inversion="inv11")

Chr18_40cM<- maxmarg(pr, map, chr=18, pos=40, return_char=TRUE) %>% 
  bind_cols(mme_qtl$covar$sex,mme_qtl$pheno[,"grow_rate"]) %>% rename(Genotype=...1,Sex=...2,'Growth Rate'=...3) %>% 
  drop_na() %>% mutate(Genotype=as.factor(Genotype)) %>% mutate(inversion="inv18")

Chr24_46cM<- maxmarg(pr, map, chr=24, pos=46, return_char=TRUE) %>% 
  bind_cols(mme_qtl$covar$sex,mme_qtl$pheno[,"grow_rate"]) %>% rename(Genotype=...1,Sex=...2,'Growth Rate'=...3) %>% 
  drop_na() %>% mutate(Genotype=as.factor(Genotype)) %>% mutate(inversion="inv24")
 
qtl_growth_effects <- bind_rows(Chr11_87cM,Chr18_40cM,Chr24_46cM) %>% 
  mutate(Genotype = recode(Genotype,SN = "NS")) %>% 
  mutate(Genotype = factor(Genotype, levels=c("SS","NS","NN")))



ggplot()+
  geom_violin(data=qtl_growth_effects,aes(x=Genotype,y=`Growth Rate`),fill="dodgerblue",alpha=0.4,linewidth =0.41)+
  geom_jitter(data=qtl_growth_effects,aes(x=Genotype,y=`Growth Rate`),shape=1,size=1.5,stroke = 0.51,width = 0.1,height=0)+
  stat_summary(data=qtl_growth_effects,
               aes(x=Genotype,y=`Growth Rate`),fun = "mean",
               geom = "crossbar", 
               width = 0.7)+
  facet_wrap(~inversion,scales="free",nrow=1)+  
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.position="none",
        panel.spacing.x = unit(0.6, "cm"),
        strip.text=element_text(color="black",size=10),
        axis.text = element_text(color="black")) +
  xlab("Genotype")+ylab("millimeter per day")


