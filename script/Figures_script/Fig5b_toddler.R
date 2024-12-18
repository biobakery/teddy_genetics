
#################
# load the data #
#################
source("~/script/packages_load.R")
library(tidyverse)
library(reshape2)
setwd("~/ddong/TEDDY/")
opts <- list()
opts$m <- "./output/microbiome_toddler.txt"
opts$g <- "./genetics.bed/filter/teddy_out8.qc.PCA_PC1-PC5.tsv"
opts$d <- "./output/metadata_sample.txt"
opts$i <- "./output/model_anova_total_toddler.txt"
opts$p <- 'bonferroni'
opts$o <- "./output/micro_toddler"
# microbiome data

microbiome_melt <- read.csv(opts$m, sep="\t", header=T) %>%
  mutate(sample_id = factor(sample_id)) %>%
  melt() %>%
  dplyr::filter(value > 0)

eps <- min(microbiome_melt$value)

microbiome <- read.csv(opts$m, sep="\t", header=T) %>%
  mutate(sample_id = factor(sample_id)) %>%
  melt() %>%
  mutate(value = log10(value + eps)) %>%
  dcast(sample_id ~ variable, value.var="value")

microbiome[1:3,1:3]

# genetics data

snps <- read.csv(opts$g, sep="\t", header=T, check.names=F)
colnames(snps)[1] <- "subject_id"
# metadata

metadata <- read.csv(opts$d, "\t", header=T)

# merge data and metadata

df_subjects <- merge(metadata, snps, by="subject_id")

dat <- merge(df_subjects, microbiome, by="sample_id") %>%
  rename(Subject=subject_id, Day=Days) %>%
  mutate(Subject=factor(Subject))
#############################################
# visualisation of significant associations #
#############################################

library(dplyr)
library(reshape2)
library(ggplot2)
library(viridis)
library(scales)
library(lmerTest)

# input

model_stats <- read.csv(opts$i, sep="\t", header=T)

head(model_stats)

# output directory

dir.create(opts$o)

####################
# p-value adjustment

df_sig <- model_stats %>%  filter(!Predictor %in% c("Country","Day","Probiotic")) %>% 
  group_by(Predictor) %>% 
  mutate(P_adj = p.adjust(P, method=opts$p)) %>%
  filter(P_adj < 0.05) %>%
  filter(grepl("PC", Predictor)) %>%
  mutate(Type = case_when(!grepl(":", Predictor) ~ "genetics",
                          grepl(":", Predictor) ~ "genetics:time")) %>%
  mutate(Coe = Coefficient ) %>% 
  mutate(Coefficient = case_when((Coefficient > 0) ~ 1, (Coefficient < 0) ~ -1)) %>%
  mutate(value = Coefficient * -log10(P_adj)) %>%
  mutate(MicroFeat = gsub("lmer\\(", "", Model)) %>%
  mutate(MicroFeat = gsub(" .*", "", MicroFeat)) %>%
  mutate(MicroFeat = gsub("_noname_unclassified", "_?_?", MicroFeat)) %>%
  mutate(MicroFeat = gsub("_unclassified", "_?", MicroFeat))  %>%
  mutate(pcn = gsub("\\:.*", "", Predictor))
df_sig_Toddler <- df_sig %>% mutate(Class = "toddler") %>% mutate(label = "Species")

df_sig1 <- df_sig %>% filter(Type == "genetics")
df_sig2 <- df_sig %>% filter(Type == "genetics:time")

df_pc <- df_sig %>%
  filter(!grepl("\\:", Type))

df_pcxday <- df_sig %>%
  filter(grepl("\\:", Type))

########################
# plot the fixed effects

combined <- list()

for (i in df_pc$Model) {
  print(i)
  
  model <- eval(parse(text=i))
  
  aug <- broom.mixed::augment(model) %>%
    rename(PC=2) %>%
    mutate(Model = i) %>%
    mutate(MicroFeat = gsub("lmer\\(| .*", "", Model)) %>%
    mutate(MicroFeat = gsub("_noname_unclassified", "_?_?", MicroFeat)) %>%
    mutate(MicroFeat = gsub("_unclassified", "_?", MicroFeat)) %>%
    mutate(pcn = gsub("*.\\~  | \\*.*", "", Model)) %>%
    mutate(pcn = gsub(".*P", "P", pcn)) %>%
    mutate(Quartile = ntile(Day, 4)) %>%
    select(Subject, MicroFeat, pcn, Quartile, .fitted, PC)
  
  combined <- rbind(combined, aug)
  
}

# map to add q value to plot

ann_sig <- df_sig %>%
  filter(!grepl(":", Predictor)) %>%
  mutate(pcn = Predictor) %>%
  mutate(P_adj = paste0("q=", scientific(P_adj))) %>%
  select(MicroFeat, pcn, P_adj, Coefficient)

######################
# plot the interaction

combined <- list()

for (i in df_pcxday$Model) {
  
  model <- eval(parse(text=i))
  
  slope <- as.data.frame(coef(model)$Subject) %>%
    tibble::rownames_to_column("Subject") %>%
    mutate(Model = i)
  
  aug <- broom::augment(model) %>%
    rename(RelAb=1, PC=2) %>%
    mutate(Model = i) %>%
    mutate(MicroFeat = gsub("lmer\\(| .*", "", Model)) %>%
    mutate(MicroFeat = gsub("_noname_unclassified", "_?_?", MicroFeat)) %>%
    mutate(MicroFeat = gsub("_unclassified", "_?", MicroFeat)) %>%
    mutate(pcn = gsub("*.\\~  | \\*.*", "", Model)) %>%
    mutate(pcn = gsub(".*P", "P", pcn)) %>%
    mutate(Quartile = ntile(PC, 4)) %>%
    select(Subject, MicroFeat, pcn, Quartile, Day, .fitted, RelAb, PC)
  
  combined <- rbind(combined, aug)	
  
}

ann_sig <- df_sig %>%
  filter(grepl(":", Predictor)) %>%
  mutate(pcn = Predictor) %>%
  mutate(pcn = gsub("\\:.*", "", pcn)) %>%
  mutate(P_adj = paste0("q=", scientific(P_adj))) %>%
  select(MicroFeat, pcn, P_adj, Coefficient)
pdf("./output/micro_toddler//interaction.pdf",10,10)

ggplot(combined, aes(x=Day, y=.fitted)) +
  facet_wrap(~ MicroFeat + pcn, scales="free", ncol=5) +
  geom_smooth(aes(colour=factor(Quartile), fill=factor(Quartile)), method="lm") +
  scale_colour_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete=TRUE) +
  theme_bw() +
  theme(
    plot.title = element_text(size=12.5, face="bold", hjust=0.5),
    axis.title = element_text(size=12.5, face="bold"),
    axis.text = element_text(size=8.75),
    legend.title = element_text(size=12.5, face="bold", vjust=0.75),
    legend.text = element_text(size=12.5),
    strip.text = element_text(size=8.75, face="bold"),
    legend.position = "top",
  ) + #aspect.ratio=0.75
  labs(title="Interaction effect(s)", 
       colour="PC Quartile",
       fill="PC Quartile",
       y=expression(paste(bold(log[10]*" abundance (%)")))) +
  geom_text(data=ann_sig, aes(x=Inf, y=Inf, label=P_adj), 
            hjust=1.5,
            vjust=1.5)

dev.off()