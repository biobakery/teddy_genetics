library(dplyr)
library(data.table)
library(ggplot2)
library(nnet)
library(survival)
library(survminer)
library(tidyverse)
library(cowplot)


##---------------------------------------------------
## Fig. 2a: Line Plot of Three Maturational Patterns
##---------------------------------------------------

### Load trajectory cluster results (subject level)
load("~/ddong/TEDDY_project/output/traj_cluster_results.Rdata")

### Load bray-curtis distance results. (sample level)
bc_distance <- "./output/bc_dists_init.tsv"
bc_distance <- data.table::fread(bc_distance, sep="\t", header=T, check.names=F) %>%
  rename(sample_id=1) %>%
  column_to_rownames("sample_id") %>%
  t(.) %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  #setNames(paste0("microbiome_", names(.))) %>%
  rename(sample_id=1)

###add cluster results to sample level. 
bc_distance <- bc_distance %>% inner_join(Cluster_results,by = c("Subject" = "subject_id"))


pdf("./figures/Fig2a_lineplot.pdf",5,3.5)
ggplot(bc_distance, aes(x= Day,y = bc_dists_init,color = Cluster))+
  geom_smooth(method = "loess") +theme_cowplot()+
  scale_color_manual(values = c( "1" = "#4c956c",
                                 "2" = "#457b9d",
                                 "3" = "#e26d5c"),
                     labels =  c("1" = "Healthy Maturation", 
                                 "2" = "Catch-up Maturation", 
                                 "3" = "Incomplete Maturation"))+
  xlim(min(bc_distance$Day),800) 
dev.off()



##---------------------------------------------------
## Fig. 2b: Cumulative Event Rate plot
##---------------------------------------------------
###Using multinomial logistic regression to calculate propensity score
set.seed(42)
ps_model <- multinom(Cluster ~ Clinical_Center + Sex + FDR + delivery_class +
                       PC1 + PC2 + PC3 + PC4 + PC5+
                       time_to_brstfed_stop + age_startsolid+
                       Antibiotics + Probiotic, 
                     data = Cluster_results)

# Extract propensity scores for each cluster
ps <- predict(ps_model, type = "probs")  
colnames(ps) <- c("ps_0", "ps_1", "ps_2")

# Adding scores to the clustering data frame
data <- cbind(Cluster_results, ps)

#Adjusted all covariants
cox_model <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ strata(Cluster) +
                     ps_0+ps_1+ps_2, data = data)

pdf("./figures/Fig2b_cumulative.pdf",width = 4.5,height = 4.5)
ggsurvplot(survfit(cox_model), 
           data = data, 
           fun = "event", 
           palette = c("#4c956c","#457b9d","#e26d5c"),
           # conf.int.alpha = 0.1,
           conf.int = T, # 添加置信区间
           xlab = "Time", 
           ylab = "Cumulative Event Rate",
           ggtheme = theme_cowplot()) # 使用简洁的主题
dev.off()
save.image("./output/Fig2a_b.Rdata")

