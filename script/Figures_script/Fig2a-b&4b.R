library(dplyr)
library(data.table)
library(nnet)
library(survival)
library(survminer)

set.seed(42)
load("~/ddong/TEDDY_project/output/traj_cluster_results.Rdata")

###Using multinomial logistic regression to calculate propensity score
ps_model <- multinom(Cluster ~ Clinical_Center + Sex+ FDR+
                       delivery  + PC1+PC2+PC3+PC4+PC5+  first_day+
                       Antibiotics +Probiotic+time_to_brstfed_stop +age_startsolid, 
                     data = Cluster800)

# Extract propensity scores for each cluster
ps <- predict(ps_model, type = "probs")  
colnames(ps) <- c("ps_0", "ps_1", "ps_2")

# Adding scores to the clustering data frame
data <- cbind(Cluster800, ps)

#adjusted all covariants
cox_model <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ strata(Cluster) +
                     ps_0+ps_1+ps_2, data = data)



###Figure 2A, smoothed trajectory line plot
pdf("./ddong/TEDDY/output/Fig2a.pdf",5,3.5)
ggplot(metadata_cluster, aes(x= Days,y = bc_dists_init,color = Cluster))+
  geom_smooth(method = "loess") +theme_cowplot()+
  scale_color_manual(values= c("#4c956c","#457b9d","#e26d5c"))
dev.off()


pdf("~/ddong/TEDDY_project/figures/Fig2b.pdf",width = 4.5,height = 4.5)
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

###########################Figure 4B################################################
# Calculate the quantiles of the PCs and divide the data into two parts (50% quantile)
quantiles <- quantile(data$PC1, probs = c(0, 1/2,1))
data$PC1_half <- cut(data$PC1, breaks = quantiles, labels = c("Q1", "Q2" ),include.lowest = T)
quantiles <- quantile(data$PC2, probs = c(0, 1/2,1))
data$PC2_half <- cut(data$PC2, breaks = quantiles, labels = c("Q1", "Q2" ),include.lowest = T)
quantiles <- quantile(data$PC3, probs = c(0, 1/2,1))
data$PC3_half <- cut(data$PC3, breaks = quantiles, labels = c("Q1", "Q2" ),include.lowest = T)
quantiles <- quantile(data$PC4, probs = c(0, 1/2,1))
data$PC4_half <- cut(data$PC4, breaks = quantiles, labels = c("Q1", "Q2" ),include.lowest = T)
quantiles <- quantile(data$PC5, probs = c(0, 1/2,1))
data$PC5_half <- cut(data$PC5, breaks = quantiles, labels = c("Q1", "Q2" ),include.lowest = T)

Cluster800PC3Q1 <- data %>% filter(PC3_half == "Q1")
Cluster800PC3Q2 <- data %>% filter(PC3_half == "Q2")
cox_modelQ1 <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ strata(Cluster) +
                       ps_0+ps_1+ps_2, data = Cluster800PC3Q1)
cox_modelQ2 <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~strata(Cluster) +
                       ps_0+ps_1+ps_2, data = Cluster800PC3Q2)

pdf("~/ddong/TEDDY_project/figures/Fig4b",width = 4.5,height = 4.5)
ggsurvplot(survfit(cox_modelQ1), 
           data = data, 
           fun = "event", 
           palette = c("#4c956c","#457b9d","#e26d5c"),
           # conf.int.alpha = 0.1,
           conf.int = T, # 添加置信区间
           xlab = "Time", 
           ylab = "Cumulative Event Rate",
           ggtheme = theme_cowplot()) # 使用简洁的主题

ggsurvplot(survfit(cox_modelQ1), 
           data = Cluster800, 
           fun = "event", 
           palette = c("#4c956c","#457b9d","#e26d5c"),
           # conf.int.alpha = 0.1,
           conf.int = T, # 添加置信区间
           xlab = "Time", 
           ylab = "Cumulative Event Rate",
           ggtheme = theme_cowplot()) # 使用简洁的主题

dev.off()



