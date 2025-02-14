library(data.table)   
library(dplyr)       
library(survival)    
library(ggplot2)
library(cowplot)

setwd("~/ddong/TEDDY_project") 

##---------------------------------------------------
## Fig4a Interaction plot
##---------------------------------------------------

### Data Loading: Clsuter results, Metadata
load("./output/traj_cluster_results.Rdata")  
meta <- fread("./TEDDY_metadata/metadata_subject.txt") 

# Merge cluster data with metadata
Cluster_results <- Cluster_results %>% left_join(meta)

# Create PC3 grouping by median split
quantiles <- quantile(Cluster_results$PC3, probs = c(0, 0.5, 1))  # Compute median
Cluster_results$PC3_half <- cut(Cluster_results$PC3, breaks = quantiles, labels = c("Q1", "Q2"), include.lowest = TRUE)
Cluster_results$PC3_half <- factor(Cluster_results$PC3_half, levels = c("Q1", "Q2"))
Cluster <- Cluster_results

### Fig4a Left >> Healthy VS Catch-up   
### Only Cluster 1
con <- data.frame()
Cluster1 <- Cluster %>% filter(Cluster ==1)
cox_model <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ Clinical_Center + Sex+ FDR+
                     delivery+ PC1+PC2+PC3_half+PC4+PC5+
                     Antibiotics +Probiotic+time_to_brstfed_stop +age_startsolid+
                     Cluster, data = Cluster1)

coef_summary <- data.frame(summary(cox_model)$coefficients)
con_99 <- data.frame(confint(cox_model, level = 0.99))
con_95 <- data.frame(confint(cox_model, level = 0.95))
con_1 <- data.frame(coef = coef_summary["PC3_halfQ2","coef"] ,
                    down99 =con_99["PC3_halfQ2",1],
                    up99 = con_99["PC3_halfQ2",2],
                    down95 =con_95["PC3_halfQ2",1],
                    up95 = con_95["PC3_halfQ2",2]) %>% 
  mutate(type = "Q2", Cluster =1)
con <- rbind(con,con_1)

#### Cluster 1 (Healthy) vs 2 (Catch-up) , PC3 Q1
Cluster2 <- Cluster %>% filter(Cluster %in% c(1,2)) %>% filter(PC3_half == "Q1")
cox_model <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ Clinical_Center + Sex+ FDR+
                     delivery+ PC1+PC2+PC3_half+PC4+PC5+
                     Antibiotics +Probiotic+time_to_brstfed_stop +age_startsolid+
                     Cluster,data = Cluster2)

coef_summary <- data.frame(summary(cox_model)$coefficients)
con_99 <- data.frame(confint(cox_model, level = 0.99))
con_95 <- data.frame(confint(cox_model, level = 0.95))
con_1 <- data.frame(coef = coef_summary["Cluster2","coef"] ,
                    down99 =con_99["Cluster2",1],
                    up99 = con_99["Cluster2",2],
                    down95 =con_95["Cluster2",1],
                    up95 = con_95["Cluster2",2]) %>%
  mutate(type = "Q1", Cluster =2)
con <- rbind(con,con_1)

######Cluster 2 (Catch-up) Q2 VS Cluster 1 (Healthy) Q1 
Cluster2 <- Cluster %>% filter((Cluster %in% c( 2) & PC3_half == "Q2") |(Cluster == 1 & PC3_half == "Q1"))

cox_model <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ Clinical_Center + Sex+ FDR+
                     delivery+ PC1+PC2+PC3_half+PC4+PC5+
                     Antibiotics +Probiotic+time_to_brstfed_stop +age_startsolid+
                     Cluster, data = Cluster2)

coef_summary <- data.frame(summary(cox_model)$coefficients)
con_99 <- data.frame(confint(cox_model, level = 0.99))
con_95 <- data.frame(confint(cox_model, level = 0.95))
con_1 <- data.frame(coef = coef_summary["PC3_halfQ2","coef"] ,
                    down99 =con_99["PC3_halfQ2",1],
                    up99 = con_99["PC3_halfQ2",2],
                    down95 =con_95["PC3_halfQ2",1],
                    up95 = con_95["PC3_halfQ2",2]) %>% 
  mutate(type = "Q2", Cluster =2)

con <- rbind(con,con_1)

################
Other <- data.frame(coef  =0 ,
                    down99= 0,
                    up99 =0 ,
                    down95 = 0 ,
                    up95 = 0,
                    type = "Q1",
                    Cluster =1)
con <- rbind(con,Other)

# con$exp..coef. <- NULL
# con$Cluster <- factor(con$Cluster)
# con.df <- melt(con)
# con.df <- con.df %>% mutate(class = paste(type,variable))
Cluster_all<- Cluster %>% filter(Cluster %in% c(1,2,3)) 
cox_model <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ Clinical_Center + Sex+ FDR+
                     delivery+PC1+PC2+PC3_half+PC4+PC5+
                     PC1:Cluster + PC2:Cluster +PC3_half:Cluster +PC4:Cluster+PC5:Cluster+
                     Antibiotics +Probiotic+time_to_brstfed_stop +age_startsolid+
                     Cluster,data = Cluster_all)

coef_summary <- summary(cox_model)
p <- signif(coef_summary$coefficients["PC3_halfQ2:Cluster2",5],3)

#### Interaction plot
plot1 <- ggplot(con,aes(x = Cluster))+
  #geom_ribbon(aes(ymax = up99,ymin = down99,fill = type),alpha = 0.5)+
  geom_ribbon(aes(ymax = up95,ymin = down95,fill = type),alpha = 0.3)+
  geom_line(aes(group = type,y = coef,color = type))+   
  scale_fill_manual(values = c("#a80022","#0053a8")) + 
  scale_color_manual(values = c("#a80022","#0053a8"))+theme_cowplot() +
  scale_x_continuous(breaks = c(1, 2), labels = c("Healthy maturation", "Catch-up maturation"))+
  ylab("Log(Hazard Ratio)") +
  annotate("text",x = 1.3, y = 0, label = paste("pvalue = ",p,sep =""),size = 4) 

##########################################################################################
######Fig4a Left >> Healthy VS Incomplete                 
# Calculated relative risk in Cluster 1 (Healthy Maturation pattern) for different PC3 groups
con <- data.frame() # data frame to save all CI
Cluster1 <- Cluster %>% filter(Cluster == 1) 
cox_model <- coxph(Surv(AB_T1D_time_to0, AB_T1D) ~ Clinical_Center + Sex + FDR +
                     delivery + PC1 + PC2 + PC3_half + PC4 + PC5 +
                     Antibiotics + Probiotic + time_to_brstfed_stop + 
                     age_startsolid + Cluster, data = Cluster1)

# Extract and store model coefficients and confidence intervals
coef_summary <- data.frame(summary(cox_model)$coefficients)
con_99 <- data.frame(confint(cox_model, level = 0.99))
con_95 <- data.frame(confint(cox_model, level = 0.95))
con_1 <- data.frame(coef = coef_summary["PC3_halfQ2", "coef"],
                    down99 = con_99["PC3_halfQ2", 1],
                    up99 = con_99["PC3_halfQ2", 2],
                    down95 = con_95["PC3_halfQ2", 1],
                    up95 = con_95["PC3_halfQ2", 2]) %>% 
  mutate(type = "Q2", Cluster = 1)

# Calculated relative risk in Cluster 1 (Healthy Maturation pattern) for different PC3 Q1
Cluster2 <- Cluster %>% filter(Cluster %in% c(1, 3)) %>% filter(PC3_half == "Q1")

cox_model <- coxph(Surv(AB_T1D_time_to0, AB_T1D) ~ Clinical_Center + Sex + FDR +
                     delivery + PC1 + PC2 + PC3_half + PC4 + PC5 +
                     Antibiotics + Probiotic + time_to_brstfed_stop + 
                     age_startsolid + Cluster, data = Cluster2)

coef_summary <- data.frame(summary(cox_model)$coefficients)
con_99 <- data.frame(confint(cox_model, level = 0.99))
con_95 <- data.frame(confint(cox_model, level = 0.95))
con_1 <- data.frame(coef = coef_summary["Cluster3","coef"] ,
                    down99 =con_99["Cluster3",1],
                    up99 = con_99["Cluster3",2],
                    down95 =con_95["Cluster3",1],
                    up95 = con_95["Cluster3",2]) %>%
  mutate(type = "Q1", Cluster =3)
con <- rbind(con,con_1)

# Calculated relative risk in Healthy Q1 VS Incomplete Q3
Cluster2 <- Cluster %>% filter((Cluster %in% c( 3) & PC3_half == "Q2") |(Cluster == 1 & PC3_half == "Q1"))
cox_model <- coxph(Surv(AB_T1D_time_to0, AB_T1D) ~ Clinical_Center + Sex + FDR +
                     delivery + PC1 + PC2 + PC3_half + PC4 + PC5 +
                     Antibiotics + Probiotic + time_to_brstfed_stop + 
                     age_startsolid + Cluster, data = Cluster2)


coef_summary <- data.frame(summary(cox_model)$coefficients)
con_99 <- data.frame(confint(cox_model, level = 0.99))
con_95 <- data.frame(confint(cox_model, level = 0.95))
con_1 <- data.frame(coef = coef_summary["PC3_halfQ2","coef"] ,
                    down99 =con_99["PC3_halfQ2",1]
                    ,up99 = con_99["PC3_halfQ2",2],
                    down95 =con_95["PC3_halfQ2",1],
                    up95 = con_95["PC3_halfQ2",2]) %>% 
  mutate(type = "Q2", Cluster =3)

con <- rbind(con,con_1)

Cluster_all<- Cluster %>% filter(Cluster %in% c(1,2,3)) 
cox_model <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ Clinical_Center + Sex+ FDR+
                     delivery+ first_day+ PC1+PC2+PC3_half+PC4+PC5+
                     PC1:Cluster + PC2:Cluster +PC3_half:Cluster +PC4:Cluster+PC5:Cluster+
                     Antibiotics +Probiotic+time_to_brstfed_stop +age_startsolid+
                     Cluster, data = Cluster_all)

coef_summary <- summary(cox_model)
p <- signif(coef_summary$coefficients["PC3_halfQ2:Cluster3",5],3)

Other <- data.frame(coef  =0 ,
                    down99= 0,
                    up99 =0 ,
                    down95 = 0 ,
                    up95 = 0,
                    type = "Q1",
                    Cluster =1)
con <- rbind(con,Other)

# interaction plot.
plot2 <- ggplot(con,aes(x = Cluster))+
  #geom_ribbon(aes(ymax = up99,ymin = down99,fill = type),alpha = 0.5)+
  geom_ribbon(aes(ymax = up95,ymin = down95,fill = type),alpha = 0.3)+
  geom_line(aes(group = type,y = coef,color = type))+   
  scale_fill_manual(values = c("#a80022","#0053a8")) + 
  scale_color_manual(values = c("#a80022","#0053a8"))+theme_cowplot() +
  scale_x_continuous(breaks = c(1, 3), labels = c("Healthy maturation", "Incomplete maturation"))+
  ylab("Log(Risk Ratio)") +
  annotate("text",x = 1.3, y = 0, label = paste("pvalue = ",p,sep =""),size = 4) 


# Save the combined plot to a PDF
pdf("./figures/Fig4a_interaction_PC3.pdf", 12, 5)
plot_grid(plot1, plot2, nrow = 1)
dev.off()

save.image("./output/Fig4a.Rdata")
