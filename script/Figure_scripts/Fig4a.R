
source("~/script/packages_load.R")
setwd("~/ddong/TEDDY_project")
load("./output/output/Cluster800.Rdata")
meta <- fread("./TEDDY_metadata/metadata_subject.txt")
Cluster800_addprs <-Cluster800 %>% left_join(meta)

quantiles <- quantile(test$PC3, probs = c(0, 1/2,1))
test$PC3_half <- cut(test$PC3, breaks = quantiles, labels = c("Q1", "Q2" ),include.lowest = T)
pc3half <- test %>% select(subject_id,PC3_half) %>% unique()

Cluster800_addprs <- Cluster800_addprs %>% left_join(pc3half)

con <- data.frame()
Cluster800_addprs$PC3_half <- factor(Cluster800_addprs$PC3_half ,levels = c("Q1","Q2"))
Cluster <- Cluster800_addprs 
Cluster1 <- Cluster %>% filter(Cluster ==1)
cox_model <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ Clinical_Center + Sex+ FDR+
                     delivery+ first_day+ PC1+PC2+PC3_half+PC4+PC5+
                     Antibiotics +Probiotic+time_to_brstfed_stop +age_startsolid+
                     Cluster, data = Cluster1)
a <- data.frame(summary(cox_model)$coefficients)
con_99 <- data.frame(confint(cox_model, level = 0.99))
con_95 <- data.frame(confint(cox_model, level = 0.95))
con_1 <- data.frame(coef = a["PC3_halfQ2","coef"] ,down99 =con_99["PC3_halfQ2",1],up99 = con_99["PC3_halfQ2",2],
                    down95 =con_95["PC3_halfQ2",1],up95 = con_95["PC3_halfQ2",2]) %>% mutate(type = "Q2", Cluster =1)



#con1 <- data.frame(t(a["HLA2",]))%>% mutate(type = "Other", Cluster =1)
con <- rbind(con,con_1)

Cluster2 <- Cluster %>% filter(Cluster %in% c(1,3)) %>% filter(PC3_half == "Q1")
cox_model <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ Clinical_Center + Sex+ FDR+
                     delivery+ first_day+ PC1+PC2+PC3_half+PC4+PC5+
                     Antibiotics +Probiotic+time_to_brstfed_stop +age_startsolid+
                     Cluster, data = Cluster2)
a <- data.frame(summary(cox_model)$coefficients)
con_99 <- data.frame(confint(cox_model, level = 0.99))
con_95 <- data.frame(confint(cox_model, level = 0.95))
con_1 <- data.frame(coef = a["Cluster3","coef"] ,down99 =con_99["Cluster3",1],up99 = con_99["Cluster3",2],
                    down95 =con_95["Cluster3",1],up95 = con_95["Cluster3",2]) %>% mutate(type = "Q1", Cluster =3)
con <- rbind(con,con_1)

# a <- summary(cox_model)$conf.int
# con1 <- data.frame(t(a["Cluster2",])) %>% mutate(type = "DR4/DR4", Cluster =2)
# con <- rbind(con,con1)
# print(ggforest(cox_model))

Cluster2 <- Cluster %>% filter((Cluster %in% c( 3) & PC3_half == "Q2") |(Cluster == 1 & PC3_half == "Q1"))

cox_model <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ Clinical_Center + Sex+ FDR+
                     delivery+ first_day+ PC1+PC2+PC3_half+PC4+PC5+
                     Antibiotics +Probiotic+time_to_brstfed_stop +age_startsolid+
                     Cluster, data = Cluster2)
print(ggforest(cox_model))
a <- data.frame(summary(cox_model)$coefficients)
con_99 <- data.frame(confint(cox_model, level = 0.99))
con_95 <- data.frame(confint(cox_model, level = 0.95))
con_1 <- data.frame(coef = a["PC3_halfQ2","coef"] ,down99 =con_99["PC3_halfQ2",1],up99 = con_99["PC3_halfQ2",2],
                    down95 =con_95["PC3_halfQ2",1],up95 = con_95["PC3_halfQ2",2]) %>% mutate(type = "Q2", Cluster =3)

con <- rbind(con,con_1)
# coef = a["HLA2","coef"] ,down99 =con_99["HLA2",1],up99 = con_99["HLA2",2],
# down99 =con_95["Cluster2",1],up99 = con_95["HLA2",2]) %>% mutate(type = "Other", Cluster =1)

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
                     delivery+ first_day+ PC1+PC2+PC3_half+PC4+PC5+
                     Antibiotics +Probiotic+time_to_brstfed_stop +age_startsolid+
                     Cluster, data = Cluster_all)
print(ggforest(cox_model))
a <- summary(cox_model)
p <- signif(a$coefficients["PC3_halfQ2:Cluster3",5],3)

plot2 <- ggplot(con,aes(x = Cluster))+
  #geom_ribbon(aes(ymax = up99,ymin = down99,fill = type),alpha = 0.5)+
  geom_ribbon(aes(ymax = up95,ymin = down95,fill = type),alpha = 0.3)+
  geom_line(aes(group = type,y = coef,color = type))+   
  scale_fill_manual(values = c("#a80022","#0053a8")) + 
  scale_color_manual(values = c("#a80022","#0053a8"))+theme_cowplot() +
  scale_x_continuous(breaks = c(1, 3), labels = c("Healthy maturation", "Incomplete maturation"))+
  ylab("Log(Hazard Ratio)") +
  annotate("text",x = 1.3, y = 0, label = paste("pvalue = ",p,sep =""),size = 4) 


#####Cluster2 vs Cluster1
con <- data.frame()
Cluster1 <- Cluster %>% filter(Cluster ==1)
cox_model <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ Clinical_Center + Sex+ FDR+
                     delivery+ first_day+ PC1+PC2+PC3_half+PC4+PC5+
                     Antibiotics +Probiotic+time_to_brstfed_stop +age_startsolid+
                     Cluster, data = Cluster1)
a <- data.frame(summary(cox_model)$coefficients)
con_99 <- data.frame(confint(cox_model, level = 0.99))
con_95 <- data.frame(confint(cox_model, level = 0.95))
con_1 <- data.frame(coef = a["PC3_halfQ2","coef"] ,down99 =con_99["PC3_halfQ2",1],up99 = con_99["PC3_halfQ2",2],
                    down95 =con_95["PC3_halfQ2",1],up95 = con_95["PC3_halfQ2",2]) %>% mutate(type = "Q2", Cluster =1)



#con1 <- data.frame(t(a["HLA2",]))%>% mutate(type = "Other", Cluster =1)
con <- rbind(con,con_1)

Cluster2 <- Cluster %>% filter(Cluster %in% c(1,2)) %>% filter(PC3_half == "Q1")
cox_model <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ Clinical_Center + Sex+ FDR+
                     delivery+ first_day+ PC1+PC2+PC3_half+PC4+PC5+
                     Antibiotics +Probiotic+time_to_brstfed_stop +age_startsolid+
                     Cluster,data = Cluster2)
a <- data.frame(summary(cox_model)$coefficients)
con_99 <- data.frame(confint(cox_model, level = 0.99))
con_95 <- data.frame(confint(cox_model, level = 0.95))
con_1 <- data.frame(coef = a["Cluster2","coef"] ,down99 =con_99["Cluster2",1],up99 = con_99["Cluster2",2],
                    down95 =con_95["Cluster2",1],up95 = con_95["Cluster2",2]) %>% mutate(type = "Q1", Cluster =2)
con <- rbind(con,con_1)

# a <- summary(cox_model)$conf.int
# con1 <- data.frame(t(a["Cluster2",])) %>% mutate(type = "DR4/DR4", Cluster =2)
# con <- rbind(con,con1)
# print(ggforest(cox_model))

Cluster2 <- Cluster %>% filter((Cluster %in% c( 2) & PC3_half == "Q2") |(Cluster == 1 & PC3_half == "Q1"))

cox_model <- coxph(Surv(AB_T1D_time_to0,AB_T1D) ~ Clinical_Center + Sex+ FDR+
                     delivery+ first_day+ PC1+PC2+PC3_half+PC4+PC5+
                     Antibiotics +Probiotic+time_to_brstfed_stop +age_startsolid+
                     Cluster, data = Cluster2)
print(ggforest(cox_model))
a <- data.frame(summary(cox_model)$coefficients)
con_99 <- data.frame(confint(cox_model, level = 0.99))
con_95 <- data.frame(confint(cox_model, level = 0.95))
con_1 <- data.frame(coef = a["PC3_halfQ2","coef"] ,down99 =con_99["PC3_halfQ2",1],up99 = con_99["PC3_halfQ2",2],
                    down95 =con_95["PC3_halfQ2",1],up95 = con_95["PC3_halfQ2",2]) %>% mutate(type = "Q2", Cluster =2)

con <- rbind(con,con_1)
# coef = a["HLA2","coef"] ,down99 =con_99["HLA2",1],up99 = con_99["HLA2",2],
# down99 =con_95["Cluster2",1],up99 = con_95["HLA2",2]) %>% mutate(type = "Other", Cluster =1)

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
                     delivery+ first_day+ PC1+PC2+PC3_half+PC4+PC5+
                     Antibiotics +Probiotic+time_to_brstfed_stop +age_startsolid+
                     Cluster,data = Cluster_all)
print(ggforest(cox_model))
a <- summary(cox_model)
p <- signif(a$coefficients["PC3_halfQ2:Cluster2",5],3)

plot1 <- ggplot(con,aes(x = Cluster))+
  #geom_ribbon(aes(ymax = up99,ymin = down99,fill = type),alpha = 0.5)+
  geom_ribbon(aes(ymax = up95,ymin = down95,fill = type),alpha = 0.3)+
  geom_line(aes(group = type,y = coef,color = type))+   
  scale_fill_manual(values = c("#a80022","#0053a8")) + 
  scale_color_manual(values = c("#a80022","#0053a8"))+theme_cowplot() +
  scale_x_continuous(breaks = c(1, 2), labels = c("Healthy maturation", "Catch-up maturation"))+
  ylab("Log(Hazard Ratio)") +
  annotate("text",x = 1.3, y = 0, label = paste("pvalue = ",p,sep =""),size = 4) 
pdf("./output/interaction_PC3_half.pdf",12,5)
plot_grid(plot1,plot2,nrow = 1)
dev.off()
