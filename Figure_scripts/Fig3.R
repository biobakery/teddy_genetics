library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)
library(viridis)

setwd("~/ddong/TEDDY_project")
##---------------------------------------------------
## Fig. 3a: BCAA, Galactose, Vitamin B, AAA 
##---------------------------------------------------

### Data Loading: QCed EC, EC name
ec.df <- fread("./output/ec_all_qc_1e6.txt")
ec_name <- fread("./output/ec_name.txt")

### Loading Manually Selected EC
ec_select <- fread("./output/select_EC.csv",header = T)
ec_select <- ec_select %>% left_join(ec_name)
ec_select <- ec_select %>% mutate(feature = paste("X",EC,sep = ""))
binnames <- c(sprintf("[%d,%d]", 100, 200), 
              sprintf("(%d,%d]", seq(200, 700, by = 100), 
                      seq(300, 800, by = 100)))

### Loading MaAslin results
data <- list()
all_EC.df <- data.frame()
setwd("./output/maaslin_EC_0323/")
for(i in 1:length(binnames)){
  b_n <- binnames[i]
  data[[i]] <- fread(paste("./", paste("Cluster",b_n,sep = "_"),"/all_results.tsv",sep = "") )
  data[[i]]  <- data[[i]]  %>% filter(metadata == "Cluster")  %>%
    mutate(Type = paste("Cluster1_Cluster",value,sep =""),
           bin = binnames[i])
  all_EC.df <- rbind(all_EC.df,data[[i]])
}
data2 <-list()
for(i in 1:length(binnames)){
  b_n <- binnames[i]
  data2[[i]] <- fread(paste("./", paste("Cluster_base2_",b_n,sep = "_"),"/all_results.tsv",sep = "") )
  data2[[i]]  <- data2[[i]]  %>% filter(metadata == "Cluster")  %>%
    mutate(Type = paste("Cluster2_Cluster",value,sep =""),
           bin = binnames[i])
  all_EC.df <- rbind(all_EC.df,data2[[i]] )
}
table(all_EC.df$Type)
ec_name <- ec_name %>% mutate(feature = paste("X",EC,sep = ""))
all_EC.df <- all_EC.df %>% left_join(ec_name)
all_EC.df <- all_EC.df %>% filter(Type != "Cluster2_Cluster1")
# save(all_EC.df,file = "./output/allEC_result.Rdata")

### Select significant results based on q-value threshold
sig <- all_EC.df %>% filter(qval < 0.25)
alldata_select_EC <- all_EC.df %>% filter(EC %in% sig$EC)
alldata_select_EC$bin  <- factor(alldata_select_EC$bin,levels = binnames)
alldata_select_EC <-alldata_select_EC %>% mutate(sig = case_when(qval < 0.001 ~ "***",
                                                                 qval < 0.05 ~ "**",
                                                                 qval < 0.25 ~ "*"))


### Load KEGG Orthology (KO) data
KO <- fread("~/Downloads/KO_03202022.csv")
KO <- KO %>% mutate(EC_number = gsub("EC","",ec)) %>% mutate(feature = paste("X",EC_number,sep = ""))
### Merge EC selection with KO data
sig_EC <- alldata_select_EC %>% select(feature,EC_name,EC) %>% unique() #367
sig_EC <- sig_EC %>% left_join(KO)
sig_EC <- sig_EC %>% select(cat1,cat2,cat3,cat4,ko, feature,EC,ec,EC_number,one_of(colnames(sig_EC)))


### Select interested KO
KO_select <- KO %>% filter(cat2 %in% c("09105 Amino acid metabolism", 
                                       "09106 Metabolism of other amino acids",
                                       "09108 Metabolism of cofactors and vitamins") |
                             cat3 %in% c("00052 Galactose metabolism [PATH:ko00052]"))

sig_EC <- alldata_select_EC %>% select(feature,EC_name,EC) %>% unique() #367
sig_EC <- sig_EC %>% left_join(KO)
sig_EC <- sig_EC %>% select(cat1,cat2,cat3,cat4,ko, feature,EC,ec,EC_number,one_of(colnames(sig_EC)))

sig_EC <- sig_EC %>% left_join(ec_select)
sig_EC_select <- sig_EC %>% filter(cat2 %in% c("09105 Amino acid metabolism", 
                                               "09106 Metabolism of other amino acids",
                                               "09108 Metabolism of cofactors and vitamins") |
                                     cat3 %in% c("00052 Galactose metabolism [PATH:ko00052]")|!is.na(EC_group_top))
sig_EC_select <- sig_EC_select %>% mutate(paste = paste(cat1,cat2,cat3,EC,sep =""))
sig_EC_select <- sig_EC_select %>% filter(!duplicated(paste))

### Select interested pathway
select_pathway <- c(
  "00260 Glycine, serine and threonine metabolism [PATH:ko00260]",           
  "00270 Cysteine and methionine metabolism [PATH:ko00270]",                
  "00280 Valine, leucine and isoleucine degradation [PATH:ko00280]",     
  "00290 Valine, leucine and isoleucine biosynthesis [PATH:ko00290]",       
  "00300 Lysine biosynthesis [PATH:ko00300]",                           
  "00310 Lysine degradation [PATH:ko00310]",                              
  "00340 Histidine metabolism [PATH:ko00340]"  ,                             
  "00360 Phenylalanine metabolism [PATH:ko00360]" ,                          
  "00380 Tryptophan metabolism [PATH:ko00380]",                            
  "00400 Phenylalanine, tyrosine and tryptophan biosynthesis [PATH:ko00400]",
  "00670 One carbon pool by folate [PATH:ko00670]",                
  "00730 Thiamine metabolism [PATH:ko00730]",                          
  "00740 Riboflavin metabolism [PATH:ko00740]",                         
  "00750 Vitamin B6 metabolism [PATH:ko00750]",
  "00760 Nicotinate and nicotinamide metabolism [PATH:ko00760]",           
  "00770 Pantothenate and CoA biosynthesis [PATH:ko00770]",               
  "00780 Biotin metabolism [PATH:ko00780]",                          
  "00790 Folate biosynthesis [PATH:ko00790]",
  "00052 Galactose metabolism [PATH:ko00052]")

KO_select <- KO %>% filter(cat3 %in% select_pathway)

ec_abundance_melt.df <- melt(ec.df,"sample_id")
colnames(ec_abundance_melt.df) <- c("sample_id","feature","abundance")

load("./output/traj_cluster_results.Rdata")
cluster_bin.df <- metadata_cluster %>% dplyr::select(subject_id,sample_id,bins,Cluster,Days)
ec_abundance_melt.df <- ec_abundance_melt.df  %>% inner_join(cluster_bin.df)
ec_abundance_melt.df <-ec_abundance_melt.df %>% left_join(ec_select)
ec_abundance_melt.df <- ec_abundance_melt.df %>% mutate(hit =1)
ec_abundance_melt.df <- ec_abundance_melt.df %>% inner_join(KO_select)
esentail_biosynthesis <- c("00290 Valine, leucine and isoleucine biosynthesis [PATH:ko00290]",
                           "00300 Lysine biosynthesis [PATH:ko00300]",
                           "00400 Phenylalanine, tyrosine and tryptophan biosynthesis [PATH:ko00400]")
metabolism <- setdiff(select_pathway[1:10],esentail_biosynthesis)

ec_abundance_melt.df <- ec_abundance_melt.df %>% mutate(EC_top = case_when(cat3 %in% select_pathway[1:10] ~ "Essential Amino acid metabolism",
                                                                           cat3 %in% select_pathway[11:18] ~ "VitaminB and co-factor",
                                                                           cat3 %in% select_pathway[19] ~ "Galactose"))


ec_abundance_melt.df <- ec_abundance_melt.df %>% mutate(EC_top = case_when(cat3 %in% c("00290 Valine, leucine and isoleucine biosynthesis [PATH:ko00290]")~ "BACC biosynthes",
                                                                           cat3 %in% c("00400 Phenylalanine, tyrosine and tryptophan biosynthesis [PATH:ko00400]") ~ "AAA biosynthes",
                                                                           cat3 %in% select_pathway[11:18] ~ "VitaminB and co-factor",
                                                                           cat3 %in% select_pathway[19] ~ "Galactose"))

ec_abundance_melt.df <- ec_abundance_melt.df %>% group_by(feature) %>% mutate(scale = scale(abundance))
ec_abundance_melt_sum_top.df <- ec_abundance_melt.df %>% ungroup() %>%  
  mutate(paste = paste(EC_number,bins,sample_id,EC_top)) %>% filter(!duplicated(paste)) %>% 
  group_by(sample_id,subject_id,bins,Cluster,EC_top,Days) %>% summarise(sum_group = sum(abundance),hit =sum(hit))

pdf("./figures/Fig3az",17,4)
ggplot(ec_abundance_melt_sum_top.df,aes(x = Days,y = sum_group,color= Cluster))+ 
  geom_smooth(method = "gam")+facet_wrap(~ EC_top,scales = "free",ncol = 5) +theme_cowplot()+
  scale_color_manual(values= c("#4c956c","#457b9d","#e26d5c"))
dev.off()

save.image('./output/Fig3a.Rdata')


##---------------------------------------------------
## Fig. 3b: Heatmap of Select ECs' Relative Abundance
##---------------------------------------------------

ec_abundance_melt_select_mean <- ec_abundance_melt.df %>% 
  group_by(feature, EC_top,EC_number,bins,Cluster,name,cat3) %>% summarise(mean = mean(scale,na.rm = T))
ec_abundance_melt_select_mean <- ec_abundance_melt_select_mean %>% group_by(feature,bins) %>% mutate(rank = rank(mean))
ec_abundance_melt_select_mean <- ec_abundance_melt_select_mean %>% left_join(ec_name,by = c("feature","EC_number" = "EC"))
ec_abundance_melt_select_mean <- ec_abundance_melt_select_mean %>% mutate(paste= paste(feature,bins,Cluster) ) %>% 
  dplyr::filter(!duplicated(paste))
ec_abundance_melt_select_mean <- ec_abundance_melt_select_mean[!duplicated(ec_abundance_melt_select_mean$paste),]

### Cluster ECs related to Vitamin, Galactose, AAA, and BCAA
### Vitamin related 
vitamin_ec.df <- ec_abundance_melt_select_mean %>% filter(EC_top == "VitaminB and co-factor") %>% filter(Cluster == 3)
vitamin_ec.df$value <- vitamin_ec.df$mean  
vitamin_ec_dcast.df <- dcast(vitamin_ec.df,"EC_number~bins")
hc.cols <- hclust(dist(vitamin_ec_dcast.df[,-1]))
vitamin_ec_dcast.df <- vitamin_ec_dcast.df %>% left_join(ec_name,by = c("EC_number" = "EC"))
cluster_vt <- as.character(vitamin_ec_dcast.df$EC_number[hc.cols$order])
cluster_vt_ecname <- as.character(vitamin_ec_dcast.df$EC_name[hc.cols$order])

### Galactose related 
Galactose_ec.df <- ec_abundance_melt_select_mean %>% filter(EC_top == "Galactose") %>% filter(Cluster == 3)
Galactose_ec.df$value <- Galactose_ec.df$mean  
Galactose_ec_dcast.df <- dcast(Galactose_ec.df,"EC_number~bins")
hc.cols <- hclust(dist(Galactose_ec_dcast.df[,-1]))
Galactose_ec_dcast.df <- Galactose_ec_dcast.df %>% left_join(ec_name,by = c("EC_number" = "EC"))
cluster_gl <- as.character(Galactose_ec_dcast.df$EC_number[hc.cols$order])
cluster_gl_ecname <- as.character(Galactose_ec_dcast.df$EC_name[hc.cols$order])

### BCAA related 
BACC_ec.df <- ec_abundance_melt_select_mean %>% filter(EC_top == "BACC biosynthes") %>% filter(Cluster == 3)
BACC_ec.df$value <- BACC_ec.df$mean  
BACC_ec_dcast.df <- dcast(BACC_ec.df,"EC_number~bins")
hc.cols <- hclust(dist(BACC_ec_dcast.df[,-1]))
cluster_bacc <- as.character(BACC_ec_dcast.df$EC_number[hc.cols$order])
BACC_ec_dcast.df <- BACC_ec_dcast.df %>% left_join(ec_name,by = c("EC_number" = "EC"))
cluster_bacc_ecname <- as.character(BACC_ec_dcast.df$EC_name[hc.cols$order])

### AAA related 
ec_abundance_melt_select_mean <- ec_abundance_melt_select_mean %>% mutate(paste= paste(feature,bins,Cluster) ) %>% filter(!duplicated(paste))
AAA_ec.df <- ec_abundance_melt_select_mean %>% filter(EC_top == "AAA biosynthes") %>% filter(Cluster == 3)
AAA_ec.df$value <- AAA_ec.df$mean  
AAA_ec_dcast.df <- dcast(AAA_ec.df,"EC_number~bins")
#vitamin_ec_dcast.df$bin  <- factor(vitamin_ec_dcast.df$bin,levels = binnames)
hc.cols <- hclust(dist(AAA_ec_dcast.df[,-1]))
cluster_AAA <- as.character(AAA_ec_dcast.df$EC_number[hc.cols$order])
AAA_ec_dcast.df <- AAA_ec_dcast.df%>% left_join(ec_name,by = c("EC_number" = "EC"))
cluster_AAA_ecname <- as.character(AAA_ec_dcast.df$EC_name[hc.cols$order])

### combine for groups 
level <- c(cluster_bacc,cluster_AAA,cluster_gl,cluster_vt)
ec_abundance_melt_select_mean$EC_number <- factor(ec_abundance_melt_select_mean$EC_number,levels = unique(rev(level)))
combine <- data.frame(EC =c(cluster_bacc,cluster_AAA,cluster_gl,cluster_vt),
                      EC_name = c(cluster_bacc_ecname,cluster_AAA_ecname,cluster_gl_ecname,cluster_vt_ecname))
combine <- combine %>% filter(combine$EC %in% sig$EC)
type <- ec_abundance_melt_select_mean %>% ungroup %>% select(EC_number,cat3,EC_top) %>% unique()
combine <- combine %>% left_join(type,by = c("EC" = "EC_number"))

ec_abundance_melt_select_mean_select <- ec_abundance_melt_select_mean_select %>% left_join(ec_name)
ec_abundance_melt_select_mean_select <- ec_abundance_melt_select_mean_select %>% 
  mutate(EC_ID_Name = paste(EC_number,EC_name,sep ="_"))
ec_abundance_melt_select_mean_select <- ec_abundance_melt_select_mean_select %>% arrange(EC_number)
ec_abundance_melt_select_mean_select$EC_ID_Name <- factor(ec_abundance_melt_select_mean_select$EC_ID_Name,
                                                          levels = unique(ec_abundance_melt_select_mean_select$EC_ID_Name))

pdf("./figures/Fig3b.pdf",15,9)
ggplot(ec_abundance_melt_select_mean_select,aes(x = bins,y = EC_ID_Name,fill = mean)) +geom_tile()+
  facet_wrap(~ Cluster,ncol = 8) +   geom_ysidetile(aes(x = "EC_group",yfill = EC_top))+
  geom_ysidetile(aes(x = "EC_group",yfill = EC_top))+ theme_cowplot()+
  scale_x_discrete(labels =c(seq(200,800,100))) +xlab("Bin(Days)")+
  scale_fill_viridis(option = "D")
dev.off()  

save.image("./output/Fig3.Rdata")