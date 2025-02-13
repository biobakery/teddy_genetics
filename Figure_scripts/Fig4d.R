##shorten the GO pathway
extract_words <- function(sentence) {
  # 使用正则表达式分割句子成单词
  words <- unlist(str_split(sentence, "\\s+"))
  
  # 提取前三个单词
  first_three <- paste(words[1:min(3, length(words))], collapse = " ")
  
  # 提取后三个单词
  last_three <- paste(words[max(1, length(words) - 1):length(words)], collapse = " ")
  
  return(paste(first_three,"...",last_three,sep = " "))
}

source("~/script/packages_load.R")
load("/TEDDY_project/output/snpfgsea_allmap.Rdata")
all_sig <- all_snpgsea %>% filter(padj< 0.1) 
#all_sig <- all %>% filter(p.adjust< 0.25) #%>% filter(PC == "sub32")
all_sig <- all_sig %>% arrange(PC,pval)
all_sig$PC <- factor(all_sig$PC,levels =  c(paste("PC",1:5,sep = ""),"PC3-PC2","PC3-PC4","PC3-PC5"))
all_sig$pathway <- factor(all_sig$pathway,levels = unique(all_sig$pathway))
all_sig_PC3 <- all_sig %>% filter(PC %in% c("PC2","PC3"))
all_sig_nonPC3 <- all_sig %>% filter(PC %in% c("PC1","PC4","PC5"))
# all_sig_PC3 <- all_sig %>% filter(PC %in% c("PC3","PC3"))
# all_sig_nonPC3 <- all_sig %>% filter(PC %in% c("PC1","PC2","PC4","PC5"))
all_sig_PC3 <- all_sig_PC3 %>% mutate(unique = !pathway %in% all_sig_nonPC3$pathway )
group <- fread("~/ddong/TEDDY_project/output/PC3_gsea_group.txt",header = F)
colnames(group ) <- c("pathway","Group")
all_sig_PC3 <- all_sig_PC3 %>% left_join(group)
all_sig_PC3$Group <- factor(all_sig_PC3$Group,levels = c("Antigen Processing and Presentation",
                                                         "Antiviral Immune Response",
                                                         "Microbial Defense",
                                                         "Immune System Regulation",
                                                         "Protein Modification and Degradation",
                                                         "Cellular Processes and Membrane Dynamics",
                                                         "DNA Damage Response"))
colors <- c("#bc4749",'#e09f3e','#a4ac86','#87bba2','#457b9d',"#5e548e","grey")

all_sig_PC3$GO_term <- unlist(lapply(all_sig_PC3$pathway, extract_words))
all_sig_PC3 <- all_sig_PC3 %>% arrange(PC,desc(Group),padj)
all_sig_PC3$GO_term <- factor(all_sig_PC3$GO_term,levels = unique(all_sig_PC3$GO_term))
library(scales)
pdf("~/ddong/TEDDY_project/output/PC_gseas_PC2&3.pdf",10.5,5.5)
ggplot(all_sig_PC3,aes(x = ES, y = GO_term)) +
  geom_point(aes(color = padj,size = size)) + #scale_color_gradient2(low = "#0053a8",mid = "white",high = "#a80022")+
  labs(title = "GSEA_PC3", x = "enrichmentScore", y = "") + 
  theme_bw()+scale_color_gradientn(colors = rev(c("#ade8f4","#0077b6","#03045e")),
                                   values = rescale(c(0,0.02,0.1)),
                                   breaks=c(0.025,0.05,0.075),
                                   limits=c(0,0.1))+
  # theme_minimal() + 
  coord_cartesian(clip = "off") +facet_grid(~PC) +
  geom_ysidetile(data = all_sig_PC3,aes(x = 1,yfill = Group))+
  scale_yfill_manual(values = colors)
dev.off()
