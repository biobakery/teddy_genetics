library(ggplot2)
library(ggside)
library(scales)
library(dplyr)
library(data.table)
library(stringr)

### Function to shorten GO pathway names
extract_words <- function(sentence) {
  # Split the sentence into words
  words <- unlist(str_split(sentence, "\\s+"))
  # Extract the first three words
  first_three <- paste(words[1:min(3, length(words))], collapse = " ")
  # Extract the last three words
  last_three <- paste(words[max(1, length(words) - 2):length(words)], collapse = " ")
  return(paste(first_three, "...", last_three, sep = " "))
}

##--------------------------------------------------------------------------------
## Fig. 4d: Genetic PC GSEA results
##--------------------------------------------------------------------------------

### Data load: GSEA results
load("./output/snpfgsea_allmap.Rdata")
all_sig <- all_snpgsea %>% filter(padj< 0.1) 
all_sig <- all_sig %>% arrange(PC,pval)

### Plot enriched pathways unique to Genetic PC2 and/or PC3 
all_sig_PC23 <- all_sig %>% filter(PC %in% c("PC2","PC3"))
all_sig_nonPC23 <- all_sig %>% filter(PC %in% c("PC1","PC4","PC5"))
all_sig_PC23 <- all_sig_PC23 %>% mutate(unique = !pathway %in% all_sig_nonPC23$pathway )
group <- fread("./output/PC3_gsea_group.txt",header = F)
colnames(group ) <- c("pathway","Group")
all_sig_PC23 <- all_sig_PC23 %>% left_join(group)
all_sig_PC23$Group <- factor(all_sig_PC23$Group,levels = c("Antigen Processing and Presentation",
                                                           "Antiviral Immune Response",
                                                           "Microbial Defense",
                                                           "Immune System Regulation",
                                                           "Protein Modification and Degradation",
                                                           "Cellular Processes and Membrane Dynamics",
                                                           "DNA Damage Response"))
colors <- c("#bc4749",'#e09f3e','#a4ac86','#87bba2','#457b9d',"#5e548e","grey")

all_sig_PC23$GO_term <- unlist(lapply(all_sig_PC23$pathway, extract_words))
all_sig_PC23 <- all_sig_PC23 %>% arrange(PC,desc(Group),padj)
all_sig_PC23$GO_term <- factor(all_sig_PC23$GO_term,levels = unique(all_sig_PC23$GO_term))

pdf("./figures/Fig4d_GSEA.pdf",10.5,5.5)
ggplot(all_sig_PC23,aes(x = ES, y = GO_term)) +
  geom_point(aes(color = padj,size = size)) + 
  labs(title = "GSEA_PC3", x = "enrichmentScore", y = "") + 
  theme_bw()+scale_color_gradientn(colors = rev(c("#ade8f4","#0077b6","#03045e")),
                                   values = rescale(c(0,0.02,0.1)),
                                   breaks=c(0.025,0.05,0.075),
                                   limits=c(0,0.1))+
  coord_cartesian(clip = "off") +facet_grid(~PC) +
  geom_ysidetile(data = all_sig_PC23,aes(x = 1,yfill = Group))+
  scale_yfill_manual(values = colors)
dev.off()
