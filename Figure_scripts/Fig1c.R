library(tidyverse)
library(data.table)
library(ggplot2)
library(cowplot)
library(viridis)

setwd("~/ddong/TEDDY_project/")
### Data Load: Genetic SNPs
genetic <- fread("./data/genetics.bed/teddy-final.bim")
colnames(genetic) <- c("Chrom","rsID","cM","Position","REF","ALT")
genetic <- genetic %>% filter(Chrom %in% c(1:22))

### Divide chromosome positions into 5 Mb bins
breaks <- seq(0,250000000,5000000)
genetic$bin <- findInterval(genetic$Position,breaks)
### Summarize the number of variants in each bin
genetic_summary <- genetic %>% group_by(Chrom,bin) %>% count()
genetic_summary$Chrom <- factor(genetic_summary$Chrom,levels = rev(c(c(1:22))))

### Find the maximum bin number for each chromosome
maxbin <- genetic %>% group_by(Chrom) %>% summarise(maxbin =max(bin))

### Creating a dataframe with all possible chromosome-bin combinations
allbin <- data.frame(Chrom = rep(maxbin$Chrom,maxbin$maxbin))
bins <- vector()
for(i in 1:nrow(maxbin)){
  bins <- c(bins,1:maxbin$maxbin[i])
}
allbin$bin <- bins
allbin$Chrom <- factor(allbin$Chrom,levels = rev(c(c(1:22))))
allbin.df <- allbin %>% left_join(genetic_summary)
allbin.df$n[is.na(allbin.df$n)]<- 0
allbin.df$log_n <- log10(allbin.df$n)


pdf("./figures/Fig1c_variants_distribution.pdf",6.5,5)
ggplot(allbin.df, aes(x =  bin, y =Chrom, fill= log_n)) + 
  geom_tile(alpha =0.9,color = "grey")+theme_cowplot()+
  ylab("Chromosome")+ scale_fill_viridis(option="magma",name = "# of variants", )+
  xlab("Position on the chromosome (5e+06)")
dev.off()
