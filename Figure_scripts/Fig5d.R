############################################################
## Fig5d:
## Toddler: SNP × time interaction plots 
############################################################

## Set project root
setwd("~/ddong/TEDDY_project/")

###########################
## 0. Load packages
###########################
library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(cowplot)

###########################
## 1. Load microbiome data (toddler) and log10 transform
###########################
microbiome_melt <- fread("./output/microbiome_toddler.txt", sep = "\t", header = TRUE) %>%
  mutate(sample_id = factor(sample_id)) %>%
  melt() %>%
  dplyr::filter(value > 0)

eps <- min(microbiome_melt$value)

microbiome <- fread("./output/microbiome_toddler.txt", sep = "\t", header = TRUE) %>%
  mutate(sample_id = factor(sample_id)) %>%
  melt() %>%
  mutate(value = log10(value + eps)) %>%
  dcast(sample_id ~ variable, value.var = "value")

###########################
## 2. Load genetics (PCs + SNP genotypes) and metadata
###########################
PCs <- fread("./data/genetics.bed/filter/teddy.qc.PCA_PC1-PC5.tsv",
             sep = "\t", header = TRUE, check.names = FALSE)
colnames(PCs)[1] <- "subject_id"

snps <- fread("./data/genetics.bed/filter/sample_genotype.txt")
subject <- fread("./data/genetics.bed/filter/subject.txt", header = FALSE)
snps$subject_id <- subject
snps <- snps %>% select(subject_id, colnames(snps))

metadata <- fread("./data/metadata/metadata_sample.txt", sep = "\t", header = TRUE)

###########################
## 3. Merge into analysis data frame
###########################
df_subjects <- merge(metadata, PCs, by = "subject_id") %>%
  inner_join(snps, by = "subject_id")

microbiome$sample_id <- as.numeric(as.character(microbiome$sample_id))

dat <- df_subjects %>%
  inner_join(microbiome, by = "sample_id") %>%
  rename(Subject = subject_id, Day = Days) %>%
  mutate(Subject = factor(Subject))

###########################
## 4. Load significant SNP–species pairs and select top 5 SNPs as example
###########################
sig <- fread("./output/toddler_species_interaction_sig.txt")
sig <- sig %>% arrange(P)
sig <- sig[1:5,]

dat2 <- dat %>%
  select(Subject, sample_id, Day, one_of(sig$Species), one_of(sig$snp))

pairs <- sig %>%
  select(snp, Species) %>%
  arrange(Species)

###########################
## 5. Plot: species ~ Day by genotype 
###########################
pdf("./figures/main/Fig5d_toddler_snp_interaction.pdf",
    width = 5.5, height = 4)

for (i in 1:nrow(pairs)) {
  
  species <- as.character(pairs[i, 2])
  s <- as.character(pairs[i, 1])
  
  dat3 <- dat2 %>%
    select(Day, one_of(species), one_of(s)) %>%
    filter(!is.na(get(s))) %>%
    filter(get(species) != (-6.69897))
  
  print(
    ggplot(dat3, aes(y = get(species), x = Day, color = as.factor(get(s)), fill = "grey")) +
      scale_color_brewer(palette = "Dark2") +
      scale_fill_manual(values = "gray60") +
      theme_cowplot() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.8)
      ) +
      labs(x = "Time (Days)", color = s) +
      geom_smooth(method = "lm") +
      ylab(species) +
      xlab(s) +
      guides(fill = "none")
  )
}

dev.off()
