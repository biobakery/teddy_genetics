#!/usr/bin/env Rscript
require(docopt)
'Usage:
   run_trajectory.R [ -g <genetics> -m <metadata> -o <output>]

Options:
   -g host genetics
   -m metadata (Including microbiome Bray-Curtis distance, antibody, and T1D information)
   -o output directory
' -> doc

opts <- docopt(doc)


library(data.table)
library(dplyr)
library(tidyr)
library(survival)
library(traj)
library(doParallel)

### Data Loading
metadata <- fread(opts$m)
geneticpc <- fread(opts$g) %>% rename(subject_id = 1)


### Data Filtering
allphenotype <- metadata %>%
  inner_join(geneticpc, by="subject_id") %>%
  filter(Days <= 800, (is.na(T1D) | T1D == FALSE) & (is.na(AB) | AB == FALSE)) %>% 
 # left_join(microbiome, by=c("sample_id","subject_id")) %>%
  filter(!is.na(bc_dists_init)) %>%
  group_by(subject_id) %>% 
  filter(n() > 3) %>%
  ungroup()

### Transform data to a matrix with columns as Day and rows as individuals
allphenotype.df <- allphenotype %>% dcast(subject_id~Days,value.var = "bc_dists_init")
allphenotype.df <- allphenotype.df %>% rename(ID = subject_id)

### Trajectory Clustering
step1 <- Step1Measures(allphenotype.df, Time = as.numeric(colnames(allphenotype.df)[-1]), ID = TRUE)
step2 <- Step2Selection(step1)
step3 <- Step3Clusters(step2, nclusters = 3)

cluster.df <- step3$partition %>% rename(subject_id = ID)

### Merge Data and Scale Genetic PCs
data_AB_T1D <- metadata %>% 
  inner_join(geneticpc, by="subject_id") %>%
  mutate(across(starts_with("PC"), ~ . / sd(.))) %>%
  left_join(cluster.df, by="subject_id") %>%
  mutate(Cluster = as.character(Cluster))

### Model for calculating relative risk
model <- coxph(Surv(AB_T1D_time, AB_T1D) ~ 
                 Clinical_Center +  Cluster + PC1 + PC2 + PC3 + PC4 + PC5 +
                 PC1:Cluster + PC2:Cluster + PC3:Cluster + PC4:Cluster + PC5:Cluster +
                 Sex + FDR + delivery_class +
                 time_to_brstfed_stop + age_startsolid + Antibiotics + Probiotic,
               data = data_AB_T1D)

model_summary <- summary(model)
model_output <- cbind(model_summary$coefficients, model_summary$conf.int[, 3:4])

### Save Results
fwrite(cluster.df, file = file.path(opts$o, "cluster_results.tsv"), sep = "\t")
fwrite(model_output, file = file.path(opts$o, "relative_risk.tsv"), sep = "\t")
