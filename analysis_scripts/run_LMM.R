#!/usr/bin/env Rscript

require(docopt)

'Usage:
   run_models.R [-i <microbiome> -g <genetics> -m <metadasta> -c <cores> -o <output>]

Options:
   -i microbiome data
   -g host genetics
   -m metadata
   -c number of cpus [default: 1]
   -o output directory

' -> doc

opts <- docopt(doc)


library(data.table)
library(dplyr)
library(tidyr)
library(survival)
library(traj)
library(doParallel)
library(foreach)
library(lmerTest)
library(performance)

###########################
##Get the smallest value greater than 0.
microbiome_melt <- fread(opts$i, sep="\t", header=T) %>%
  mutate(sample_id = factor(sample_id)) %>%
  melt() %>%
  dplyr::filter(value > 0)
eps <- min(microbiome_melt$value)

##Add the smallest value greater than 0 to the relative abundance.
microbiome <- fread(opts$i, sep="\t", header=T) %>%
  mutate(sample_id = as.character(sample_id)) %>%
  melt() %>%
  mutate(value = log(value + eps)) %>%
  dcast(sample_id ~ variable, value.var="value")

# Load genetic data
genetics <- fread(opts$g, sep="\t", header=T, check.names=F)

## Load metadata
metadata <- fread(opts$m, header=T)

###combine those three data
df_sample<- merge(metadata, genetics, by="subject_id")
microbiome$sample_id <- as.numeric(microbiome$sample_id )
dat <- merge(df_sample, microbiome, by="sample_id") %>%
  rename(Subject=subject_id, Day=Days) %>%
  mutate(Subject=factor(Subject))
###########################
# make a list of the models

formulae <- vector()
for (i in colnames(microbiome)[-1]) {
  
  for (j in colnames(genetics)[-c(1)]) {
    
    formula <- data.frame(formula = paste0("lmer(", i, " ~ ", j, " * ", "Day + Clinical_Center + ( 1 + Day | Subject ), data=dat)"))
    
    formulae <- rbind(formulae, formula)	
  }
}


formulae <- levels(as.factor(formulae$formula))

###################################
# loop through the list in parallel

library(doParallel)

registerDoParallel(cores=opts$c)


system.time(model_anova_total <- foreach(i = formulae, 
                                         .combine = rbind,  # Combine results by row
                                         .packages = c("lmerTest", "tidyverse", "performance"),
                                         .errorhandling = 'remove'  # Skip failed iterations
) %dopar% {
  
  # Fit the mixed-effects model using the given formula
  model <- eval(parse(text = i))
  
  # Extract fixed effect coefficients
  model_coef <- as.data.frame(coef(model)$Subject) %>%
    select(contains("PC"))
  
  # Perform ANOVA on the model
  model_anova <- data.frame(anova(model)) %>%
    tibble::rownames_to_column("Predictor") %>%
    select(Predictor, `Pr..F.`) %>%
    rename(P = `Pr..F.`) %>%
    mutate(Type = gsub("PC.", "PC", Predictor)) %>%
    mutate(Predictor = factor(Predictor)) %>%
    mutate(Model = i) %>%
    mutate(Coefficient = c(.subset2(model_coef, 1)[1], NA, NA, .subset2(model_coef, 2)[1]))
  
  return(model_anova)
}
)


model_anova_total <- model_anova_total %>% filter(grepl("PC",Predictor))
fwrite(model_anova_total, file = file.path(opts$o, "model_anova_total.tsv"), sep="\t")
