
##### R - script to test diab

library ('base')
library ('stats')
library ("vctrs", lib.loc = "/n/huttenhower_lab/data/teddy_genetics/mondher_analysis/R_lib")
library ("lmerTest", lib.loc = "/n/huttenhower_lab/data/teddy_genetics/mondher_analysis/R_lib")

## Data Preparation 
## diab subjects

diab = read.csv ('cleaned_mp201_t1d3_case_control.tsv', sep = '\t')
## Need two column from the diab: the subject ID and wether or not they developed diabetes
subject = diab$Mask.Id
diabetes = diab$Outcome


## Need the subject - sample table to calculate mean shannon diversity indices
sub_sample = read.csv('metadata_sample-subject-day.tsv', sep = '\t')



## Need microbiome diversity column
diver = read.csv ('versity_indices.tsv', sep = '\t')



## The model to use

model = lmer ("diab ~ microbiome + genes + metadata + (1|Subject)")






