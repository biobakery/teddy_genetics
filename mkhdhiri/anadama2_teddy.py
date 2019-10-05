### To import Anadama2 (for python2)
import sys
sys.path.insert (1, '/n/huttenhower_lab/data/teddy_genetics/mondher_analysis/py3venev/lib_py2')

"""
This script is written to run the Teddy2 analysis as an Anadama2 workflow
"""

import anadama2
from anadama2 import Workflow

### Initializing a workflow instance with Anadama2
#workflow = Workflow(version="0.0.1", description="A workflow to run for the Teddy2 project", remove_options=["input","output"])

workflow = Workflow(version="0.0.1", description="A workflow to run for the Teddy2 project")

#### Steps of the workflow to be run using Anadama2
#___________________________________________________
### Step 1 : preprocessing of the metagenomics data.
#---------------------------------------------------
# Meta_data cleaning

## The inputs are : all_phlan.tsv
## Command to run : ./prep_species.py all_phlan.tsv Day
## 
## The output is : major_species.tsv, major_species_log10.tsv, mp201_sample_mapping.tsv


# 1- Keep major species 
workflow.add_task("./prep_species.py [depends[0]] Day", depends=["all_phlan.tsv"])
# -----<> Outputs : major_species.tsv | major_species_log10.tsv | species.tsv

# 2- Calculating bray Curtis distance
#workflow.add_task ('./prep_bc_trends.py species.tsv init')
## ???? What is the difference between the init and the prev ??
# ----<> Outputs : bc_dists_init.details.tsv | bc_dists_init.tsv


# 3- Compute Species acquisition dates
#workflow.add_task ('./prep_acquisitions.py species.ts')
# ----<> Outputs : acquisitions.tsv

# 4- Compute the diversity indices
#workflow.add_task ('./drop_rows.py -i species.tsv -rm Subject Day -o species_no_meta.tsv')
#workflow.add_task ('./diversity.R -i species_no_meta.tsv -o Diversity_indices.tsv')
# ----<> Outputs : Diversity_indices.tsv (species_no_meta.tsv will be removed)


#___________________________________________________
### Step 2 : preprocessing of the genetics data    .
#---------------------------------------------------
# Cleaning

#workflow.add_task ('./prep_genetics.py snps.ped snps.map species.tsv')
# ----<> Outputs : genetics.tsv
####!!!!!!! only 877 subject are present in the genetics.tsv file while the species.tsv contains a list of 12152 subject   ===> Not all; the subjects have genetics data 
#### !!!!! Why the use of the species at this step ! ! 

### The PC was generated using the Plink 

#workflow.add_task ('plink --file snps --pca header tabs var-wts') 



#___________________________________________________
### Step 3 : Modeling                              .
#---------------------------------------------------
#

# ! The modeling is based on the bray curtis distance rather than the genetics file per se : where is the script to generate the PC from the genetics data or should I generate my own    === > Bc_trends script or bug script (the model could function witn either 

# ! The model is in the format of model = lmer.lmer( "Bug ~ Day + SNP + Day * SNP + (1|Subject)" )  do not understand the argument (1|Subject)
# !  
#!!! Why the model was run in chunks ? the use of results.py script what it will add , to what purpose it was used ???

# workflow.add_task ('./model_rpy2.py species.tsv genetics_pca.tsv results')

#!!! The results  little expl 


#___________________________________________________
### Step 4 : Result representation                 .
#---------------------------------------------------
#


#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



### Work on generic Anadama2 script during the lab meet 

# 1  Focus on the parsing of the data generated through the model in step 3 in order to generate the needed results.

# 2  Repeat the model and rerun the script of the model step by step in order to check what is done
# 3  See the presentation and the results that were produced by Eric and understand them well and understand how to reproduce them

# 4  Integrate the figure generation into the Anadama2 workflow

# 5  Have one script (variant of a workflow) fully functional using Anadam2 for row data to the production of results and figures.

# 6  Test one model (with the PCs at this point) {the use of direct genetics data will be included for the following meeting - need more time to process Eric code} 



# ! How the results were parsed -- Specific parsing script that has already been written -- Or I will write my own ?
# ! what are the scripts used for the results??


#### Other questions : Enrich PCA : script enrich.py : for what purpose was it used 



### Section #5: Run tasks (Required)
workflow.go()









