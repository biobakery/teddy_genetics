### To import Anadama2 (for python2)
import sys
sys.path.insert (1, '/n/huttenhower_lab/data/teddy_genetics/mondher_analysis/py3venev/lib_py2')
# pip install --target= /n/huttenhower_lab/data/teddy_genetics/mondher_analysis/py3venev/lib_py2 rpy2
"""
This script is written to run the Teddy2 analysis as an Anadama2 workflow
"""
import anadama2
from anadama2 import Workflow

### Initializing a workflow instance with Anadama2
workflow = Workflow(version="0.0.1", description="A workflow to run for the Teddy2 project", remove_options=["input","output"])

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
workflow.do("./prep_species.py all_phlan.tsv Day")
# -----<> Outputs : major_species.tsv | major_species_log10.tsv | species.tsv
# 2- Calculating bray Curtis distance
workflow.do ('./prep_bc_trends.py species.tsv init')
# ----<> Outputs : bc_dists_init.details.tsv | bc_dists_init.tsv

# 3- Compute Species acquisition dates
workflow.do ('./prep_acquisitions.py species.tsv')
# ----<> Outputs : acquisitions.tsv

# 4- Compute the diversity indices
workflow.do ('./drop_rows.py -i species.tsv -rm Subject Day -o species_no_meta.tsv')
workflow.do ('./diversity.R -i species_no_meta.tsv -o Diversity_indices.tsv')
# ----<> Outputs : Diversity_indices.tsv (species_no_meta.tsv will be removed)
#___________________________________________________
### Step 2 : preprocessing of the genetics data    .
#---------------------------------------------------
# Cleaning
workflow.do ('./prep_genetics.py snps.ped snps.map species.tsv')
# ----<> Outputs : genetics.tsv
### The PC was generated using the Plink 
workflow.do ('plink --file snps --pca header tabs var-wts') 
### Clean The raw plink output
workflow.do ('./clean_plink.py -i plink.eigenvec -o genetics_pca.tsv')

#___________________________________________________
### Step 3 : Modeling                              .
#---------------------------------------------------
#
workflow.do ('./model_rpy2.py species.tsv genetics_pca.tsv results')
#___________________________________________________
### Step 4 : Result representation                 .
#---------------------------------------------------
#


#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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









