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
workflow.add_task( "./prep_species.py [depends[0]] Day", depends=["all_phlan.tsv"], targets= ['species.tsv'])
#workflow.do("./prep_species.py all_phlan.tsv Day")

# -----<> Outputs : major_species.tsv | major_species_log10.tsv | species.tsv
##

## 2- Calculating bray Curtis distance
workflow.add_task ('./prep_bc_trends.py [depends[0]] init', depends=["species.tsv"])
# ----<> Outputs : bc_dists_init.details.tsv | bc_dists_init.tsv

# 3- Compute Species acquisition dates
workflow.add_task ('./prep_acquisitions.py [depends[0]]', depends = ['species.tsv'])
# ----<> Outputs : acquisitions.tsv

# 4- Compute the diversity indices
workflow.add_task ('./drop_rows.py -i [depends[0]] -rm Subject Day -o [targets[0]]', depends = ['species.tsv'], targets = ['species_no_meta.tsv'])
workflow.add_task ('./diversity.R -i [depends[0]] -o [targets[0]]', depends = ['species_no_meta.tsv'], targets = ['Diversity_indices.tsv'] )
# ----<> Outputs : Diversity_indices.tsv (species_no_meta.tsv will be removed)
#___________________________________________________
### Step 2 : preprocessing of the genetics data    .
#---------------------------------------------------
# Cleaning
workflow.add_task ('./prep_genetics.py [depends[0]] [depends[1]] [depends[2]]', depends = ['snps.ped', 'snps.map', 'species.tsv'], targets = ['genetics.tsv'])
# ----<> Outputs : genetics.tsv
### The PC was generated using the Plink 
workflow.add_task ('plink --file snps --pca header tabs var-wts', depends = ['snps.ped','snps.map'], targets = ['plink.eigenval', 'plink.eigenvec', 'plink.eigenvec.var', 'plink.log', 'plink.nosex']) 
### Clean The raw plink output
workflow.add_task ('./clean_plink.py -i [depends[0]] -o [targets[0]]', depends=['plink.eigenvec'], targets=['genetics_pca.tsv'])


#___________________________________________________
### Step 3 : Modeling                              .
#---------------------------------------------------
#

workflow.add_task ('source /n/huttenhower_lab/data/teddy_genetics/mondher_analysis/py3venev/env/bin/activate')


workflow.add_task (' source /n/huttenhower_lab/data/teddy_genetics/mondher_analysis/py3venev/env/bin/activate | ./model_rpy2.py [depends[0]] [depends[1]] [targets[0]]', depends=['species.tsv', 'genetics_pca.tsv'], targets=['results'])

#workflow.add_task ('deactivate')

#workflow.add_task_gridable ('./model_rpy2.py [depends[0]] [depends[1]] [targets[0]]', depends=['species.tsv', 'genetics_pca.tsv'], targets=['results'], mem = 160000, cores = 1, time = 120)
#___________________________________________________
### Step 4 : Result representation                 .
#---------------------------------------------------
#

### Section #5: Run tasks (Required)
workflow.go()









