"""
This script is written to run the Teddy analysis as an Anadama2 workflow
"""

import anadama2
from anadama2 import Workflow

### Initializing a workflow instance with Anadama2
workflow = Workflow(version="0.0.2", description="A workflow to run for the Teddy2 project", remove_options=["input","output"])


#### Steps of the workflow to be run using Anadama2

#___________________________________________________
### Step 1 : preprocessing of the metagenomics data.
#---------------------------------------------------
# Meta_data cleaning

## The inputs are : all_phlan.tsv
## Command to run : ./prep_species.py all_phlan.tsv Day
## 
## The output is : major_species.tsv, major_species_log10.tsv, mp201_sample_mapping.tsv




# Add custom arguments and parse arguments (Optional)
#workflow.add_argument("phlan_file", desc="The MetaPhlan2 file", default="tsv")
#args = workflow.parse_args()

#workflow.add_argument("meta", desc="The MetaData colummn name embedded within the MetaPhlan2 file")
#args = workflow.parse_args()



### Section #3: Get input/output file names (Optional)
#in_files = workflow.get_input_files(extension=args.input_extension)
#out_files = workflow.name_output_files(name=in_files, tag="metaphlan2_preprocessing")


### Section #4: Add tasks (Required)


workflow.add_task("prep_species.py [depends[0] depends[1]]", depends=["all_phlan.tsv", "Day"], targets=[])



#workflow.add_task(
 #   "metaphlan2.py [depends[0]] --input_type fastq --no_map > [targets[0]]",
  #  depends=["sample1.fastq"], targets=["sample1_profile.txt"])


### Section #5: Run tasks (Required)
workflow.go()
















