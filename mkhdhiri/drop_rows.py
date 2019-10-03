#!/net/rcnfs03/srv/export/huttenhower_lab/share_root/data/teddy_genetics/mondher_analysis/py3venev/anaconda3/bin/python
"""
 This script is written to clean/remove row data from input file : species.tsv for Anadama2 Teddy2 project
"""

############# Loading Libraries #####
print ('\n Loading libraries... \n')

import argparse
import pandas as pd
print ('Libraries has been loaded \n')
###################################

parser = argparse.ArgumentParser(description = 'Parse Arguments')

parser.add_argument ('-i', '--input', help = 'The input file', required = True, dest = 'i')
parser.add_argument ('-o', '--output', help = 'The output file', required = True, dest = 'o')
parser.add_argument ('-rm', '--remove', nargs = '+', help = 'The rows to be removed', required = True, dest = 'r')
args = parser.parse_args()

### Read the input data file
df = pd.read_csv (args.i, sep = '\t')
#############################

### Create a list of index for the desired to be removed rows 

index_list = []
for r in args.r:
	ind = df[df['ID']==r].index.values.astype(int)[0]
	index_list.append (ind)

for ind in index_list:
	df = df.drop(ind , inplace=False)
df = df.reset_index(drop = True)

df.to_csv (args.o ,index = False , sep = '\t')	










