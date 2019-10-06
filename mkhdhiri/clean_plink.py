#!/net/rcnfs03/srv/export/huttenhower_lab/share_root/data/teddy_genetics/mondher_analysis/py3venev/anaconda3/bin/python
"""
This script is written to parse the raw output of PLINK
"""

__author__ = 'Mondher Khdhiri <mkhdhiri@hsph.harvard.edu>'
__date__ = ''

############# Loading Libraries #####
print ('\n Loading libraries... \n')

import argparse
import pandas as pd
print ('Libraries has been loaded \n')
###################################

parser = argparse.ArgumentParser(description = 'Parse Arguments')

parser.add_argument ('-i', '--input', help = 'The input file', required = True, dest = 'i')
parser.add_argument ('-o', '--output', help = 'The output file', required = True, dest = 'o')
args = parser.parse_args()

### Read the input data file
df = pd.read_csv (args.i, sep = '\t')
#############################

## Remove the FID column
df = df.drop(columns =['FID'])
## Transpose the data without index
df = df.set_index('IID').T
df.index.names = ['IID']
df.to_csv (args.o ,index = True , sep = '\t')

##################################

