#!/n/helmod/apps/centos7/Core/R_core/3.5.1-fasrc01/bin/Rscript


# This script was written by Mondher Khdhiri <mkhdhiri@broadinstitute.org>
## Date : 27 Sep 2019

# This script is written to run preprocessed data from the MetaPhlan2 in order
# to generate different diversity indices. 

## The typical input for this script is a MetaPhlan2 file
## The typical output is a file with the diversity indices within

############# Loading Libraries #####
print ('Loading libraries...')
library ('vegan')
library ('argparse', lib.loc = '/net/rcnfs03/srv/export/huttenhower_lab/share_root/data/teddy_genetics/mondher_analysis/R_lib')
print ('Libraries has been loaded')
######################################

############## Parse required Arguments ####
## Potential modification is to add a parameter for the choice of the 
## diversity parameter to be calculated.

parser = ArgumentParser(description = 'Parse Arguments')
parser$add_argument('-in', metavar='i', type='character', help='The input file', required= TRUE)
parser$add_argument('-out', metavar='o', type='character', help='The output file')#, required = TRUE)

args = parser$parse_args()

############################################
#### Read the input data

df = read.csv(args$i, sep = '\t')

##### Processing data 

datalist = list()
j = length (df)
for (i in 2:j) {

    x = diversity (df[i], index = 'shannon')
    y = diversity (df[i], index = 'simpson')

    dat = data.frame(Shannon = x, Simpson = y)
    dat$i = i  #  keep track of which iteration number
    datalist[[i]] = dat # add it to the list
 
}

# Generating the dataframe with the needed index
big_data = do.call(rbind, datalist)

# Changing the index list into sample names 

for (i in 1:j-1){
k = i+1
big_data$i[i] = colnames (df[k])
}

### Changing col names and order
colnames(big_data)[colnames(big_data)=='i'] =  'Subject_ID'
big_data = big_data[c('Subject_ID', 'Shannon', 'Simpson')]

## Print into the output file 

print ('Genetratin the output file... ')

write.table(big_data, file = args$o, row.names = FALSE, dec = '.', sep = '\t', quote = FALSE)

print (' ______________________DONE!____________________________')

#######################################################################################################
