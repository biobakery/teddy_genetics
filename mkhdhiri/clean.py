#!/usr/bin/python3

################################################
__author__= "mkhdhirii@broadinstitute.org"
__date__ = "21/10/2019"
__version__ = '0.0.1'

################################################
import numpy as np
import pandas as pd
import statistics
import plotly.express as px
from functools import reduce

################################################
def read_data (data_file):
	x = pd.read_csv (data_file, sep = '\t')
	return (x)
################################################
# read the needed data
diab = read_data ('cleaned_mp201_t1d3_case_control.tsv')
subject_sample = read_data ('metadata_sample-subject-day.tsv')
diversity = read_data ('Diversity_indices.tsv')
genetics = read_data ('genetics_pca.tsv')
metadata = read_data ('clean_mp201_child_medications.tsv')

# Cleaning the diversity matrix
diversity['Sample_ID'] = diversity['Sample_ID'].map(lambda x: x.lstrip('X'))
# Convert first column into int64
diversity = diversity.astype ({'Sample_ID': 'int64'})
# Transpose the genetics matrix and keep the first row as header
new_gen = genetics.T.reset_index(drop = False)
new_header = new_gen.iloc[0] #grab the first row for the header
gen = new_gen[1:] #take the data less the header row
gen.columns = new_header #set the header row as the df header

# Change NAN in the medication to None value
metadata = metadata.replace (np.nan,'None', regex = True)

## Calculate mean diversity and per subject
div_sub = diversity.merge (subject_sample, how = 'left', left_on = diversity['Sample_ID'], right_on = subject_sample['Sample_ID'])

dic = {}
sub_shan = {}
for subject in div_sub ['Subject']:
	subset =  div_sub.loc[div_sub['Subject']== subject]
	dic[subject] = {}
	dic[subject]['Day'] = subset['Day'].tolist()
	dic[subject]['Sample_ID'] = subset['key_0'].tolist()
 	dic[subject]['Shannon'] = subset['Shannon'].tolist()
	dic[subject]['Simpson'] = subset['Simpson'].tolist()
	dic[subject]['Shannon_mean'] = subset['Shannon'].mean()
	dic[subject]['Shannon_std'] = subset['Shannon'].std()
	sub_shan[subject] = subset['Shannon'].tolist()

## plot the Shannon indices for all the subjects 
df = div_sub[['Subject', 'Shannon']]
fig = px.box (df, x = 'Subject', y = 'Shannon')
# write the image into png file
fig.write_image('Fig1.png')

#######################################################
#
## Selecting data to be merged
#
gen = [gen['IID'], gen['PC1']]
 




db = diab[['Mask_Id', 'Outcome']]
med = metadata[['Mask_Id', 'Medname']]
diab_med = db.merge (med, how = 'left', left_on = 'Mask_Id', right_on = 'Mask_Id')

diab_med_gen = diab_med.merge (gen, how= 'left', left_on = 'Mask_Id', right_index = True)



gen = genetics[['IID','PC1']]


daf.set_index('Attribute',inplace=True)


data_frame = data_frame.T
### Merge all necessary data for the model into one dataframe
#
#- Data that I need to merge: SubjectID->diab->Diversity->genetics_PCA->metadata

df_all = [data_frame, diab, genetics, metadata]










