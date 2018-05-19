#
# separate-tumor-types.py
# script to separate tumor types based on the stage of cancer from the
# panTCGA RNA-seq tumor dataset
#

import pandas as pd 
import numpy as np 

NUM_TUMORS = 11092



if __name__ == '__main__':

	print('reading data files...')

	# retrieve labels worksheet containing labels
	sample_info = pd.read_excel('./PanTCGA_Expression_Data/PanTCGA_Worksheet-v1-20180514.xlsx', sheetname='SampleInfo')

	# retrieve worksheet containing the specific types of cancer
	clinical = pd.read_excel('./PanTCGA_Expression_Data/PanTCGA_Worksheet-v1-20180514.xlsx', sheetname='Clinical')

	# retrieve panTCGA data
	pan_kidney_data = pd.read_csv('./PanTCGA_Expression_Data/Kidney_FPKM_Quantile_No-Outliers_5_17_18.tab', sep='\t')

	print('creating datasets...')
	# create dictionary for indexes of kirc, kirp, and kich
	idxs = {}
	idxs['kirc']  = sample_info.index[sample_info['Project ID'] == 'TCGA-KIRC'].tolist()
	idxs['kirp']  = sample_info.index[sample_info['Project ID'] == 'TCGA-KIRP'].tolist()
	idxs['kich']  = sample_info.index[sample_info['Project ID'] == 'TCGA-KICH'].tolist()

	# get rid of label indices that are htseqcount labels, only want FPKM
	for k in idxs:
		idxs[k] = [x for x in idxs[k] if x < NUM_TUMORS]

	# partition healthy vs cancerous
	healthy_list = list(sample_info.index[sample_info['Sample Type'] == 'Solid Tissue Normal'])
	healthy_list = [x for x in healthy_list if x < NUM_TUMORS]

	healthy = {}
	for k in idxs:
		healthy[k] = []
		for h in healthy_list:
			if h in idxs[k]:
				healthy[k].append(h)

	unhealthy = {}
	for k in idxs:
		unhealthy[k] = []
		for h in idxs[k]:
			if h not in healthy_list:
				unhealthy[k].append(h)

	healthy_labels = {}
	for k in idxs:
		healthy_labels[k] = []
		for i in healthy[k]:
			healthy_labels[k].append(str(sample_info['GEM_Header'][i]))

	unhealthy_labels = {}
	for k in idxs:
		unhealthy_labels[k] = []
		for i in unhealthy[k]:
			unhealthy_labels[k].append(str(sample_info['GEM_Header'][i]))


	healthy_df = pd.DataFrame()
	for k in healthy_labels:
		for l in healthy_labels[k]:
			healthy_df[l] = pan_kidney_data[l]

	unhealthy_df = pd.DataFrame()
	for k in unhealthy_labels:
		for l in unhealthy_labels[k]:
			if l in pan_kidney_data.columns:
				unhealthy_df[l] = pan_kidney_data[l]
			else:
				print (l + ' not found')

	print('writing files...')
	healthy_df.to_csv('./healthy_kidney_data.csv', sep='\t')
	unhealthy_df.to_csv('./unhealthy_kidney_data.csv', sep='\t')

	print('done')
















