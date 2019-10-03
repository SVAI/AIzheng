"""
Extract an expression matrix and label file for a set of tumor types from the
PanTCGA dataset. The metadata spreadsheet should contain the following sheets:

- SampleInfo
- Clinical
- Exposure
- File-Case-Mapping
"""
import argparse
import numpy as np
import pandas as pd



def select_rows_by_values(df, column, values):
	return pd.DataFrame().append([df[df[column] == v] for v in values], sort=False)



if __name__ == "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--data", help="input expression matrix", default="PanTCGA_Expression_Data/Kidney_FPKM_Quantile_No-Outliers_5_17_18.tab")
	parser.add_argument("--metadata", help="metadata spreadsheet", default="PanTCGA_Expression_Data/PanTCGA_Worksheet-v1-20180514.xlsx")
	parser.add_argument("--tumor-types", help="list of tumor types", nargs="*", default=["TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP"])
	parser.add_argument("--output-data", help="output expression matrix", default="kidney.emx.txt")
	parser.add_argument("--output-labels", help="output labels", default="kidney.labels.txt")

	args = parser.parse_args()

	print("loading input files...")

	# load expression data
	data = pd.read_csv(args.data, sep="\t")

	# load sample info spreadsheet
	sample_info = pd.read_excel(args.metadata, sheet_name="SampleInfo", index_col=0)
	sample_info = sample_info.iloc[0:11093]
	sample_info = select_rows_by_values(sample_info, "Project ID", args.tumor_types)

	# append columns from clinical data to sample info
	clinical = pd.read_excel(args.metadata, sheet_name="Clinical", index_col=0)
	file_case_mapping = pd.read_excel(args.metadata, sheet_name="File-Case-Mapping", index_col=0)

	sample_info["case_id"] = file_case_mapping.loc[sample_info["File Name"], "cases.0.case_id"].values
	sample_info["tumor_stage"] = clinical.loc[sample_info["case_id"], "tumor_stage"].values

	print("extracting data and labels...")

	# extract common samples
	sample_names = set(sample_info.index).intersection(set(data.columns))

	# extract samples that have labels
	data = data[sample_names]

	# create labels dataframe
	labels = pd.DataFrame({
		"tumor_type": sample_info["Project ID"],
		"tumor_state": (sample_info["Sample Type"] == "Solid Tissue Normal").astype(int),
		"tumor_stage": sample_info["tumor_stage"]
	}, index=sample_info.index)

	# map tumor stages to numeric values
	mapping = {
		"stage i": 1,
		"stage ii": 2,
		"stage iii": 3,
		"stage iv": 4,
		"not reported": np.nan
	}

	for key in mapping:
		labels.loc[labels["tumor_stage"] == key, "tumor_stage"] = mapping[key]

	# save output dataframes
	print("saving output files...")

	data.to_csv(args.output_data, sep="\t", na_rep="NA", float_format="%.8f")
	labels.to_csv(args.output_labels, sep="\t", na_rep="NA")
