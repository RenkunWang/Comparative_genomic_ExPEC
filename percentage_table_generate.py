import pandas as pd
import os


#This script generate the ExPEC prevalence table from the ABRicate BLAST result summaries. 

strain_number = "SGL231"
working_folder = r".\alignment\analysis"
strain_folder = os.path.join(working_folder, strain_number)


#Summaries files should be named initially
summary_files = [
    "summary_EB1-36_KO178CA_90_90_11.3.24.tab",
    "summary_EB37-92_KO178CA_90_90_20.3.24.tab",
    "summary_EB93-136_KO178CA_90_90_20.3.24.tab"
]


dfs = []
for filename in summary_files:
    file_path = os.path.join(strain_folder, filename)
    df = pd.read_csv(file_path, sep='\t')
    dfs.append(df)



#process the summary based on the coverage of BLAST > 90% -> 1 else 0
 
merged_df = pd.concat(dfs, axis=0, ignore_index=True)
merged_df.fillna(".", inplace=True)
def process_cell(cell):
    if isinstance(cell, str):
        if ";" in cell:
            max_value = max(float(value) for value in cell.split(";"))
        elif cell == ".":
            return 0
        else:
            try:
                max_value = float(cell)
            except ValueError:
                return 0
    else:
        max_value = cell

    return 1 if max_value > 90 else 0


print(f"Processed table is generating...")


first_two_columns = merged_df.iloc[:, :2]
df_rest = merged_df.iloc[:, 2:]
df_rest = df_rest.apply(lambda col: col.map(process_cell))
df = pd.concat([first_two_columns, df_rest], axis=1)
columns_ordered = list(df.columns[:2]) + sorted(df.columns[2:])
df = df[columns_ordered]
summary_file = os.path.join(strain_folder, f"{strain_number}_processed_summary.csv")
df.to_csv(summary_file, index=False)


print(f"Processed table is done and saved to {summary_file}")
#  merge the summary with the metadata file based on the barcode of genomes
print(f"Merged table is generating...")
metadata_file = os.path.join(working_folder, "EB_PG_v.8_20241010.xlsx")
Metadata = pd.read_excel(metadata_file)
summary = pd.read_csv(summary_file)

Metadata.rename(columns={'Barcode': 'barcode'}, inplace=True)
summary.rename(columns={'#FILE': 'barcode'}, inplace=True)

summary['barcode'] = summary['barcode'].str.replace('.fa$', '', regex=True)

merged_df = pd.merge(Metadata, summary, on='barcode', how='left')

merged_table = os.path.join(strain_folder, f"{strain_number}_merged_table.csv")
merged_df.to_csv(merged_table, index=False)

print(f"Merged table is done and saved to {merged_table}")
# calculate the genes prevalence
print(f"Percentage table is generating...")

merged_df = pd.read_csv(merged_table)
category_column = 'Source.curated'
category_values = ['ExPEC', 'IPEC', 'ND']
result = pd.DataFrame()
for column in merged_df.columns[merged_df.columns.get_loc('NUM_FOUND') + 1:]:
    merged_df[column] = pd.to_numeric(merged_df[column], errors='coerce')

    category_totals = {cat: len(merged_df[merged_df[category_column] == cat]) for cat in category_values}

    for cat in category_values:
        matching_rows = merged_df[(merged_df[column] == 1) & (merged_df[category_column] == cat)]

        percentage = len(matching_rows) / category_totals[cat] if category_totals[cat] > 0 else 0

        result.loc[column, f"{cat}_percentage"] = percentage

    expec_percentage = result.loc[column, 'ExPEC_percentage']
    ipec_percentage = result.loc[column, 'IPEC_percentage']
    nd_percentage = result.loc[column, 'ND_percentage']

    result.loc[column, 'ExPEC_vs_IPEC_difference'] = (expec_percentage - ipec_percentage) if ipec_percentage > 0 else 0
    result.loc[column, 'ExPEC_vs_ND_difference'] = (expec_percentage - nd_percentage) if nd_percentage > 0 else 0


result_percentage = result.copy()
for col in result.columns:
    if pd.api.types.is_numeric_dtype(result[col]):
        result_percentage[col] = (result[col] * 100).round(1)


Percentage_table = os.path.join(strain_folder, f"{strain_number}_percentage_output.csv")
result_percentage.to_csv(Percentage_table)

print(f"Percentage table is done and saved to {Percentage_table}")





