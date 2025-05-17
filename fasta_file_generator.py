import pandas as pd
import os

#define the strain number, based on the information in the record table, generate corresponding fasta file
strain_number = "SGL623"

alignment_folder = r".\alignment"

excel_file = os.path.join(alignment_folder, f"{strain_number}.xlsx")
txt_file = os.path.join(alignment_folder, f"{strain_number}_sequence.txt")
output_fasta = os.path.join(alignment_folder, f"{strain_number}_elements.fasta")


df = pd.read_excel(excel_file)
#the whole genome read from the sequence txt and connect the different contigs
whole_genome = ""

df['Annotation'] = df['Annotation'].str.replace(' ', '_', regex=False)

with open(txt_file, "r", encoding="utf-8", errors="ignore") as file:
    for line in file:
        
        stripped_line = line.strip()
        
        if stripped_line and stripped_line[0].isdigit():
            
            cleaned_line = ''.join(stripped_line.split()[1:])  
            whole_genome += cleaned_line

# generate the nucleotide sequence based on the genome position and genome sequence
with open(output_fasta, "w") as fasta_output:
    for index, row in df.iterrows():
        position_range = row['Position'].replace('..', '.').split('.')
        start = int(position_range[0]) - 1  
        stop = int(position_range[1])
        element_sequence = whole_genome[start:stop]
        label = f">{row['ROD']}_{row['Strain']}_{row['Annotation']}_{row['Position']}"
        fasta_output.write(f"{label}\n{element_sequence}\n")

print(f"FASTA file has been created: {output_fasta}")
