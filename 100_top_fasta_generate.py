import pandas as pd

# This script read the percentage table and generate the top genes list and corresponding sequence fasta file. The sequence in fasta file are searched in the combined fasta for the selected genes. 
csv_path = r".\alignment\analysis\Gene_selection\csv\combined.csv"
txt_fasta_path = r".\alignment\analysis\Gene_selection\csv\combined.txt"
top_genes_txt_path = r".\alignment\analysis\Gene_selection\csv\100_top_genes.txt"
top_genes_fasta_path = r".\alignment\analysis\Gene_selection\csv\100_top_genes.fasta"


df = pd.read_csv(csv_path)
df_sorted = df.sort_values(by='ExPEC_vs_ND_difference', ascending=False)


top_100_genes = df_sorted.iloc[:100, 0].astype(str).tolist()


with open(top_genes_txt_path, 'w') as f:
    for gene in top_100_genes:
        f.write(gene + '\n')


with open(txt_fasta_path, 'r', encoding='utf-8', errors='ignore') as f:
    lines = f.readlines()


matched_sequences = []
current_header = ''
current_seq_lines = []

def save_sequence(header, seq_lines):
    if header in top_100_genes:
        matched_sequences.append('>' + header + '\n' + ''.join(seq_lines))

for line in lines:
    if line.startswith('>'):
        if current_header:
            save_sequence(current_header, current_seq_lines)
        current_header = line[1:].strip()
        current_seq_lines = []
    else:
        current_seq_lines.append(line)
if current_header:
    save_sequence(current_header, current_seq_lines)

with open(top_genes_fasta_path, 'w') as f:
    f.writelines(matched_sequences)
