import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
# this script read the information of COG categories from the excel table and generate the cog chart to present the proportions of different categories. COG classes from NCBI

df = pd.read_excel(r".\alignment\analysis\Gene_selection\ecoli_cogclassifier\COG.xlsx", engine='openpyxl')

# cog for prediction from genomic data can be possibly different categories
all_cogs = []
for cat in df['COG_category'].dropna():
    all_cogs.extend(list(str(cat).strip())) 


cog_counts = Counter(all_cogs)
total = sum(cog_counts.values()) 


cog_map = {
    'E': 'Amino acid transport and metabolism',
    'G': 'Carbohydrate transport and metabolism',
    'C': 'Energy production and conversion',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'K': 'Transcription',
    'J': 'Translation, ribosomal structure and biogenesis',
    'P': 'Inorganic ion transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'R': 'General function prediction only',
    'T': 'Signal transduction mechanisms',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'L': 'Replication, recombination and repair',
    'I': 'Lipid transport and metabolism',
    'S': 'Function unknown',
    'F': 'Nucleotide transport and metabolism',
    'N': 'Cell motility',
    'V': 'Defense mechanisms',
    'X': 'Mobilome: prophages, transposons',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'W': 'Extracellular structures',
    'B': 'Chromatin structure and dynamics',
    'A': 'RNA processing and modification',
    'Z': 'Cytoskeleton'
}


percentages = {k: (v / total) * 100 for k, v in cog_counts.items()}
sorted_keys = sorted(percentages.keys())  
labels = [f"{k} ({cog_map.get(k, 'Unknown')})" for k in sorted_keys]
values = [percentages[k] for k in sorted_keys]


plt.figure(figsize=(12, 8))
plt.barh(labels, values, color='lightseagreen')
plt.xlabel("Percentage of Genes (%)")
plt.title("60 TOP shared COG Accessory")
plt.tight_layout()
plt.show()
