import pandas as pd
count_file = pd.read_csv(r"GSE130325_Cirovic_RNA_counts.txt",delim_whitespace=True)
count_file.to_csv (r'GSE130325_Cirovic_RNA_counts.csv', index=None)