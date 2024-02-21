
import pandas as pd


test = pd.read_csv('/mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output6/Amplicon_design_output/SNP_inclusion-31-400.csv')

print(test['Amplicon_ID'].value_counts())