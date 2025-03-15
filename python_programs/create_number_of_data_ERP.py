#!/usr/bin/env python

import pandas as pd

# file = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/rmg_erp_output/fiveSpecies_2_dmel6_ujc_er_vs_dmel_data_2_dmel6_ujc_noMultiGene_read_per_ERP.csv"
file = "B:/sex_specific_splicing/rmg_erp_output/fiveSpecies_2_dmel6_ujc_er_vs_dmel_data_2_dmel6_ujc_noMultiGene_read_per_ERP.csv"

dfr = pd.read_csv(file, low_memory=False)

grpDfr = dfr.groupby(['sampleID', 'geneID'])[
    'ERP_plus'].nunique().reset_index()
grpDfr = grpDfr.rename({'ERP_plus': 'numERP'}, axis=1)

perSampleDfr = grpDfr.groupby('sampleID')['numERP'].sum().reset_index()

perGeneDfr = dfr[['geneID', 'ERP_plus']]
perGeneDfr = perGeneDfr.drop_duplicates()

print(f"Number of unique ERP across samples: {len(perGeneDfr)}")

DIR="Z:/SHARE/McIntyre_Lab/kinfe/thesis/supplementary"

perSampleDfr.to_csv(
    f"{DIR}/number_of_data_ERP.csv", index=False)
