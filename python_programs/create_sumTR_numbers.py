# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 01:57:25 2025

@author: knife
"""

import pandas as pd

DIR = "Z:/SHARE/McIntyre_Lab/kinfe/thesis/supplementary"

readAndUJCFile = f"{DIR}/number_of_data_ujc_and_read.csv"
erpFile = f"{DIR}/number_of_data_ERP.csv"

readAndUJCDfr = pd.read_csv(readAndUJCFile, low_memory=False)
readAndUJCDfr.drop(readAndUJCDfr.tail(3).index,
                   inplace=True)  # drop last n rows
erpDfr = pd.read_csv(erpFile, low_memory=False)
erpDfr.drop(erpDfr.tail(3).index, inplace=True)  # drop last n rows

mergeDfr = pd.merge(readAndUJCDfr, erpDfr, how='outer', on='sampleID')


mergeDfr[['species', 'sex', 'rep', 'TR']
         ] = mergeDfr['sampleID'].str.split('_', expand=True)
mergeDfr[['empty', 'repNum']] = mergeDfr['rep'].str.split('rep', expand=True)
mergeDfr['repNum'] = mergeDfr['repNum'].astype(int)

rep456 = mergeDfr[mergeDfr['repNum'] > 3].copy()

rep456['sample'] = rep456['species'] + \
    '_' + rep456['sex'] + '_' + rep456['rep']

outDfr = rep456[['sample', 'numRead', 'numUJC', 'numERP']
                ].groupby('sample').agg(sum)

outDfr.loc['Total:'] = outDfr.sum()

outDfr.to_csv(f'{DIR}/number_of_data_sumTR_read_ujc_ERP.csv')
