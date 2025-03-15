#!/usr/bin/env python
import pandas as pd
import re

DIR = "Z:/SHARE/McIntyre_Lab/kinfe/thesis/supplementary"

ujcFile = f"{DIR}/datafile_ujc.csv"
erpFile = f"{DIR}/datafile_ERP.csv"

sexDetFile = f"{DIR}/sex_det.csv"

inUJCDfr = pd.read_csv(ujcFile, low_memory=False)
inERPDfr = pd.read_csv(erpFile, low_memory=False)

inUJCDfr['UJC'] = inUJCDfr['jxnHash']
inERPDfr['ERP'] = inERPDfr['ERP_plus']

print("There are", inUJCDfr['geneID'].nunique(),
      "annotated genes with at least one read in the data.")

# For both datafiles...
loopDct = {
    'UJC': inUJCDfr,
    'ERP': inERPDfr
}

for name, dfr in loopDct.items():

    # ...Remove unneeded columns...
    dfr = dfr.loc[:, ~dfr.columns.str.contains(
        "rep1|rep2|rep3|logUQ|size|ttest|Folded|anno")]
    if name == 'UJC':
        dfr = dfr.loc[:, :'sexClass'].join(dfr[name])

    # ...Change the name of the read count column...
    pattern = r"^.*Cnts_.*_[MF]_rep[456]$"
    for col in dfr.columns:
        if re.match(pattern, col):
            sex = col.split('_')[2]
            rep = col.split('rep')[1]
            dfr = dfr.rename(columns={col: sex+'_rep'+rep+'_readCnt'})

    # ...Count the number of total reads...
    femaleColLst = ['F_rep4_readCnt', 'F_rep5_readCnt', 'F_rep6_readCnt',]
    maleColLst = ['M_rep4_readCnt', 'M_rep5_readCnt', 'M_rep6_readCnt',]
    dfr['F_total'] = dfr[femaleColLst].sum(axis=1)
    dfr['M_total'] = dfr[maleColLst].sum(axis=1)

    # ...Remove UJC/ERP that have 0 reads in reps 4/5/6...
    print("There are", len(dfr), 'data', name, 'in the datafile.')
    dfr = dfr[dfr['F_total'] + dfr['M_total'] != 0]
    print("There are", len(dfr), 'data', name,
          'in the datafile that have at least 1 read in reps 4/5/6.')

    geneDfr = dfr.groupby('geneID').agg(
        num=(name, 'nunique'),
        flag_analyzable=('flag_analyzable', 'sum'),
        unbiased=('sexClass', lambda x: (x == 'unbiased').sum()),
        F_limited=('sexClass', lambda x: (x == 'F_limited').sum()),
        M_limited=('sexClass', lambda x: (x == 'M_limited').sum()),
        F_bias=('sexClass', lambda x: (x == 'F_bias').sum()),
        M_bias=('sexClass', lambda x: (x == 'M_bias').sum()),
        F_rep4_readCnt=('F_rep4_readCnt', 'sum'),
        F_rep5_readCnt=('F_rep5_readCnt', 'sum'),
        F_rep6_readCnt=('F_rep6_readCnt', 'sum'),
        M_rep4_readCnt=('M_rep4_readCnt', 'sum'),
        M_rep5_readCnt=('M_rep5_readCnt', 'sum'),
        M_rep6_readCnt=('M_rep6_readCnt', 'sum'),
    )

    # Only using reps 4/5/6
    geneDfr['flag_F_readCnt0'] = (
        geneDfr[femaleColLst] > 0).any(axis=1).astype(int)
    geneDfr['flag_M_readCnt0'] = (
        geneDfr[maleColLst] > 0).any(axis=1).astype(int)

    geneDfr['F_total_readCnt'] = geneDfr[femaleColLst].sum(axis=1)
    geneDfr['M_total_readCnt'] = geneDfr[maleColLst].sum(axis=1)

    gnExpressRCThresh = str(10)

    geneDfr['flag_gene_expressed_readCnt'+gnExpressRCThresh] = geneDfr.apply(
        lambda row: 1 if row['F_total_readCnt'] +
        row['M_total_readCnt'] >= int(gnExpressRCThresh) else 0,
        axis=1)

    # read counts per gene are equal in UJC/ERP datafile (verified)
    geneDfr = geneDfr.drop(femaleColLst + maleColLst, axis=1)

    colToRename = [col for col in geneDfr.columns if col !=
                   name and 'rep' not in col and 'readCnt' not in col
                   and 'geneID' not in col]

    geneDfr = geneDfr.rename(columns={'num': 'num_'+name})
    geneDfr = geneDfr.rename(
        columns={col: 'num_' + name + '_' + col for col in colToRename})

    geneDfr = geneDfr.reset_index()

    # This number is currently wrong
    print("There are", geneDfr['geneID'].nunique(),
          "dmel6 genes with at least one read in reps 4/5/6.")

    outDfr = geneDfr.sort_values(by='geneID').reset_index(drop=True)

    outDfr.to_csv(
        f'{DIR}/{name}_sex_bias_per_gene.csv', index=False)
