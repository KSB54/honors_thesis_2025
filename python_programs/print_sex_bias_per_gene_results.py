# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 19:58:41 2025

@author: knife
"""

import pandas as pd
import matplotlib.pyplot as plt
import re

# Printing per gene results

DIR = "Z:/SHARE/McIntyre_Lab/kinfe/thesis/"

# grab files...

ujcFile = f"{DIR}/supplementary/UJC_sex_bias_per_gene.csv"
erpFile = f"{DIR}/supplementary/ERP_sex_bias_per_gene.csv"

ujcDfr = pd.read_csv(ujcFile, low_memory=False)
erpDfr = pd.read_csv(erpFile, low_memory=False)

# Number of unique UJC/ERP in the data
ujcDfr['num_UJC'].sum()
erpDfr['num_ERP'].sum()

# GENE EXPRESSION IN LONG READ DATA

# Print number of genes with at least one analyzable and/or sex-limited jxnHash
numAnlyzblOrSxLmtGn = len(ujcDfr[(ujcDfr['num_UJC_flag_analyzable'] >= 1) | (
    ujcDfr['num_UJC_F_limited'] >= 1) | (ujcDfr['num_UJC_M_limited'] >= 1)])

print("There are", numAnlyzblOrSxLmtGn, "genes with at least one "
      "analyzable and/or sex-limited jxnHash.")
print()

# Print number of expressed genes
numExpressed = ujcDfr['flag_gene_expressed_readCnt10'].sum()
print("There are", numExpressed, " expressed genes.")
print()

sexLimitGeneLst = ujcDfr[
    ((ujcDfr['num_UJC_F_limited'] > 0) |
     (ujcDfr['num_UJC_M_limited'] > 0))
]['geneID']

sexBiasGeneLst = ujcDfr[
    ((ujcDfr['num_UJC_F_bias'] > 0)) |
    (ujcDfr['num_UJC_M_bias'] > 0)
]['geneID']

print(
    f"There are {len(sexBiasGeneLst)} unique genes with at least one biased UJC.")
print()

print(
    f"There are {len(sexLimitGeneLst)} unique genes with at least one sex-limited UJC.")

len(set(sexLimitGeneLst) - set(sexBiasGeneLst))
print(f"{len(set(sexLimitGeneLst) - set(sexBiasGeneLst))} of these genes are exclusively sex-limited.")


# Print a statistical summary of the num jxnHash per gene
print("Statistical summary of number of jxnHash per gene:")
print(ujcDfr['num_UJC'].describe())
print()

# The same but exclusively male/female
femaleDfr = ujcDfr[ujcDfr['flag_F_readCnt0'] == 1]
maleDfr = ujcDfr[ujcDfr['flag_M_readCnt0'] == 1]

print("Statistical summary of number of jxnHash per gene (at least 1 female read):")
print(femaleDfr['num_UJC'].describe())
print()

print("Statistical summary of number of jxnHash per gene (at least 1 male read):")
print(maleDfr['num_UJC'].describe())
print()

# Boxplot of num UJC!!
fig, axLst = plt.subplots(2, 1, figsize=(7, 7))

outlierLst = []

for dfr in [ujcDfr, femaleDfr, maleDfr]:
    data = dfr['num_UJC']
    q1 = data.quantile(0.25)
    q3 = data.quantile(0.75)
    iqr = q3 - q1
    lowerBoundary = q1 - 1.5 * iqr
    upperBoundary = q3 + 1.5 * iqr
    outlierLst.append(
        str(len(data[(data < lowerBoundary) | (data > upperBoundary)])))

xLabelLst = ['all\nNum Outliers: '+outlierLst[0], 'ujc with at least one\n female read\nNum Outliers: ' +
             outlierLst[1], 'ujc with at least one\n male read\nNum Outliers: '+outlierLst[2]]

axLst[0].boxplot(
    [ujcDfr['num_UJC'], femaleDfr['num_UJC'], maleDfr['num_UJC']])
axLst[0].set_xticklabels(xLabelLst)

axLst[1].boxplot([ujcDfr['num_UJC'], femaleDfr['num_UJC'],
                 maleDfr['num_UJC']], showfliers=False)
axLst[1].get_xaxis().set_ticks([])

axLst[1].set_title('No Outliers', y=-0.15)

plt.subplots_adjust(hspace=0.3)

fig.suptitle('Number of UJCs per Gene')

# fig.tight_layout()
# plt.show()
plt.savefig(f"{DIR}/figures/boxplot_num_ujc_per_gene.png",
            dpi=600, format="png", bbox_inches='tight')


# EVALUATING SEX BIAS IN GENES

ujcDfr.columns
genesWMixBiasSet = set()
# Print sex biased
for name, dfr in zip(['UJC', 'ERP'], [ujcDfr, erpDfr]):

    femaleBiasOrLimitLst = dfr[
        ((dfr['num_'+name+'_F_limited'] > 0)) |
        (dfr['num_'+name+'_F_bias'] > 0)
    ]['geneID']

    maleBiasOrLimitLst = dfr[
        ((dfr['num_'+name+'_M_limited'] > 0)) |
        (dfr['num_'+name+'_M_bias'] > 0)
    ]['geneID']

    mixedBiasLst = dfr[
        (((dfr['num_'+name+'_M_limited'] > 0)) |
         (dfr['num_'+name+'_M_bias'] > 0))
        &
        (((dfr['num_'+name+'_F_limited'] > 0)) |
         (dfr['num_'+name+'_F_bias'] > 0))
    ]['geneID']

    numBothBiasOrLimit = len(
        mixedBiasLst)

    genesWMixBiasSet.update(mixedBiasLst)

    print("There are", len(femaleBiasOrLimitLst),
          f"genes that have at least 1 {name} that is female-biased or female-limited.")
    print("There are", len(maleBiasOrLimitLst),
          f"genes that have at least 1 {name} that is male-biased or male-limited.")
    print("There are", numBothBiasOrLimit,
          f"genes that have a mix of both male/female-biased or male/female-limited {name}s.")

    allBiasedLimitedSet = set(
        femaleBiasOrLimitLst + maleBiasOrLimitLst + mixedBiasLst)
    print(
        f"There are {len(allBiasedLimitedSet)} unique genes with at least 1 {name} that are biased/limited in some way.")


# Do the same but for only one UJC/ERP in gene

oneUJCDfr = ujcDfr[ujcDfr['num_UJC'] == 1]
oneERPDfr = erpDfr[erpDfr['num_ERP'] == 1]

print(f"There are {len(oneUJCDfr)} genes with one UJC.")
print(f"There are {len(oneERPDfr)} genes with one ERP.")

for name, dfr in zip(['UJC', 'ERP'], [oneUJCDfr, oneERPDfr]):

    femaleBiasOrLimitLst = dfr[
        ((dfr['num_'+name+'_F_limited'] > 0)) |
        (dfr['num_'+name+'_F_bias'] > 0)
    ]['geneID'].tolist()

    maleBiasOrLimitLst = dfr[
        ((dfr['num_'+name+'_M_limited'] > 0)) |
        (dfr['num_'+name+'_M_bias'] > 0)
    ]['geneID'].tolist()

    mixedBiasLst = dfr[
        (((dfr['num_'+name+'_M_limited'] > 0)) |
         (dfr['num_'+name+'_M_bias'] > 0))
        &
        (((dfr['num_'+name+'_F_limited'] > 0)) |
         (dfr['num_'+name+'_F_bias'] > 0))
    ]['geneID'].tolist()

    print(f"{name}")

    print("There are", len(femaleBiasOrLimitLst),
          "genes: female-biased or female-limited (and expressed).")
    print("There are", len(maleBiasOrLimitLst),
          "genes: male-biased or male-limited (and expressed).")
    print("There are", len(mixedBiasLst),
          "genes: both male/female-biased or male/female-limited (and expressed).")

    allBiasedLimitedSet = set(
        femaleBiasOrLimitLst + maleBiasOrLimitLst + mixedBiasLst)
    print(
        f"There are {len(allBiasedLimitedSet)} unique genes with at least 1 {name} that are biased/limited in some way.")

    print()
    print()

    outDct = dict()
    for geneID in allBiasedLimitedSet:
        if geneID in mixedBiasLst:
            outDct[geneID] = "mixed_bias"
        elif geneID in maleBiasOrLimitLst:
            outDct[geneID] = "male_bias_or_limit"
        elif geneID in femaleBiasOrLimitLst:
            outDct[geneID] = "female_bias_or_limit"
        else:
            outDct[geneID] = "???"

    outDfr = pd.DataFrame(pd.Series(outDct)).rename_axis(
        'geneID').reset_index()
    outDfr = outDfr.rename({0: 'sex_bias'}, axis=1)
    outDfr.to_csv(
        f"{DIR}/supplementary/biased_genes_w_one_{name}.csv", index=False)

# Do the same but only for sex det genes (list of 22 curated genes)
sexDetFile = f"{DIR}/supplementary/sex_det.csv"
sexDetGnDfr = pd.read_csv(sexDetFile, low_memory=False)

sexDetUJCDfr = ujcDfr[ujcDfr['geneID'].isin(
    sexDetGnDfr['primary_fbgn'].tolist())]
sexDetERPDfr = erpDfr[erpDfr['geneID'].isin(
    sexDetGnDfr['primary_fbgn'].tolist())]

sexDetGn2Symbol = sexDetGnDfr.set_index('primary_fbgn')['symbol'].to_dict()

outDfrLst = []
for name, dfr in zip(['UJC', 'ERP'], [sexDetUJCDfr, sexDetERPDfr]):

    femaleBiasOrLimitLst = dfr[
        ((dfr['num_'+name+'_F_limited'] > 0)) |
        (dfr['num_'+name+'_F_bias'] > 0)
    ]['geneID'].tolist()

    maleBiasOrLimitLst = dfr[
        ((dfr['num_'+name+'_M_limited'] > 0)) |
        (dfr['num_'+name+'_M_bias'] > 0)
    ]['geneID'].tolist()

    mixedBiasLst = dfr[
        (((dfr['num_'+name+'_M_limited'] > 0)) |
         (dfr['num_'+name+'_M_bias'] > 0))
        &
        (((dfr['num_'+name+'_F_limited'] > 0)) |
         (dfr['num_'+name+'_F_bias'] > 0))
    ]['geneID'].tolist()

    print(f"{name}")

    print("There are", len(femaleBiasOrLimitLst),
          "genes: female-biased or female-limited (and expressed).")
    print("There are", len(maleBiasOrLimitLst),
          "genes: male-biased or male-limited (and expressed).")
    print("There are", len(mixedBiasLst),
          "genes: both male/female-biased or male/female-limited (and expressed).")

    allBiasedLimitedSet = set(
        femaleBiasOrLimitLst + maleBiasOrLimitLst + mixedBiasLst)
    print(
        f"There are {len(allBiasedLimitedSet)} unique genes with at least 1 {name} that are biased/limited in some way.")

    print()
    print()

    outDct = dict()
    for geneID in sexDetGn2Symbol.keys():
        if geneID in mixedBiasLst:
            outDct[geneID] = name+"_mixed_bias"
        elif geneID in maleBiasOrLimitLst:
            outDct[geneID] = name+"_male_bias_or_limit"
        elif geneID in femaleBiasOrLimitLst:
            outDct[geneID] = name+"_female_bias_or_limit"
        else:
            outDct[geneID] = "unbiased"

    outDfr = pd.DataFrame(pd.Series(outDct)).rename_axis(
        'geneID').reset_index()
    outDfr['symbol'] = outDfr['geneID'].apply(lambda x: sexDetGn2Symbol[x])
    outDfr = outDfr.rename({0: name+'_sex_bias'}, axis=1)
    outDfr = outDfr[['symbol', 'geneID', name+'_sex_bias']]

    outDfrLst.append(outDfr)

# merge is checked and all are both
mergeDfr = pd.merge(outDfrLst[0], outDfrLst[1], on=[
                    'symbol', 'geneID'], how='outer')

mergeDfr.columns = ['Gene Symbol',
                    'FlyBase Gene Identifier', 'UJC Bias', 'ERP Bias']

mergeDfr.to_csv(
    f"{DIR}/tables/results_table_sexdet_gene_bias.csv", index=False)
