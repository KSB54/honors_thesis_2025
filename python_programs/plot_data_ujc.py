#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 15:42:11 2024

@author: k.bankole
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import ScalarFormatter
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable


def readExonData(gtfFile):
    # Read in a GTF. Create geneID and transcriptID column based on attributes.

    columns = ['seqname', 'source', 'feature', 'start',
               'end', 'score', 'strand', 'frame', 'attribute']

    data = pd.read_csv(gtfFile, sep='\t', comment='#',
                       header=None, names=columns,
                       low_memory=False)

    data = data[data['feature'] == "exon"]

    data['start'] = data['start'].astype(int)
    data['end'] = data['end'].astype(int)
    data['geneID'] = data['attribute'].str.extract(r'gene_id "([^"]+)"')
    data['transcriptID'] = data['attribute'].str.extract(
        r'transcript_id "([^"]+)"')

    data.drop(columns=['score', 'frame', 'attribute'], inplace=True)

    return data


def createDataSplitPanel(dfr, panel, panelNum, minStart, maxEnd, wholeFigureXLim, longestModelLngth, pValCol, colorLst, exonHeight, alpha, includeFigIndex=False):

    # If the panel is empty, display NO UJCs
    if dfr.empty:

        panel.axis('off')
        panel.text(0.5, 0.5, '\nNO UJCS', fontsize=20,
                   ha='center', va='center', transform=panel.transAxes)
        panel.set_xlim(wholeFigureXLim)

        colorLegend = ScalarMappable()
        colorLegend.set_array([])
        xPadding = 0.05 * longestModelLngth
    # Otherwise
    else:
        # Create a figure y-position for jxnHash based on sorted ERP/jxnHash
        dfr = dfr.sort_values(by=['ERP_plus', 'jxnHash'], ignore_index=True)
        dfr = dfr.reset_index(names='figurePos')
        dfr['figurePos'] = dfr['figurePos'] + 1

        # Scale y padding based on the number of UJC in gene
        # Scale x padding based on longest UJC/geneModel
        xPadding = 0.05 * longestModelLngth
        yPadding = .005 * len(dfr)

        # Loop through every UJC in dataframe
        for row in dfr.to_dict('records'):

            # Extract variables from row in dataframe
            figurePos = row['figurePos']
            exonLst = row['coords']

            dataProp = row['prop_F']
            pVal = row[pValCol]

            # Set up variables for coloring UJCs based on proportion of reads that are F
            colorMap = LinearSegmentedColormap.from_list(
                "female_v_male_cmap", colorLst)
            normalize = Normalize(vmin=0, vmax=1)

            # Loop through every exon for that UJC
            prevEnd = None
            for exon in exonLst:
                start = exon[0]
                end = exon[1]

                # Set the color of the rectangles and lines.
                # Red (100% female) to Purple (50/50) to Blue (100% male)
                if dataProp == np.nan:
                    color = 'black'
                else:
                    color = colorMap(normalize(dataProp))

                # Add a rectangle to figure for every exon using start/end pos
                panel.add_patch(Rectangle(
                    xy=(start, figurePos - exonHeight / 2),
                    width=end - start,
                    height=exonHeight,
                    facecolor=color,
                    edgecolor=color,
                    alpha=alpha
                ))

                # After the first exon, plot a line from the start of the last exon to the beginning of the first
                if prevEnd is not None:
                    panel.plot(
                        [prevEnd, start],
                        [figurePos, figurePos],
                        color=color,
                        alpha=alpha
                    )

                prevEnd = end

            # If a figure index is being output -> add text to the right of
            # transcripts with the figure position that matches what it will
            # look like in the index
            if includeFigIndex:
                panel.text(
                    (maxEnd + .5 * xPadding),  # x
                    figurePos,  # y
                    f"{panelNum}-" + str(figurePos),  # text
                    fontsize=10,
                    ha='center',
                    va='center'
                )

            # If no pval, mark with an x on the left of transcript
            if np.isnan(pVal):
                panel.plot(
                    minStart - (.5 * xPadding),
                    figurePos,
                    marker="x",
                    color='black',
                    alpha=alpha,
                )

            # If pval significant, mark with an star on the left of transcript
            elif pVal < 0.05:
                panel.plot(
                    minStart - (.5 * xPadding),
                    figurePos,
                    marker="*",
                    color='black',
                    alpha=alpha,
                    markersize=10
                )

        # Create an object for the legend for the range of colors (will be added to plot later after everything is added)
        colorLegend = ScalarMappable(norm=normalize, cmap=colorMap)
        colorLegend.set_array([])

        # y-axis
        # Set y-axis limits to go slightly past the number of UJC in the dfr
        panel.set_ylim(len(dfr) + yPadding + 1, -yPadding)

        # Make sure transcripts are set to their assigned y-position (=figurePos, which is based on sorted ERP)
        panel.set_yticks(dfr['figurePos'])

        # Use ERP as the label for y-axis
        panel.set_yticklabels(dfr['ERP_plus'])

        # Set y-axis label
        panel.set_ylabel("ERP")

        # x-axis
        # Set x-axis limits
        panel.set_xlim(wholeFigureXLim)

        # Assure that scientific notation is not used and that the x-axis can fit labels
        panel.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        panel.ticklabel_format(style="plain", axis='x')

        # Set x-axis label
        panel.set_xlabel(
            "Genomic Position")

        return colorLegend, xPadding


def createGeneModelPanel(erCoordLst, esCoordLst, geneModelPanel,
                         exonHeight, alpha, color, edgecolor, segmentColor,
                         xlim, xPadding, strand):
    prevEnd = None
    # Loop through every ER
    for er in erCoordLst:

        start = er[0]
        end = er[1]

        # Add a rectangle to figure for every exon region using start/end pos
        geneModelPanel.add_patch(Rectangle(
            xy=(start, 1 - exonHeight / 2),
            width=end - start,
            height=exonHeight,
            facecolor=color,
            edgecolor=edgecolor,
            alpha=alpha
        ))

        # After the first exon, plot a line from the start of the last exon to the beginning of the first
        if prevEnd is not None:
            geneModelPanel.plot([prevEnd, start], [1, 1],
                                color=color, alpha=alpha)

        prevEnd = end

    # Loop through every ES and plot a vertical line at the end of every segment
    # (except the last, which is just the end of the gene)
    for es in esCoordLst[:-1]:
        end = es[1]
        geneModelPanel.vlines(
            x=end,
            ymin=1 - exonHeight / 2,
            ymax=1 + exonHeight / 2,
            color=segmentColor
        )

    # y-axis
    # Set y-axis limits to large enough to fit gene model
    geneModelPanel.set_ylim(1 - exonHeight - .3, 1 + exonHeight + .3)

    # Make y-axis label blank
    geneModelPanel.axes.get_yaxis().set_ticks([])

    # x-axis
    # Set x-axis limits
    geneModelPanel.set_xlim(xlim)

    # Assure that scientific notation is not used and that the x-axis can fit labels
    geneModelPanel.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    geneModelPanel.ticklabel_format(style="plain", axis='x')

    # Set x-axis label
    geneModelPanel.set_xlabel("Genomic Position")

    # Title this panel of the figure
    geneModelPanel.title.set_text(
        "Gene Model Exon Region and Segments\n"
        f"Number of Exon Regions: {len(erCoordLst)} "
        f"Number of Exon Segments: {len(esCoordLst)}"
    )

    # If it is the first exon, add the indicator of 5'/3' depending on strand
    geneModelPanel.text(
        xlim[0] - .5 * xPadding,  # x
        1,  # y
        ("5'" if strand == "+" else "3'"),  # text
        fontsize=20,
        va='center',
        ha='center',
    )

    # If it is the last exon (loop is over), add the indicator of 5'/3' depending on strand
    geneModelPanel.text(
        xlim[1] + .5 * xPadding,  # x
        1,  # y
        ("3'" if strand == "+" else "5'"),  # text
        fontsize=20,
        va='center',
        ha='center',
    )


# List of dmel genes that can be plot
geneDct = {
    'tra': 'FBgn0003741',
    'tra2': 'FBgn0003742',
    'mub': 'FBgn0262737',
    'her': 'FBgn0001185',
    'fru': 'FBgn0004652',
    'fl(2)d': 'FBgn0000662',
    'Sxl': 'FBgn0264270',
    'spf45': 'FBgn0086683',
    'Psi': 'FBgn0014870',
    'dsx': 'FBgn0000504'
}

DIR = "Z:/SHARE/McIntyre_Lab/kinfe/thesis/"
# DIR = "C:/Users/knife/Desktop/thesis/"

symbol = "Sxl"
inputGn = geneDct[symbol]
genome = "dmel6"
dataName = "dmel_data"

dataGTFFile = f"B:/rmg_lmm_dros_data/dmel_data_2_dmel6_ujc_noMultiGene.gtf"
dataAnnoFile = f"{DIR}/supplementary/datafile_ujc.csv"
dataERPFile = f"{DIR}/supplementary/datafile_erp.csv"
dataCvgFile = f"Z:/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary/datafiles/prop_gene_cvrg_datafile_dmel6.csv"

erFile = f"B:/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_er.gtf"
esFile = f"B:/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_es.gtf"

dataInfoERPFile = f"B:/sex_specific_splicing/rmg_erp_output/fiveSpecies_2_dmel6_ujc_er_vs_dmel_data_2_dmel6_ujc_noMultiGene_infoERP.csv"

dataCol = ['rawCnts_dmel_F_rep4', 'rawCnts_dmel_F_rep5', 'rawCnts_dmel_F_rep6',
           'rawCnts_dmel_M_rep4', 'rawCnts_dmel_M_rep5', 'rawCnts_dmel_M_rep6']

pValCol = 'pval_ttest_Equal'
threshType = 'num'  # num or pct, only two options

# optional, default false
includeFigIndex = True

# num or pct (0-1)
cvgThreshWERP = 5
cvgThreshNoERP = 5

# Read in all input files

print("Reading files...")
inDataGtfDfr = readExonData(dataGTFFile)
inDataAnnoDfr = pd.read_csv(dataAnnoFile, low_memory=False)
inDataERPDfr = pd.read_csv(dataERPFile, low_memory=False)
inDataCvgDfr = pd.read_csv(dataCvgFile, low_memory=False)
inERDfr = readExonData(erFile)
inESDfr = readExonData(esFile)
inDataInfoERPDfr = pd.read_csv(dataInfoERPFile, low_memory=False)

print("Processing...")
# Subset every input to the desired input gene
# TODO: check that these aren't empty
dataGTFDfr = inDataGtfDfr[inDataGtfDfr['geneID'] == inputGn].copy().rename(
    columns={'transcriptID': 'jxnHash'})

print(f"Number of data UJC in {inputGn}:", dataGTFDfr['jxnHash'].nunique())

annoDfr = inDataAnnoDfr[inDataAnnoDfr['geneID'] == inputGn].copy()
dataInfoERPDfr = inDataInfoERPDfr[inDataInfoERPDfr['geneID'] == inputGn].copy(
)
cvgDfr = inDataCvgDfr[inDataCvgDfr['geneID'] ==
                      inputGn].copy().drop(inDataCvgDfr.columns[0], axis=1)
erDfr = inERDfr[inERDfr['geneID'] == inputGn].copy()
esDfr = inESDfr[inESDfr['geneID'] == inputGn].copy()

erpData = inDataERPDfr[inDataERPDfr['geneID'] == inputGn].copy()

# Grab strand
strand = dataGTFDfr['strand'].iloc[0]

# make gtf unique on jxnHash based on GTF with exon coords stored
dataGTFDfr['coords'] = dataGTFDfr.apply(
    lambda row: (row['start'], row['end']), axis=1)
uniqJxnHashDfr = dataGTFDfr.groupby(['geneID', 'jxnHash']).agg(
    list)['coords'].reset_index().sort_values(by='jxnHash')

# Merge in ERP_plus to dataGTFDfr
print(uniqJxnHashDfr.columns)
uniqJxnHashDfr = pd.merge(uniqJxnHashDfr, dataInfoERPDfr[[
                          'jxnHash', 'ERP_plus']], on='jxnHash', how='outer', indicator='merge_check')

if (uniqJxnHashDfr['merge_check'] != "both").all():
    raise Exception(
        "There are jxnHash in the UJC GTF that are not in the infoERP or vice versa...")
else:
    uniqJxnHashDfr.drop('merge_check', axis=1, inplace=True)

# Based on the type of threshold used, compare flag each noAnno jxnHash for
# whether or not it is >= to threshold in either M or F
if threshType == 'num':
    compCol = "ujc_rawCnts"
elif threshType == 'pct':
    compCol = "prop_gene_coverage"

cvgDfr['aboveERPThresh'] = cvgDfr.apply(
    lambda row:
    row[compCol + '_M'] >= cvgThreshWERP or
    row[compCol + '_F'] >= cvgThreshWERP, axis=1)

cvgDfr['aboveNoERPThresh'] = cvgDfr.apply(
    lambda row:
    row[compCol + '_M'] >= cvgThreshNoERP or
    row[compCol + '_F'] >= cvgThreshNoERP, axis=1)


# split datafile:
#      annoDfr_wAnno:  where flag_jxnHash_in_fiveSpecies_full_anno = 1
#      annoDfr_noAnno: where flag_jxnHash_in_fiveSpecies_full_anno = 1
print(erpData.columns)
print(annoDfr.columns)

annoDfr_wAnno = annoDfr[annoDfr['flag_jxnHash_in_fiveSpecies_full_anno'] == 1].copy()
annoDfr_noAnno = annoDfr[annoDfr['flag_jxnHash_in_fiveSpecies_full_anno'] == 0].copy()
print(len(annoDfr))
print(len(annoDfr_noAnno))

# Flag all annotated jxnHash for keeping
if not annoDfr_wAnno.empty:
    annoDfr_wAnno.loc[:, 'keepHash'] = True

# keep only needed cols from ERP datafile
keep = ['ERP_plus', 'flag_in_full_anno', 'geneID']
add_ERP = erpData[keep].drop_duplicates()
before = len(erpData)
after = len(add_ERP)
print(f"before = {before}, after = {after}")  # 84 both

# merge to data without Anno
noAnno_w_erp_flags = pd.merge(
    annoDfr_noAnno,
    add_ERP,
    on=['ERP_plus', 'geneID'],
    how="outer",
    indicator='merge_check'
)
print(noAnno_w_erp_flags['merge_check'].value_counts(
    dropna=False).sort_index())

wERPDfr = noAnno_w_erp_flags[(noAnno_w_erp_flags['flag_in_full_anno'] == 1) & (
    noAnno_w_erp_flags['merge_check'] == 'both')]
wERPDfr = wERPDfr.drop(columns=['merge_check'])

noERPDfr = noAnno_w_erp_flags[(noAnno_w_erp_flags['flag_in_full_anno'] != 1) & (
    noAnno_w_erp_flags['merge_check'] == 'both')]
noERPDfr = noERPDfr.drop(columns=['merge_check'])

if not wERPDfr.empty:
    # Flag noAnno jxnHash that are above the threshold for keeping
    wERPDfr['keepHash'] = wERPDfr.apply(
        lambda row: cvgDfr.loc[cvgDfr['jxnHash'] == row['jxnHash'], 'aboveERPThresh'].iloc[0], axis=1)

if not noERPDfr.empty:
    noERPDfr['keepHash'] = noERPDfr.apply(
        lambda row: cvgDfr.loc[cvgDfr['jxnHash'] == row['jxnHash'], 'aboveNoERPThresh'].iloc[0], axis=1)

# Loop through all three dataframes and
# merge in coords from GTF and add a prop_F column calculated from the dataCol

dfrWCoordsLst = []
for dfr in [annoDfr_wAnno, wERPDfr, noERPDfr]:

    if not dfr.empty:

        dfr = dfr[dfr['keepHash'] == True]

        try:
            dfr = dfr[['geneID', 'jxnHash', pValCol] + dataCol].copy()
            dfr[pValCol] = dfr[pValCol].apply(
                lambda x: 0 if x == '<.0001' else float(x))

        except KeyError as e:

            dfr = dfr[['geneID', 'jxnHash'] + dataCol].copy()
            dfr[pValCol] = np.nan

            if pValCol in str(e):

                print(f"WARNING: A DATAFILE DOES NOT HAVE THE INPUT PVAL COLUMN ({pValCol})\n"
                      "UJCS WILL BE TREATED AS IF THEY DO NOT HAVE A PVAL")
            else:
                raise

        # Get average number of F read across desired samples
        dfr['total_F_read'] = dfr[dfr.filter(
            like='_F_rep').columns].sum(axis=1)
        dfr['avg_F_read'] = dfr['total_F_read'] / \
            len(dfr.filter(like='_F_rep').columns)

        # Get average number of M read across desired samples
        dfr['total_M_read'] = dfr[dfr.filter(
            like='_M_rep').columns].sum(axis=1)
        dfr['avg_M_read'] = dfr['total_M_read'] / \
            len(dfr.filter(like='_M_rep').columns)

        # Create a proportion of female reads using above averages
        dfr['sumAvgRead'] = dfr['avg_F_read'] + \
            dfr['avg_M_read']

        dfr['prop_F'] = dfr.apply(
            lambda row: row['avg_F_read'] / row['sumAvgRead'] if row['sumAvgRead'] > 0 else None, axis=1)

        wPropDfr = dfr[['geneID', 'jxnHash', pValCol,
                        'prop_F']].copy().reset_index(drop=True)

        # Merge in coords from uniq on jxnHash GTF info Dfr
        mergeDfr = pd.merge(uniqJxnHashDfr[['jxnHash', 'coords', 'ERP_plus']],
                            wPropDfr, on='jxnHash', how='outer', indicator='merge_check')

        if (mergeDfr['merge_check'] == 'right_only').any():
            raise Exception("Merge error, there are jxnHash that are only in the "
                            "datafile (and not in the data UJC GTF).")
        else:
            wCoordsDfr = mergeDfr[mergeDfr['merge_check'] == 'both'].copy()
            wCoordsDfr = wCoordsDfr.drop('merge_check', axis=1)

        # Remove any data UJC that have no reads
        # only occurs if dataCol is not all reps
        wCoordsDfr = wCoordsDfr.dropna(subset=['prop_F'])
        dfrWCoordsLst.append(wCoordsDfr)
    else:
        dfrWCoordsLst.append(pd.DataFrame())

# Count the number of data UJCs that have been kept (are above threshold,
# have read counts)
totalUJCAboveThresh = sum([dfr['jxnHash'].nunique()
                          for dfr in dfrWCoordsLst if not dfr.empty])

# Check that the number of UJCs in the dfrWCoords match the amount that should
# have been kept based on thresholds

if threshType == 'num':
    threshStringWERP = f"{cvgThreshWERP} reads"
    threshStringNoERP = f"{cvgThreshNoERP} reads"

elif threshType == 'pct':
    threshStringWERP = f"{cvgThreshWERP:0.2%} of reads in the gene"
    threshStringNoERP = f"{cvgThreshNoERP:0.2%} of reads in the gene"

print(f"Number of data UJC in {inputGn} that "
      "have reads in the input dataCol, have at least "
      + threshStringWERP + " for unannotated UJC with annotated ERP and at least "
      + threshStringNoERP + " for unannotated UJC with no annotated ERP:",
      totalUJCAboveThresh)

# Get exon regions and exon segments for gene model panel
# make er and es gtf unique on geneID based on GTF with exons stored
erDfr['ER_coords'] = erDfr.apply(
    lambda row: (row['start'], row['end']), axis=1)
geneERDfr = erDfr.groupby(['geneID']).agg(
    list)['ER_coords'].reset_index()

esDfr['ES_coords'] = esDfr.apply(
    lambda row: (row['start'], row['end']), axis=1)
geneESDfr = esDfr.groupby(['geneID']).agg(
    list)['ES_coords'].reset_index()

# Create dict with ER_coords and ES_coords
geneERESCoordDct = pd.merge(geneERDfr, geneESDfr, on='geneID',
                            how='outer').set_index('geneID').iloc[0].to_dict()

# Use dict to create a list of ERs and ESs (use for plotting gene model)
erCoordLst = geneERESCoordDct['ER_coords']
esCoordLst = geneERESCoordDct['ES_coords']

# Setup Plot and Panels

# Make the height of each panel match the ratio of the number of UJCs each panel has
height_ratios = [
    len(dfrWCoordsLst[0]),
    len(dfrWCoordsLst[1]),
    len(dfrWCoordsLst[2]),]

# Height for the gene model panel, min 1, max .75*the smallest ratio from above
height_ratios.append(max([1, 0.75 * min(height_ratios)]))

# Create a 4 panel figure (4 row one column) using the height ratios
fig, axLst = plt.subplots(4, 1, figsize=(
    15, 15), height_ratios=height_ratios)

# Create the limits for the x-axis of the panel
# Determine the minimum start and maximum end for every transcript across all dataframes
minStart = float('inf')
maxEnd = float('-inf')
for dfr in dfrWCoordsLst:
    if not dfr.empty:
        allCoordLst = sum(dfr['coords'], [])
        allStartLst = [x[0] for x in allCoordLst]
        allEndLst = [x[1] for x in allCoordLst]

        minStart = min([minStart] + allStartLst)
        maxEnd = max([maxEnd] + allEndLst)

# Check and see if the gene model has an earlier start/later end than any of the transcripts
geneModelStart = min([x[0] for x in erCoordLst])
geneModelEnd = max([x[1] for x in erCoordLst])

# True minimum and maximum after including gene model
wholeMinStart = min([minStart, geneModelStart])
wholeMaxEnd = max([maxEnd, geneModelEnd])

# Longest tr/gene model
longestModelLngth = wholeMaxEnd - wholeMinStart

# Calcluate limits of the x axis based on min start/min end and longest model
wholeFigureXLim = (wholeMinStart - 0.05 * longestModelLngth,
                   wholeMaxEnd + 0.05 * longestModelLngth)

# Loop through every dfr in the list (anno, wERP, noERP) and plot them on their separate panels
for index, dfr in enumerate(dfrWCoordsLst):

    if dfr.empty:
        createDataSplitPanel(
            dfr=dfr,
            panel=list(axLst)[index],
            panelNum=index + 1,
            minStart=minStart,
            maxEnd=maxEnd,
            longestModelLngth=longestModelLngth,
            wholeFigureXLim=wholeFigureXLim,
            pValCol=pValCol,
            colorLst=['blue', 'purple', 'red'],
            exonHeight=0.6,
            alpha=0.6,
            includeFigIndex=includeFigIndex
        )
    else:
        # Create and number the panel
        colorLegend, xPadding = (
            createDataSplitPanel(
                dfr=dfr,
                panel=list(axLst)[index],
                panelNum=index + 1,
                minStart=minStart,
                maxEnd=maxEnd,
                longestModelLngth=longestModelLngth,
                wholeFigureXLim=wholeFigureXLim,
                pValCol=pValCol,
                colorLst=['blue', 'purple', 'red'],
                exonHeight=0.6,
                alpha=0.6,
                includeFigIndex=includeFigIndex
            )
        )

# Plot gene model underneath data panels (axLst index 3)
createGeneModelPanel(
    erCoordLst=erCoordLst,
    esCoordLst=esCoordLst,
    geneModelPanel=axLst[3],
    exonHeight=0.6,
    alpha=0.6,
    color='lightblue',
    edgecolor='navy',
    segmentColor='black',
    xlim=wholeFigureXLim,
    xPadding=xPadding,
    strand=strand
)

# Create separate panel for color bar and move it into a good position
# left, bottom, width, height
cbPanel = fig.add_axes([.6, 0.175, 0.5, 0.715], frameon=False)
cbPanel.set_xticks([])
cbPanel.set_yticks([])
colorBar = plt.colorbar(colorLegend, ax=cbPanel, orientation="vertical")
labelText = ("Proportion of Female Reads in Data\n"
             "(0 = 100% M, 1 = 100% F,\n "
             "Black = Not Supported by Data)")
colorBar.set_label(labelText)
# Center align the label text
colorBar.ax.xaxis.label.set_horizontalalignment('center')

# Title each panel (title for whole figure is in the title for the first panel)
if symbol:
    firstLine = f"{genome} - {inputGn} - {symbol} (strand={strand})\n"
else:
    firstLine = f"{genome} - {inputGn} (strand={strand})\n"

if threshType == 'num':
    threshStringWERP = f"{cvgThreshWERP}"
    threshStringNoERP = f"{cvgThreshNoERP}"

elif threshType == 'pct':
    threshStringWERP = f"{cvgThreshWERP:0.2%} of"
    threshStringNoERP = f"{cvgThreshNoERP:0.2%} of"

panelTitleLst = [firstLine +
                 f"Total Number of UJCs: {totalUJCAboveThresh}\n\n"

                 'Annotated Data UJCs\n'
                 f'Number of UJCs: {len(dfrWCoordsLst[0])}',

                 f'Unannotated Data UJCs with annotated ERPs that have >= {threshStringWERP} reads in the gene\n'
                 f'Number of UJCs: {len(dfrWCoordsLst[1])}',

                 f'Unannotated Data UJCs with no annotated ERPs that have >= {threshStringNoERP} reads in the gene\n'
                 f'Number of UJCs: {len(dfrWCoordsLst[2])}']

for index in range(3):
    axLst[index].set_title(panelTitleLst[index])

# For dev, to display plots in IDE without saving
# fig.tight_layout()
# plt.show()

# Output files!!!

if symbol:
    symbolSection = f"{inputGn}_{symbol}"
else:
    symbolSection = f"{inputGn}"

if threshType == 'num':
    tFileStringWERP = f"{int(cvgThreshWERP)}"
    tFileStringNoERP = f"{int(cvgThreshNoERP)}"
elif threshType == 'pct':
    tFileStringWERP = f"{int(cvgThreshWERP * 100)}"
    tFileStringNoERP = f"{int(cvgThreshNoERP * 100)}"

plotOutFile = (f"{DIR}/figures/"
               f"plot_data_ujc_erpLabel_splitPanel_w_{threshType}Thresh"
               f"_wERP_{tFileStringWERP}_noERP_{tFileStringNoERP}"
               f"_{symbolSection}_{dataName}_2_{genome}.svg")

fig.subplots_adjust(hspace=0.4, wspace=0.3)
fig.savefig(plotOutFile, dpi=600, format="svg", bbox_inches='tight')
print(f"Saved plot: {plotOutFile}")

plt.close(fig)

# If figure index is output:
if includeFigIndex:
    figIndxOutFile = (f"{DIR}/supplementary/"
                      f"table_data_ujc_erpLabel_splitPanel_w_{threshType}Thresh"
                      f"_wERP_{tFileStringWERP}_noERP_{tFileStringNoERP}"
                      f"_{symbolSection}_{dataName}_2_{genome}_figureIndex.csv")

    figIndexDct = {'geneID': [], 'jxnHash': [], 'figurePos': []}

    loopNum = 0
    for dfr in dfrWCoordsLst:

        loopNum += 1
        if dfr.empty:
            continue
        else:

            # Create a figure y-position for jxnHash based on sorted ERP
            dfr = dfr.sort_values(
                by=['ERP_plus', 'jxnHash'], ignore_index=True)
            dfr = dfr.reset_index(names='figurePos')
            dfr['figurePos'] = dfr['figurePos'] + 1

            for row in dfr.to_dict('records'):
                figIndexDct['geneID'].append(row['geneID'])
                figIndexDct['jxnHash'].append(row['jxnHash'])
                figIndexDct['figurePos'].append(
                    f'{loopNum}-{row["figurePos"]}')

    figIndexDfr = pd.DataFrame(figIndexDct)
    figIndexDfr.to_csv(figIndxOutFile, index=False)
    print(f"Saved figure index: {figIndxOutFile}")
