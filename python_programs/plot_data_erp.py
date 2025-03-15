#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
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


def createERPRectanglePlot(dfr, ax, plotNum, boxWidth, boxHeight, xSpace, yScale):
    dfr = dfr.reset_index(names='figurePos')
    dfr['figurePos'] = dfr['figurePos'] + 1

    # Set up variables for coloring ERPs based on proportion of reads that are F
    colorLst = ['blue', 'purple', 'red']
    colorMap = LinearSegmentedColormap.from_list(
        "female_v_male_cmap", colorLst)
    normalize = Normalize(vmin=0, vmax=1)

    if not dfr.empty:

        for row in dfr.to_dict('records'):

            erpBinary = row['erpBinary']
            figurePos = row['figurePos']
            dataProp = row['prop_F']
            pVal = row['pVal']

            xCoord = 0

            for binaryFlag in erpBinary:
                fill = bool(binaryFlag)

                # Set the color of the rectangles and lines.
                # Red (100% female) to Purple (50/50) to Blue (100% male)
                if pd.isna(dataProp):
                    color = 'black'
                else:
                    color = colorMap(normalize(dataProp))

                ax.add_patch(Rectangle(
                    xy=(xCoord, (figurePos - boxHeight / 2)),
                    width=boxWidth,
                    height=boxHeight,
                    facecolor=color if fill else 'white',
                    edgecolor=color,
                    # alpha=alpha
                ))

                xCoord = xCoord + boxWidth + xSpace

            # y-axis
            # Set y-axis limits to go slightly past the number of UJC in the dfr
            ax.set_ylim(len(dfr) + yScale + 1, -yScale)

            # Make sure transcripts are set to their assigned y-position (=figurePos, which is based on sorted ERP)
            ax.set_yticks(dfr['figurePos'])

            # Use ERP as the label for y-axis
            ax.set_yticklabels(dfr['ERP_plus'])

            # Set y-axis label
            ax.set_ylabel("ERP")

            # turn off x-axis labels
            ax.xaxis.set_visible(False)

            # Create an object for the legend for the range of colors (will be added to plot later after everything is added)
            colorLegend = ScalarMappable(norm=normalize, cmap=colorMap)
            colorLegend.set_array([])

            finalX = xCoord + boxWidth/2

            # Figure index
            ax.text(
                finalX,  # x
                figurePos+yScale,  # y
                f"{plotNum}-" + str(figurePos) + "\n",  # text
                fontsize=10,
                ha='center',
                va='center'
            )

            # If no pval, mark with an x on the left of transcript
            if np.isnan(pVal):
                ax.plot(
                    -0.5,
                    figurePos,
                    marker="x",
                    color='black',
                )

            # If pval significant, mark with an star on the left of transcript
            elif pVal <= 0.05:
                ax.plot(
                    -0.5,
                    figurePos,
                    marker="*",
                    color='black',
                    markersize=10
                )
    else:
        ax.axis('off')
        ax.text(0.5, 0.5, '\n\n\nNO EXPRESSED ERPs', fontsize=20,
                ha='center', va='center', transform=ax.transAxes)

        colorLegend = ScalarMappable()
        colorLegend.set_array([])

    ax.relim()
    ax.autoscale_view()

    return colorLegend


# Test genes for dmel6
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

symbol = "Sxl"
# symbol = "dsx"
inputGn = geneDct[symbol]

DIR = "Z:/SHARE/McIntyre_Lab/kinfe/thesis/"
# DIR = "C:/Users/knife/Desktop/thesis/"

dataName = 'dmel_data'
genome = "dmel6"
dataFile = f"{DIR}/supplementary/datafile_erp.csv"
erFile = f"B:/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_er.gtf"
noAnnoThresh = 5
# outdir = "/nfshome/k.bankole/Desktop"

inDataFileDfr = pd.read_csv(dataFile, low_memory=False)
inERDfr = readExonData(erFile)

annoDfr = inDataFileDfr[inDataFileDfr['flag_in_full_anno'] == 1]
noAnnoDfr = inDataFileDfr[inDataFileDfr['flag_in_full_anno'] == 0]

annoDfr = annoDfr[annoDfr['geneID'] == inputGn].copy()
noAnnoDfr = noAnnoDfr[noAnnoDfr['geneID'] == inputGn].copy()
numTotalER = len(inERDfr[inERDfr['geneID'] == inputGn].copy())
strand = annoDfr['ERP_plus'].iloc[0][0]

plotInfoDfrLst = []

for dfr in [annoDfr, noAnnoDfr]:

    # dfr = annoDfr.copy()
    isAnno = False
    if dfr.equals(annoDfr):
        isAnno = True

    # Get average number of F and M read across desired samples
    dfr['total_F_read'] = dfr[dfr.filter(
        regex=("^erpCnts_.*_F_rep[4-6]$")).columns].sum(axis=1)
    dfr['avg_F_read'] = dfr['total_F_read'] / \
        len(dfr.filter(regex=("^erpCnts_.*_F_rep[4-6]$")).columns)

    dfr['total_M_read'] = dfr[dfr.filter(
        regex=("^erpCnts_.*_M_rep[4-6]$")).columns].sum(axis=1)
    dfr['avg_M_read'] = dfr['total_M_read'] / \
        len(dfr.filter(regex=("^erpCnts_.*_M_rep[4-6]$")).columns)

    dfr = dfr[dfr['total_F_read'] +
              dfr['total_M_read'] > 0].copy()
    if noAnnoThresh and not isAnno:
        dfr = dfr[dfr['total_F_read'] +
                  dfr['total_M_read'] > noAnnoThresh].copy()

    if not dfr.empty:
        # Create a proportion of female reads using above averages
        dfr['sumAvgRead'] = dfr['avg_F_read'] + dfr['avg_M_read']
        dfr['prop_F'] = dfr.apply(
            lambda row: row['avg_F_read'] / row['sumAvgRead'] if row['sumAvgRead'] > 0 else None, axis=1)

        plotInfoDfr = dfr[['ERP_plus', 'sumAvgRead', 'prop_F',
                           'total_F_read', 'total_M_read', 'pval_ttest_Equal']].copy().reset_index(drop=True)

        plotInfoDfr['pVal'] = plotInfoDfr['pval_ttest_Equal']

        plotInfoDfr['info'] = plotInfoDfr['ERP_plus'].str.split('|')

        plotInfoDfr[['ERP', 'flagDataOnlyExon', 'flagIR']] = pd.DataFrame(
            plotInfoDfr['info'].tolist(), index=plotInfoDfr.index)

        plotInfoDfr['erpBinary'] = plotInfoDfr['ERP'].apply(
            lambda x: list(int(i) for i in x.split('_')[1]))

        # Drop ERP if there are no reads only if there is no annotated ERP
        if not isAnno:
            plotInfoDfr = plotInfoDfr[plotInfoDfr['prop_F'].notna()]

        plotInfoDfr = plotInfoDfr.sort_values(
            by=['ERP_plus'], ignore_index=True, ascending=False)

        plotInfoDfrLst.append(plotInfoDfr)
    else:
        plotInfoDfrLst.append(pd.DataFrame())

# Set up figure/axes

totalNumERP = len(plotInfoDfrLst[0]) + len(plotInfoDfrLst[1])

# Make the height of each panel match the ratio of the number of ERPs each panel has
height_ratios = [len(plotInfoDfrLst[0]), len(plotInfoDfrLst[1])]

yScale = 0.3

fig, axLst = plt.subplots(2, 1, figsize=(
    numTotalER, max(totalNumERP*yScale, 5)), height_ratios=height_ratios)

plotNum = 0
for plotInfoDfr, ax in zip(plotInfoDfrLst, axLst):

    ax.set_frame_on(False)  # Removes the frame

    # Plot rectangles for anno and noAnno ERPs
    plotNum += 1
    colorLegend = (
        createERPRectanglePlot(
            dfr=plotInfoDfr,
            plotNum=plotNum,
            ax=ax,
            boxWidth=1,
            boxHeight=0.8,
            xSpace=0.2,
            yScale=yScale
        )
    )

# Add the indicator of 5'/3' depending on strand
fontSize = 0.06460296096904442 * \
    ((numTotalER * max(totalNumERP*yScale, 5))) + 22.70794078061911

fig.text(
    0.15,  # x
    1,  # y
    ("5'" if strand == "+" else "3'"),  # text
    fontsize=fontSize,
    va='center',
    ha='center',
    transform=fig.transFigure
)

fig.text(
    1,  # x
    1,  # y
    ("3'" if strand == "+" else "5'"),  # text
    fontsize=fontSize,
    va='center',
    ha='center',
    transform=fig.transFigure
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

if noAnnoThresh:
    axTitleLst = [firstLine +
                  f"Total Number of ERPs: {totalNumERP}\n\n"

                  'Annotated ERPs\n'
                  f'Number of ERPs: {len(plotInfoDfrLst[0])}',

                  f'Unannotated Data ERPs\n'
                  f'with num reads > {noAnnoThresh}\n'
                  f'Number of ERPs: {len(plotInfoDfrLst[1])}\n']

else:
    axTitleLst = [firstLine +
                  f"Total Number of ERPs: {totalNumERP}\n\n"

                  'Annotated Data ERPs\n'
                  f'Number of ERPs: {len(plotInfoDfrLst[0])}',

                  f'Unannotated Data ERPs\n'
                  f'with num reads > {noAnnoThresh}\n'
                  f'Number of ERPs: {len(plotInfoDfrLst[1])}\n']

for index in range(2):
    axLst[index].set_title(axTitleLst[index])

fig.subplots_adjust(hspace=0.7, wspace=0.5)

plotOutFile = (
    f"{DIR}/figures/"
    f"plot_data_erpRectangle_w_noAnno_rc{noAnnoThresh}"
    f"_read_{inputGn}_{symbol}_{dataName}_2_{genome}.svg"
)

fig.savefig(plotOutFile, dpi=600, format="svg", bbox_inches='tight')
print(f"Saved plot: {plotOutFile}")

# fig.tight_layout()
# plt.show()
plt.close()

# Figure Index
figIndexDct = {'geneID': [], 'ERP_plus': [],
               'figurePos': [], 'total_F_read': [], 'total_read': []}

loopNum = 0
for dfr in plotInfoDfrLst:

    loopNum += 1
    if dfr.empty:
        continue
    else:
        erpNum = 0
        for row in dfr.to_dict('records'):
            erpNum += 1
            figIndexDct['geneID'].append(inputGn)
            figIndexDct['ERP_plus'].append(row['ERP_plus'])
            figIndexDct['figurePos'].append(
                f'{loopNum}-{erpNum}')
            figIndexDct['total_F_read'].append(row['total_F_read'])
            figIndexDct['total_read'].append(
                row['total_F_read'] + row['total_M_read'])

figIndxOutFile = (
    f"{DIR}/supplementary/"
    f"table_data_erpRectangle_w_noAnno_rc{noAnnoThresh}"
    f"_read_{inputGn}_{symbol}_{dataName}_2_{genome}_figureIndex.csv"
)

figIndexDfr = pd.DataFrame(figIndexDct)
figIndexDfr.to_csv(figIndxOutFile, index=False)
print(f"Saved figure index: {figIndxOutFile}")
