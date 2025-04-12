import base64
import io
import logomaker
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import seaborn as sns
import sys
from wordcloud import WordCloud


# Figure parameters
labelSizeTitle = 18
labelSizeAxis = 16
labelSizeTicks = 13
lineThickness = 1.5
tickLength = 4
figSize = (12, 9)
figBorders = [0.882, 0.075, 0.05, 0.98]

# Experimental parameters
AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
      'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# Set matplotlib backend as non-interactive
matplotlib.use('Agg')


def processData(substrates, entropyMin, NSelect, enzymeName, defaultSubs):
    if defaultSubs:
        subLen = len(substrates[0])

        # Count: Substrates
        subsCounts = {}
        for sub in substrates:
            if len(sub) == subLen:
                keepSub = True
                for aa in sub:
                    if aa not in AA:
                        keepSub = False
                        break
                if keepSub:
                    if sub in subsCounts.keys():
                        subsCounts[sub] += 1
                    else:
                        subsCounts[sub] = 1
        substrates = subsCounts
    else:
        subLen = len(next(iter(substrates)))

    # Sort the dictionary by counts from highest to lowest
    substrates = dict(sorted(substrates.items(), key=lambda item: item[1], reverse=True))

    # Define: Substrate positions
    pos = [f'R{index + 1}' for index in range(subLen)]

    # Count AAs
    totalSubs = 0
    countedAA = pd.DataFrame(0, index=AA, columns=pos)
    for sub, count in substrates.items():
        totalSubs += count
        for index, aa in enumerate(sub):
            countedAA.loc[aa, pos[index]] += count

    # Evaluate: AA probability
    probAA = countedAA / totalSubs

    # Figure: Entropy
    figProb = plotProbabilities(probAA, totalSubs, enzymeName)
    entropy, entropyMax, figEntropy = plotEntropy(probAA, AA, enzymeName)
    figLogo = plotWeblogo(probAA, entropy, entropyMax, totalSubs, enzymeName)
    NBinSubs, figBarCounts, figBarProb, figWords = binSubstrates(
        substrates, entropy, entropyMin, NSelect, enzymeName)

    # Create dataset
    dataset = {}
    dataset['N'] = totalSubs
    dataset['NBinSubs'] = NBinSubs
    dataset['probability'] = figProb
    dataset['entropy'] = figEntropy
    dataset['pLogo'] = figLogo
    dataset['barCounts'] = figBarCounts
    dataset['barProb'] = figBarProb
    dataset['words'] = figWords

    return dataset



def subsDefault():
    substrates = {
        'VVLQAGTK': 19076,
        'VILQSVGA': 18383,
        'LALQSACW': 17331,
        'LNLQGILD': 16201,
        'IYLQALMP': 16010,

        'TSLQARKS': 15449,
        'AFLQAHFT': 15175,
        'VTLQCTYS': 14330,
        'VLLQAKQL': 14201,
        'LVLQANPC': 14010,

        'IMLQGVIW': 13449,
        'AGLQASAH': 12175,
        'KHLQSENE': 11330,
        'MVLQGDVN': 10201,
        'YGLQCNEV': 10010,

        'VCMQCQVQ': 9449,
        'GELQSWHF': 9005,
        'VDMQCMWG': 7330,
        'VWMQCSII': 7013,
        'VVMQCCSM': 6276,

        'VLIQCCPR': 1991,
        'VVMQSGSM': 1754,
        'VVFQCHNR': 1576,
        'KRFQCKLR': 1234,
        'PEFQCCLQ': 900
        }

    return substrates



def createCustomColorMap(colorType):
    if colorType == 'Probability':
        colors = ['#FFFFFF', '#ABFF9B', '#39FF14', '#2E9418', '#2E9418', '#005000']
    elif colorType == 'Word Cloud':
        colors = ['#A7A7A7', '#39FF14', '#2E9418', 'black']
    else:
        colors = ['#008631', '#39E75F','#CC5500', '#F79620', 'black']

    # Create colormap
    if len(colors) == 1:
        colorList = [(0, colors[0]), (1, colors[0])]
    else:
        colorList = [(i / (len(colors) - 1), color) for i, color in enumerate(colors)]
        
    return LinearSegmentedColormap.from_list('custom_colormap', colorList)



def plotProbabilities(probAA, totalSubs, enzymeName):
    # Create heatmap
    cMapCustom = createCustomColorMap(colorType='Probability')


    # Plot the heatmap with numbers centered inside the squares
    fig, ax = plt.subplots(figsize=figSize)
    heatmap = sns.heatmap(probAA, annot=True, fmt='.3f', cmap=cMapCustom,
                          cbar=True, linewidths=lineThickness-1,
                          linecolor='black', square=False, center=None,
                          annot_kws={'fontweight': 'bold'})
    ax.set_xlabel('Substrate Position', fontsize=labelSizeAxis)
    ax.set_ylabel('Residue', fontsize=labelSizeAxis)
    ax.set_title(f'\n{enzymeName}\nProbability Distribution',
                 fontsize=labelSizeTitle, fontweight='bold')
    plt.subplots_adjust(top=figBorders[0], bottom=figBorders[1],
                        left=figBorders[2]+0.01, right=figBorders[3]+0.085)

    # Set the thickness of the figure border
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(lineThickness)

    # Set tick parameters
    ax.tick_params(axis='both', which='major', length=tickLength,
                   labelsize=labelSizeTicks, width=lineThickness)
    ax.tick_params(axis='y', labelrotation=0)

    # Set x-ticks
    xTicks = np.arange(len(probAA.columns)) + 0.5
    ax.set_xticks(xTicks)
    ax.set_xticklabels(probAA.columns)

    # Set y-ticks
    yTicks = np.arange(len(probAA.index)) + 0.5
    ax.set_yticks(yTicks)
    ax.set_yticklabels(probAA.index)


    for _, spine in ax.spines.items():
        spine.set_visible(True)

    # Modify the colorbar
    cbar = heatmap.collections[0].colorbar
    cbar.ax.tick_params(axis='y', which='major', labelsize=labelSizeTicks,
                        length=tickLength, width=lineThickness)
    cbar.outline.set_linewidth(lineThickness)
    cbar.outline.set_edgecolor('black')

    # Convert figure to base64
    imgStream = io.BytesIO()
    fig.savefig(imgStream, format='png')
    imgStream.seek(0)
    figProb = base64.b64encode(imgStream.read()).decode('utf-8')

    return figProb



def plotEntropy(probAA, AA, enzymeName):
    # Evaluate: Entropy
    entropy = pd.DataFrame(0.0, index=probAA.columns, columns=['ΔS'])
    entropyMax = np.log2(len(AA))
    for indexColumn in probAA.columns:
        S = 0  # Reset entropy total for a new position
        for indexRow, probRatio in probAA.iterrows():
            prob = probRatio[indexColumn]
            if prob == 0:
                continue
            else:
                S += -prob * np.log2(prob)
        entropy.loc[indexColumn, 'ΔS'] = entropyMax - S


    # Map entropy values to colors using the colormap
    colors = [(0, 'navy'),
              (0.3 / entropyMax, 'navy'),
              (0.7 / entropyMax, 'dodgerblue'),
              (0.97 / entropyMax, 'white'),
              (0.98 / entropyMax, 'white'),
              (1.0 / entropyMax, 'white'),
              (1.65 / entropyMax, 'red'),
              (3 / entropyMax, 'firebrick'),
              (1, 'darkred')]
    colorBar = LinearSegmentedColormap.from_list('custom_colormap', colors)

    # Map entropy values to colors using the colormap
    normalize = Normalize(vmin=0, vmax=entropyMax)  # Normalize the entropy values
    cMap = [colorBar(normalize(value)) for value in entropy['ΔS'].astype(float)]

    # Plotting the entropy values as a bar graph
    fig, ax = plt.subplots(figsize=(figSize[0], 7))
    plt.bar(entropy.index, entropy['ΔS'], color=cMap,
            edgecolor='black', linewidth=lineThickness, width=0.8)
    plt.xlabel('Substrate Position', fontsize=labelSizeAxis)
    ax.set_ylabel('ΔS', fontsize=labelSizeAxis, rotation=0, labelpad=15)
    plt.title(f'{enzymeName}', fontsize=labelSizeTitle, fontweight='bold')
    plt.subplots_adjust(top=figBorders[0], bottom=figBorders[1]+0.02,
                        left=figBorders[2]+0.015, right=figBorders[3]-0.035)

    # Set tick parameters
    ax.tick_params(axis='both', which='major', length=tickLength,
                   labelsize=labelSizeTicks)

    # Set x-ticks
    xTicks = np.arange(0, len(entropy.iloc[:, 0]), 1)
    ax.set_xticks(xTicks)
    ax.set_xticklabels(entropy.index, rotation=0, ha='center')
    for tick in ax.xaxis.get_major_ticks():
        tick.tick1line.set_markeredgewidth(lineThickness)  # Set tick width

    # Set y-ticks
    yMax = int(np.ceil(entropyMax))
    yTicks = range(0, yMax+1)
    yTickLabels = [f'{tick:.0f}' if tick != 0 else f'{int(tick)}' for tick in yTicks]
    ax.set_yticks(yTicks)
    ax.set_yticklabels(yTickLabels)
    for tick in ax.yaxis.get_major_ticks():
        tick.tick1line.set_markeredgewidth(lineThickness)  # Set tick width

    # Set edge thickness
    for spine in ax.spines.values():
        spine.set_linewidth(lineThickness)

    # Set axis limits
    ax.set_ylim(0, yMax)

    # Add color bar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=normalize, cmap=colorBar), cax=cax)
    cbar.ax.tick_params(axis='y', which='major', labelsize=labelSizeTicks,
                        length=tickLength, width=lineThickness)
    for tick in cbar.ax.yaxis.get_major_ticks():
        tick.tick1line.set_markeredgewidth(lineThickness)  # Set tick width
    cbar.outline.set_linewidth(lineThickness)

    # Convert figure to base64
    imgStream = io.BytesIO()
    fig.savefig(imgStream, format='png')
    imgStream.seek(0)
    figEntropy = base64.b64encode(imgStream.read()).decode('utf-8')

    return entropy, entropyMax, figEntropy



def plotWeblogo(probAA, entropy, entropyMax, N, enzymeName):
    # AA Types: Aliphatic, Acidic, Basic, Hydroxyl, Amide, Aromatic, Sulfur
    letterColors = ['darkgreen','firebrick','deepskyblue','pink','indigo','black','gold']
    
    # Set local parameters
    bigLettersOnTop = False
    if bigLettersOnTop:
        stackOrder = 'big_on_top'
    else:
        stackOrder = 'small_on_top'
    colors = {
        'A': letterColors[0],
        'R': letterColors[2],
        'N': letterColors[4],
        'D': letterColors[1],
        'C': letterColors[6],
        'E': letterColors[1],
        'Q': letterColors[4],
        'G': letterColors[0],
        'H': letterColors[2],
        'I': letterColors[0],
        'L': letterColors[0],
        'K': letterColors[2],
        'M': letterColors[6],
        'F': letterColors[5],
        'P': letterColors[0],
        'S': letterColors[3],
        'T': letterColors[3],
        'W': letterColors[5],
        'Y': letterColors[5],
        'V': letterColors[0]
    }
    
    # Calculate: Letter heights
    heights = probAA.copy()
    for pos in heights.columns:
        heights.loc[:, pos] = heights.loc[:, pos] * entropy.loc[pos, 'ΔS']
    

    # Rename column headers
    colHeights = heights.columns
    heights.columns = range(len(heights.columns))

    # Set -inf to zero
    if heights.isin([np.inf, -np.inf]).any().any():
        heights.replace([np.inf, -np.inf], 0, inplace=True)

    # Plot the sequence motif
    fig, ax = plt.subplots(figsize=figSize)
    motif = logomaker.Logo(heights.transpose(), ax=ax, color_scheme=colors,
                           width=0.95, stack_order=stackOrder)
    motif.ax.set_title(f'\n{enzymeName}\nTotal Substrates = {int(N):,}\n',
                       fontsize=labelSizeTitle, fontweight='bold')

    # Set tick parameters
    ax.tick_params(axis='both', which='major', length=tickLength,
                   labelsize=labelSizeTicks)

    # Set borders
    motif.style_spines(visible=False)
    motif.style_spines(spines=['left', 'bottom'], visible=True)
    for spine in motif.ax.spines.values():
        spine.set_linewidth(lineThickness)

    # Set x-ticks
    motif.ax.set_xticks([pos for pos in range(len(colHeights))])
    motif.ax.set_xticklabels(colHeights, fontsize=labelSizeTicks,
                             rotation=0, ha='center')
    for tick in motif.ax.xaxis.get_major_ticks():
        tick.tick1line.set_markeredgewidth(lineThickness) # Set tick width

    # Set y-ticks
    yMax = int(np.ceil(entropyMax))
    yTicks = range(0, yMax+1)
    yTickLabels = [f'{tick:.0f}' if tick != 0 else f'{int(tick)}' for tick in yTicks]
    ax.set_yticks(yTicks)
    ax.set_yticklabels(yTickLabels)
    for tick in ax.yaxis.get_major_ticks():
        tick.tick1line.set_markeredgewidth(lineThickness) # Set tick width
    for spine in ax.spines.values():
        spine.set_linewidth(lineThickness) # Set edge thickness
    ax.set_ylim(0, yMax) # Set axis limits
    plt.subplots_adjust(top=figBorders[0], bottom=figBorders[1],
                        left=figBorders[2], right=figBorders[3])

    # Label the axes
    motif.ax.set_xlabel('Position', fontsize=labelSizeAxis)
    motif.ax.set_ylabel('Bits', fontsize=labelSizeAxis)

    # Evaluate dataset for fixed residues
    spacer = np.diff(motif.ax.get_xticks())  # Find the space between each tick
    spacer = spacer[0] / 2

    # Use the spacer to set a grey background to fixed residues
    for index, pos in enumerate(colHeights):
        if 1.0 in probAA.loc[:, pos]:
            # Plot grey boxes on each side of the xtick
            motif.ax.axvspan(index - spacer, index + spacer,
                             facecolor='darkgrey', alpha=0.2)

    # Convert figure to base64
    imgStream = io.BytesIO()
    fig.savefig(imgStream, format='png')
    imgStream.seek(0)
    figLogo = base64.b64encode(imgStream.read()).decode('utf-8')

    return figLogo


def binSubstrates(substrates, entropy, entropyMin, NSelect, enzymeName):
    # Identify the locations of specificity
    posSpecific = []
    entropyMin = float(entropyMin)
    for pos in entropy.index:
        if entropy.loc[pos, 'ΔS'] >= entropyMin:
            posSpecific.append(pos)

    # Bin substrates
    indexStart = entropy.index.get_loc(posSpecific[0])
    indexEnd = entropy.index.get_loc(posSpecific[-1]) + 1
    binnedSubs = {}
    countTotalSubs = 0
    countUniqueSubs = 0
    for substrate, count in substrates.items():
        countTotalSubs += count
        sub = substrate[indexStart:indexEnd]
        if sub in binnedSubs.keys():
            binnedSubs[sub] += count
        else:
            binnedSubs[sub] = count
            countUniqueSubs += 1

    # Sort the dictionary by counts from highest to lowest
    binnedSubs = dict(sorted(binnedSubs.items(), key=lambda item: item[1], reverse=True))

    # Plot: Binned substrates
    NBinSubs, figBinCounts, figBinProb, figWords = plotBinnedSubstrates(
        binnedSubs, countTotalSubs, NSelect, enzymeName)


    return NBinSubs, figBinCounts, figBinProb, figWords



def plotBinnedSubstrates(binnedSubs, N, NSelect, enzymeName):
    def plotBarGraph(xValues, yValues, yMax, yLabel, title):
        barColor = '#23FF55'
        barWidth = 0.75
        yMin = 0

        # Plot the data
        fig, ax = plt.subplots(figsize=figSize)
        bars = plt.bar(xValues, yValues, color=barColor, width=barWidth)
        plt.ylabel(yLabel, fontsize=labelSizeAxis)
        plt.title(title, fontsize=labelSizeTitle, fontweight='bold')
        plt.axhline(y=0, color='black', linewidth=lineThickness)
        plt.ylim(yMin, yMax)

        # Integer y-ticks if appropriate
        if yLabel == 'Counts' and yMax < 100:
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        # Set edge color
        for bar in bars:
            bar.set_edgecolor('black')

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=tickLength,
                       labelsize=labelSizeTicks, width=lineThickness)
        plt.xticks(rotation=90, ha='center')

        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(lineThickness)

        # Convert figure to base64
        imgStream = io.BytesIO()
        fig.savefig(imgStream, format='png')
        imgStream.seek(0)
        figBarGraph = base64.b64encode(imgStream.read()).decode('utf-8')

        return figBarGraph

    # Convert to integer
    NSelect = int(NSelect)

    # Evaluate: Counts
    NBinSubs = 0
    xCount, yCount, xProb, yProb = [], [], [], []
    for substrate, count in binnedSubs.items():
        NBinSubs += 1
        xCount.append(str(substrate))
        yCount.append(count)
        xProb.append(str(substrate))
        yProb.append(count / N)
        if NBinSubs == NSelect:
            break
    maxCounts, maxProb = np.ceil(max(yCount)), max(yProb)


    # Evaluate: Counts
    magnitude = np.floor(np.log10(maxCounts))
    unit = 10**(magnitude-1)
    if magnitude == 0:
        yMaxCount = maxCounts + 1
    else:
        yMaxCount = np.ceil(maxCounts / unit) * unit
    if yMaxCount < max(yCount):
        increaseValue = unit / 2
        while yMaxCount < max(yCount):
            yMaxCount += increaseValue

    # Evaluate: Prob
    magnitude = np.floor(np.log10(maxProb))
    adjustedMax = maxProb * 10**abs(magnitude)
    yMaxProb = np.ceil(adjustedMax) * 10**magnitude
    if yMaxProb == maxProb:
        yMaxProb += 10**magnitude
    adjVal = 5 * 10 ** (magnitude - 1)
    yMaxAdjusted = yMaxProb - adjVal
    if yMaxAdjusted > maxProb:
        yMaxProb = yMaxAdjusted


    # Make: Figures
    figBinCounts = plotBarGraph(xCount, yCount, yMaxCount, 'Counts', enzymeName)
    figBinProb = plotBarGraph(xProb, yProb, yMaxProb, 'Probability', enzymeName)

    # Evaluate: Word cloud
    figWords = plotWordCloud(binnedSubs)

    return NBinSubs, figBinCounts, figBinProb, figWords



def plotWordCloud(binnedSubs):
    cmap = createCustomColorMap(colorType='Word Cloud')

    # Create word cloud
    wordcloud = (WordCloud(
        width=950,
        height=800,
        background_color='white',
        min_font_size=10,  # Minimum font size
        max_font_size=100,  # Maximum font size
        scale=5,  # Increase scale for larger words
        colormap=cmap  # cool, hsv, plasma, _
    ).generate_from_frequencies(binnedSubs))


    # Create a figure
    fig, ax = plt.subplots(figsize=figSize)
    ax.imshow(wordcloud, interpolation='bilinear')
    ax.axis('off') # Turn off axis

    # Convert figure to base64
    imgStream = io.BytesIO()
    fig.savefig(imgStream, format='png', bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)  # Close figure to free memory
    imgStream.seek(0)
    figWords = base64.b64encode(imgStream.read()).decode('utf-8')

    return figWords
