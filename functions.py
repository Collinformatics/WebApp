import base64
import io
import logomaker
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from wordcloud import WordCloud


# Figure parameters
labelSizeTitle = 18
labelSizeAxis = 16
labelSizeTicks = 13
lineThickness = 1.5
tickLength = 4
figSize = (12, 9)
figBorders = [0.882, 0.075, 0.05, 0.98]
dpi = 300

# Experimental parameters
AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
      'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# Set matplotlib backend as non-interactive
matplotlib.use('Agg')


def createCustomColorMap(colorType):
    if colorType == 'Probability':
        colors = ['#FFFFFF', '#ABFF9B', '#39FF14', '#2E9418', '#2E9418', '#005000']
    elif colorType == 'Word Cloud':
        colors = ['#9BFF9B', '#39FF14', '#2FB513', '#006000', 'black']
    else:
        colors = ['#008631', '#39E75F', '#CC5500', '#F79620', 'black']

    # Create colormap
    if len(colors) == 1:
        colorList = [(0, colors[0]), (1, colors[0])]
    else:
        colorList = [(i / (len(colors) - 1), color) for i, color in enumerate(colors)]

    return LinearSegmentedColormap.from_list('custom_colormap', colorList)



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
    N = 0
    countedAA = pd.DataFrame(0, index=AA, columns=pos)
    for sub, count in substrates.items():
        N += count
        for index, aa in enumerate(sub):
            countedAA.loc[aa, pos[index]] += count

    # Evaluate: AA probability
    probAA = countedAA / N

    # Figure: Entropy
    figProb = plotProbabilities(probAA, N, enzymeName)
    entropy, entropyMax, figEntropy = plotEntropy(probAA, AA, enzymeName)
    figLogo = plotWeblogo(probAA, entropy, entropyMax, N, enzymeName)
    NBinSubs, figBarCounts, figBarProb, figWords, figTrie = binSubstrates(
        substrates, entropy, entropyMin, NSelect, enzymeName)

    # Create dataset
    dataset = {}
    dataset['N'] = N
    dataset['NBinSubs'] = NBinSubs
    dataset['probability'] = figProb
    dataset['entropy'] = figEntropy
    dataset['pLogo'] = figLogo
    dataset['barCounts'] = figBarCounts
    dataset['barProb'] = figBarProb
    dataset['words'] = figWords
    dataset['trie'] = figTrie

    return dataset



def subsDefault():
    substrates = {
        'VVLQSGTK': 19076,
        'VILQSVGA': 18383,
        'LALQAACW': 17331,
        'LNLQGILD': 16201,
        'IYLQALMP': 16010,

        'TSLQARKS': 15449,
        'AFLQAHFT': 15175,
        'VTLQCTYS': 14330,
        'VLLQGKQL': 14201,
        'LVLQANPC': 14010,

        'IMLQGVIW': 13449,
        'AGLQASAH': 12175,
        'IHLQSENE': 11330,
        'MVLQGDVN': 10201,
        'YGLQCNEV': 10010,

        'VCMQCQVQ': 9449,
        'AELQSWHF': 7405,
        'VDMQCMWG': 7330,
        'LWMQCSII': 7013,
        'VVMQCCSM': 6276,

        'VLMQCCPR': 1991,
        'VVMQSGSM': 1754,
        'VVFQCHNR': 1576,
        'KRFQCKLR': 1234,
        'PEFQCCLQ': 900
        }

    return substrates



def plotProbabilities(probAA, N, enzymeName):
    # Create heatmap
    cMapCustom = createCustomColorMap(colorType='Probability')


    # Plot the heatmap with numbers centered inside the squares
    fig, ax = plt.subplots(figsize=figSize, dpi=dpi)
    heatmap = sns.heatmap(probAA, annot=True, fmt='.3f', cmap=cMapCustom,
                          cbar=True, linewidths=lineThickness-1,
                          linecolor='black', square=False, center=None,
                          annot_kws={'fontweight': 'bold'})
    ax.set_xlabel('Substrate Position', fontsize=labelSizeAxis)
    ax.set_ylabel('Residue', fontsize=labelSizeAxis)
    plt.title(f'\n{enzymeName}\nProbability Distribution\nN = {N:,}',
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
    letterColors = ['#00D500','#CE0000','#08FFFF','#FF0088','#662695','#151515','#FFEA00']

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
    fig, ax = plt.subplots(figsize=figSize, dpi=dpi)
    motif = logomaker.Logo(heights.transpose(), ax=ax, color_scheme=colors,
                           width=0.95, stack_order=stackOrder)
    plt.title(f'\n{enzymeName}\nN = {int(N):,}\n',
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

    # Plot: Binned substrates NBinSubs, figBinCounts, figBinProb, figWords, figTrie
    output = plotBinnedSubstrates(
        binnedSubs, countTotalSubs, NSelect, entropy, entropyMin, enzymeName)

    return output



def plotBinnedSubstrates(binnedSubs, N, NSelect, entropy, entropyMin, enzymeName):
    def plotBarGraph(xValues, yValues, yMax, yLabel, title):
        barColor = '#23FF55'
        barWidth = 0.75
        yMin = 0

        # Plot the data
        fig, ax = plt.subplots(figsize=figSize, dpi=dpi)
        bars = plt.bar(xValues, yValues, color=barColor, width=barWidth)
        plt.ylabel(yLabel, fontsize=labelSizeAxis)
        plt.title(f'{title}\nN = {N:,}', fontsize=labelSizeTitle, fontweight='bold')
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


    # Evaluate: Motifs
    motifs = {}
    motifsList, yCount, yProb = [], [], []
    for NBinSubs, (substrate, count) in enumerate(binnedSubs.items()):
        NBinSubs += 1
        motifs[substrate] = count
        motifsList.append(str(substrate))
        yCount.append(count)
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
    figBinCounts = plotBarGraph(motifsList, yCount, yMaxCount, 'Counts', enzymeName)
    figBinProb = plotBarGraph(motifsList, yProb, yMaxProb, 'Probability', enzymeName)

    # Evaluate: Word cloud
    figWords = plotWordCloud(motifs, N, enzymeName)

    # Evaluate: Suffix tree
    figTrie, NBinSubs = plotSuffixTree(motifs, N, NSelect, entropy, entropyMin, enzymeName)


    return NBinSubs, figBinCounts, figBinProb, figWords, figTrie



def plotWordCloud(motifs, N, enzymeName):
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
    ).generate_from_frequencies(motifs))


    # Create a figure
    fig, ax = plt.subplots(figsize=figSize, dpi=dpi)
    plt.title(f'\n{enzymeName}\nN = {N:,}\n',
              fontsize=labelSizeTitle-5, fontweight='bold')
    ax.imshow(wordcloud, interpolation='bilinear')
    ax.axis('off') # Turn off axis

    # Set edge thickness
    for spine in ax.spines.values():
        spine.set_linewidth(lineThickness)

    # Convert figure to base64
    imgStream = io.BytesIO()
    fig.savefig(imgStream, format='png', bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)  # Close figure to free memory
    imgStream.seek(0)
    figWords = base64.b64encode(imgStream.read()).decode('utf-8')

    return figWords



class TrieNode:
    # Create nodes
    def __init__(self):
        self.children = {}
        self.endOfWord = False



class Trie:
    def __init__(self):
        # Create the root as an empty node
        self.root = TrieNode()

    def insert(self, word: str) -> None:
        # Add words to the trie
        current =  self.root

        for AA in word:
            if AA not in current.children:
                current.children[AA] = TrieNode()
            current = current.children[AA]

    def search(self, word):
        # Search for a word in the trie
        current = self.root

        for AA in word:
            if AA not in current.children:
                return False
            current = current.children[AA]
        return True

    def startsWith(self, prefix):
        # Find the prefix in the trie
        current = self.root

        for AA in prefix:
            if AA not in current.children:
                return False
        return True



def evaluateSubtrees(trie, motifTrie):
    def subtreeTable(subtreeFreq):
        # Sort motifs by length
        sortedMotifs = sorted(subtreeFreq.keys(), key=len)

        # Organize motifs by their length and sort by frequency (highest first)
        motifGroups = {}
        for motif in sortedMotifs:
            length = len(motif)
            if length not in motifGroups:
                motifGroups[length] = []
            motifGroups[length].append((motif, subtreeFreq[motif]))

        # Sort motifs in each length group by frequency (descending)
        for length in motifGroups:
            motifGroups[length].sort(key=lambda x: x[1], reverse=True)

        # Convert motifs back to formatted strings
        for length in motifGroups:
            motifGroups[length] = [f"{motif}: {round(freq, 5)}"
                                   for motif, freq in motifGroups[length]]

        # Find the max number of motifs in any length group
        maxRows = max(len(motifs) for motifs in motifGroups.values())

        # Construct the table row by row
        tableData = []
        for i in range(maxRows):
            row = []
            for length in sorted(motifGroups.keys()):
                motifs = motifGroups[length]
                row.append(motifs[i] if i < len(
                    motifs) else "")  # Fill missing values with empty strings
            tableData.append(row)

        # Convert to DataFrame
        motifTable = pd.DataFrame(tableData,
                                  index=range(1, len(motifTrie.keys()) + 1),
                                  columns=[str(length)
                                           for length in sorted(motifGroups.keys())])

        return motifTable


    # Count: Total substrates
    motifsTotal = 0
    for motif, count in motifTrie.items():
        motifsTotal += count

    # Evaluate: Partial sequence counts
    subtreeCount = {}
    motifLength = len(next(iter(motifTrie)))
    for index in range(motifLength):
        for motif, count in motifTrie.items():
            subSeq = motif[0:index+1]
            if subSeq in subtreeCount.keys():
                subtreeCount[subSeq] += count
            else:
                subtreeCount[subSeq] = count

    # Evaluate: Partial sequence frequency
    subtreeFreq = {}
    for subSeq, count in subtreeCount.items():
        subtreeFreq[subSeq] = count / motifsTotal
    prevSeqLen = 1
    for subSeq, count in subtreeFreq.items():
        if len(subSeq) != prevSeqLen:
            prevSeqLen = len(subSeq)
    motifTable = subtreeTable(subtreeFreq)

    return motifTable



def plotSuffixTree(motifs, N, NSelect, entropy, entropyMin, enzymeName):
    # Figure Parameters: Suffix Tree
    inOffset = 100
    inNodeSizeMax = 1700
    inNodeSizeMin = 500
    inFontSize = 16
    inScaleX = 2
    inScaleY = 1
    inNodeCoordSpacer = 100 # Increase space between clusters


    # Evaluate: Specificity
    motifPos = entropy.copy()
    for indexPos, position in enumerate(entropy.index):
        if entropy.loc[position, 'ΔS'] < entropyMin:
            motifPos.drop(position, inplace=True)

    # Sort the frame
    motifPos = motifPos.sort_values(by='ΔS', ascending=False).copy()

    # Initialize Trie
    trie = Trie()

    # Find motif positions based on entropy threshold
    indexPos = []
    levelLabels = ['Level: Pos']
    for index in motifPos.index:
        posEntropy = motifPos.loc[index, 'ΔS']
        if posEntropy >= entropyMin:
            levelLabels.append(index)
            indexPos.append(int(index.replace('R', '')) - 1)

    # Adjust indices
    motifLen = len(next(iter(motifs.keys()))) - 1
    if max(indexPos) > motifLen:
        adjustment = max(indexPos) - motifLen
        indexPos = [x - adjustment for x in indexPos]


    def addMotif(motif, count):
        # Extract important AAs from the motif
        motif = ''.join(motif[index] for index in indexPos)

        # Add motif to the trie
        if motif in motifTrie.keys():
            motifTrie[motif] += count
        else:
            motifTrie[motif] = count
            trie.insert(motif)

    # Extract the motifs
    motifTrie = {}
    NUniqueTrieMotifs = 0 # Determine the number of unique sequences in the trie
    for motif, count in motifs.items():
        # Add the motif to the tree
        addMotif(motif, count)
        NUniqueTrieMotifs = len(motifTrie.keys())
        if NUniqueTrieMotifs >= NSelect:
            break
    motifTrie = dict(sorted(motifTrie.items(), key=lambda item: item[1],
                            reverse=True))

    # Calculate: RF
    motifTable = evaluateSubtrees(trie=trie, motifTrie=motifTrie)


    def computeNodeSizes(motifTable, inNodeSizeMax, inNodeSizeMin):
        # Initialize data structures
        nodes = [inNodeSizeMax]
        nodeSizesDF = pd.DataFrame('',
                                   index=motifTable.index,
                                   columns=motifTable.columns)

        # Calculate: Node size
        for col in motifTable.columns:
            # Extract all RF values
            rfValues = []
            for entry in motifTable[col].dropna():
                if ": " in entry:
                    _, rf = entry.split(": ")
                    rfValues.append(float(rf))

            if not rfValues:
                return nodeSizesDF  # return empty if no RFs

            rfMin = min(rfValues)
            rfMax = max(rfValues)

            for index, entry in enumerate(motifTable[col].dropna()):
                if ": " in entry:
                    motif, rf = entry.split(": ")
                    rf = float(rf)

                    # Normalize and scale size
                    if rfMax == rfMin:
                        scaled = 1  # avoid division by zero
                    else:
                        scaled = (rf - rfMin) / (rfMax - rfMin)

                    nodeSize = inNodeSizeMin + scaled * (inNodeSizeMax - inNodeSizeMin)

                    if len(motif) > 2:
                        motif = motif[-2:]

                    nodeSizesDF.loc[index + 1, col] = f'{motif}: {nodeSize:.2f}'
                    nodes.append(nodeSize)
        # # Inspect node sizes
        # print(f'Node Size:\n{nodeSizesDF}\n')

        return nodes


    def addNodesToGraph(node, graph, coordXMin, scaleY, offset, clusterSpacer):
        coords = {} # Stores node positions
        coordXMin = 0
        nodeCountLevel = {} # Track all nodes per level
        queue = [(node, None, '', 0)] # (node, parent, char, level)

        # Pass 1: Collect nodes for each level before positioning
        while queue:
            nodeCurrent, parent, char, level = queue.pop(0)
            nodeID = f"{char}-{level}-{id(nodeCurrent)}"

            if level not in nodeCountLevel:
                nodeCountLevel[level] = []
            nodeCountLevel[level].append((nodeCurrent, parent, char, nodeID))

            # Add children to the queue
            for child_char, nodeChild in nodeCurrent.children.items():
                queue.append((nodeChild, nodeID, child_char, level + 1))

        # Pass 2: Assign positions level by level
        for level, nodes in nodeCountLevel.items():
            nodeNumber = len(nodes)
            for i, (nodeCurrent, parent, char, nodeID) in enumerate(nodes):
                if parent is None:
                    coords[nodeID] = (0, 0) # Root node position
                else:
                    parentX, parentY = coords[parent]

                    # Calculate the x position
                    clusterSpacing = offset + (level / clusterSpacer)
                    if nodeNumber % 2 == 1:
                        # Odd number of nodes: Center one on parentX
                        posX = parentX + (i - nodeNumber // 2) * clusterSpacing
                    else:
                        # Even number of nodes: Spread symmetrically
                        posX = parentX + (i - (nodeNumber / 2 - 0.5)) * clusterSpacing

                    coords[nodeID] = (posX, parentY - scaleY)
                    if posX < coordXMin:
                        coordXMin = posX

                # Add node and edge to graph
                graph.add_node(nodeID, label=char)
                if parent is not None:
                    # graph.add_edge(parent, nodeID)
                    graph.add_edge(parent, nodeID, arrowstyle='->')

        return coords, coordXMin, nodeCountLevel.keys()



    # Evaluate: Node sizes
    nodes = computeNodeSizes(motifTable, inNodeSizeMax, inNodeSizeMin)

    # Build the graph
    graph = nx.DiGraph()
    coords, coordXMin, levels = addNodesToGraph(
        trie.root, graph, inScaleX, inScaleY, inOffset, inNodeCoordSpacer)

    # Get node labels
    labels = {node: data['label'] for node, data in graph.nodes(data=True)}

    # Plot the data
    fig, ax = plt.subplots(figsize=figSize, dpi=dpi)
    nx.draw(graph, coords, with_labels=True, labels=labels, node_size=nodes,
            node_color="#23FF55", font_size=inFontSize, font_weight="bold",
            edge_color="#101010", ax=ax, arrows=False) # Draw graph
    plt.title(f'\n{enzymeName}\nN = {N:,}\nUnique Sequences: {NUniqueTrieMotifs:,}',
              fontsize=16, fontweight='bold')
    plt.tight_layout()

    # Label each level
    xMin, xMax = ax.get_xlim()
    posX = xMin + (coordXMin - xMin)/2
    for index, level in enumerate(levels):
        if index == 0:
            label = levelLabels[index]
        else:
            xMin, xMax = ax.get_xlim()
            posX = xMin + (coordXMin - xMin - 66) / 2
            label = f'{index}: {levelLabels[index]}'
        ax.text(posX, -level, label, ha='right', va='center', fontsize=inFontSize,
                fontweight='bold', color='black', transform=ax.transData)

    # Convert figure to base64
    imgStream = io.BytesIO()
    fig.savefig(imgStream, format='png', bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)  # Close figure to free memory
    imgStream.seek(0)
    figTrie = base64.b64encode(imgStream.read()).decode('utf-8')

    return figTrie, NUniqueTrieMotifs
