import base64
import io
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd


def calcEntropy(probAA, AA):
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

    print(f'Positional Entropy:\n{entropy}\n\nMax Entropy: '
          f'{entropyMax.round(4)}\n')
    
    return entropy, entropyMax


def plotEntropy(prob, AA, enzymeName):
    entropy, entropyMax = calcEntropy(prob, AA)
    
    # Figure parameters
    titleSize = 18
    labelSizeTitle = 18
    labelSizeAxis = 16
    labelSizeTicks = 13
    lineThickness = 1.5
    tickLength = 4

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
    cMap = [colorBar(normalize(value)) for value in entropy['ΔEntropy'].astype(float)]

    # Plotting the entropy values as a bar graph
    fig, ax = plt.subplots(figsize=(9.5, 5))
    plt.bar(entropy.index, entropy['ΔEntropy'], color=cMap,
            edgecolor='black', linewidth=lineThickness, width=0.8)
    plt.xlabel('Substrate Position', fontsize=labelSizeAxis)
    plt.ylabel('Positional Entropy (ΔS)', fontsize=labelSizeAxis)
    plt.title(f'{enzymeName}', fontsize=titleSize, fontweight='bold')

    # Set borders
    plt.subplots_adjust(top=0.923, bottom=0.124, left=0.096, right=0.943)
    # plt.subplots_adjust(top=0.93, bottom=0.11, left=0.096, right=0.943)

    # Set tick parameters
    ax.tick_params(axis='both', which='major', length=tickLength,
                   labelsize=labelSizeTicks)

    # Set xticks
    xTicks = np.arange(0, len(entropy.iloc[:, 0]), 1)
    ax.set_xticks(xTicks)
    ax.set_xticklabels(entropy.index, rotation=0, ha='center')
    for tick in ax.xaxis.get_major_ticks():
        tick.tick1line.set_markeredgewidth(lineThickness)  # Set tick width

    # Set yticks
    yMax = np.ceil(entropyMax)
    print(f'Y Max: {yMax}')
    yTicks = range(0, yMax)
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

    # Convert the plot to a base64-encoded string
    imgStream = io.BytesIO()
    fig.savefig(imgStream, format='png')
    imgStream.seek(0)
    imgBase64 = base64.b64encode(imgStream.read()).decode('utf-8')

    return entropy, imgBase64
