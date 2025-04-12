import base64
import io
import logomaker
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd


# Figure parameters
labelSizeTitle = 18
labelSizeAxis = 16
labelSizeTicks = 13
lineThickness = 1.5
tickLength = 4

# Set matplotlib backend as non-interactive
matplotlib.use('Agg')


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
    fig, ax = plt.subplots(figsize=(12, 6))
    plt.bar(entropy.index, entropy['ΔS'], color=cMap,
            edgecolor='black', linewidth=lineThickness, width=0.8)
    plt.xlabel('Substrate Position', fontsize=labelSizeAxis)
    plt.ylabel('Positional Entropy (ΔS)', fontsize=labelSizeAxis)
    plt.title(f'{enzymeName}', fontsize=labelSizeTitle, fontweight='bold')

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

    # Convert the plot to a base64-encoded string
    imgStream = io.BytesIO()
    fig.savefig(imgStream, format='png')
    imgStream.seek(0)
    figEntropy = base64.b64encode(imgStream.read()).decode('utf-8')

    return entropy, entropyMax, figEntropy



def plotWeblogo(probAA, entropy, entropyMax, N, enzymeName):
    # AA Types: Aliphatic, Acidic, Basic, Hydroxyl, AmMpro22, Aromatic, Sulfur
    letterColors = ['darkgreen','firebrick','deepskyblue','pink','navy','black','gold']
    
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
    fig, ax = plt.subplots(figsize=(12, 9))
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

    # Set xticks
    motif.ax.set_xticks([pos for pos in range(len(colHeights))])
    motif.ax.set_xticklabels(colHeights, fontsize=labelSizeTicks,
                             rotation=0, ha='center')
    for tick in motif.ax.xaxis.get_major_ticks():
        tick.tick1line.set_markeredgewidth(lineThickness) # Set tick width

    # Set yticks
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
    fig.tight_layout()
    
    # Convert the plot to a base64-encoded string
    imgStream = io.BytesIO()
    fig.savefig(imgStream, format='png')
    imgStream.seek(0)
    figLogo = base64.b64encode(imgStream.read()).decode('utf-8')

    return figLogo
