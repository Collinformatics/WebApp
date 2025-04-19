import base64
import io



class WebApp:
    def __init__(self):
        self.buttonState = False
        self.message = ''



    def pressButton(self):
        if self.buttonState:
            self.buttonState = False
        else:
            self.buttonState = True
        print(f'Button State: {self.buttonState}')



    def makeFigure(data):
        import matplotlib
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd

        # Set matplotlib backend as non-interactive
        matplotlib.use('Agg')

        # Figure parameters
        labelSizeTitle = 18
        labelSizeAxis = 16
        labelSizeTicks = 13
        lineThickness = 1.5
        tickLength = 4
        figSize = (12, 9)
        figBorders = [0.882, 0.075, 0.05, 0.98]
        dpi = 300


        # Plot the heatmap with numbers centered inside the squares
        fig, ax = plt.subplots(figsize=figSize, dpi=dpi)
        # =============== Define figure here ===============
        ax.set_xlabel('X Label', fontsize=labelSizeAxis)
        ax.set_ylabel('Y Label', fontsize=labelSizeAxis)
        plt.title('Title', fontsize=labelSizeTitle, fontweight='bold')
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
        xTicks = np.arange(len(data.columns)) + 0.5
        ax.set_xticks(xTicks)
        ax.set_xticklabels(data.columns)

        # Set y-ticks
        yTicks = np.arange(len(data.index)) + 0.5
        ax.set_yticks(yTicks)
        ax.set_yticklabels(data.index)

        # Set spines
        for _, spine in ax.spines.items():
            spine.set_visible(True)

        # Convert figure to base64
        imgStream = io.BytesIO()
        fig.savefig(imgStream, format='png')
        imgStream.seek(0)
        figure = base64.b64encode(imgStream.read()).decode('utf-8')

        return figure
