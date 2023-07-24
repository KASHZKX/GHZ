from cmath import log10
import numpy as np
import math
import os
import matplotlib.pyplot as plt
import matplotlib.transforms
# import latex
import matplotlib
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

directory_path = "encode_math/"

class ChartGenerator:
    # data檔名 Y軸名稱 X軸名稱 Y軸要除多少(10的多少次方) Y軸起始座標 Y軸終止座標 Y軸座標間的間隔
    def __init__(self, dataName, Ylabel, Xlabel):
        filename = directory_path + dataName
        print("start generate", filename)

        if not os.path.exists(filename):
            print("file doesn't exist")
            return
        
        with open(filename, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            
        
        
        # Ydiv, Ystart, Yend, Yinterval
        Ypow = 0
        Xpow = 0

        # color = [
        #     "#000000",
        #     "#006400",
        #     "#FF1493",
        #     "#7FFF00",   
        #     "#900321",
        # ]
        color = [
            "#AA3311",   
            "#FF0000",
            "#FF6600",
            "#901231",
        ]
        # matplotlib.rcParams['text.usetex'] = True

        fontsize = 30
        Xlabel_fontsize = fontsize
        Ylabel_fontsize = fontsize
        Xticks_fontsize = fontsize
        Yticks_fontsize = fontsize
            
        # matplotlib.rcParams['text.usetex'] = True
        # fig, ax = plt.subplots(figsize=(8, 6), dpi=600) 
        
        andy_theme = {
        # "axes.grid": True,
        # "grid.linestyle": "--",
        # "legend.framealpha": 1,
        # "legend.facecolor": "white",
        # "legend.shadow": True,
        # "legend.fontsize": 14,
        # "legend.title_fontsize": 16,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "axes.labelsize": 20,
        "axes.titlesize": 20,
        "font.family": "Times New Roman",
        "mathtext.it": "Times New Roman:italic",
        # "mathtext.default": "regular",
        "mathtext.fontset": "custom"
        # "mathtext.fontset": "custom"
        # "figure.autolayout": True
        # "text.usetex": True,
        # "figure.dpi": 100,
        }
        
        matplotlib.rcParams.update(andy_theme)
        fig, ax1 = plt.subplots(figsize = (7, 6), dpi = 600)
        # ax1.spines['top'].set_linewidth(1.5)
        # ax1.spines['right'].set_linewidth(1.5)
        # ax1.spines['bottom'].set_linewidth(1.5)
        # ax1.spines['left'].set_linewidth(1.5)
        ax1.tick_params(direction = "in")
        ax1.tick_params(bottom = True, top = True, left = True, right = True)
        ax1.tick_params(pad = 10)


        ##data start##
        x = []
        x_data = []
        _y = []
        numOfData = 0
        counter = 0
        for line in lines:
            line = line.replace('\n','')
            data = line.split(' ')
            numOfline = len(data)
            numOfData += 1
            for i in range(numOfline):
                if i == 0:
                    x.append(data[i])
                    x_data.append(counter)
                elif i != numOfline-1:
                    _y.append(data[i])
            # print(_y)
            counter += 1
        numOfAlgo = len(_y) // numOfData

        y = [[] for _ in range(numOfAlgo)]
        for i in range(numOfData * numOfAlgo):
            y[i % numOfAlgo].append(_y[i])

        # print(x)
        # print(y)

        maxData = 0
        minData = math.inf

        for i in range(-10, -1, 1):
            if float(x[numOfData - 1]) <= 10 ** i:
                Xpow = (i - 2)
        
        #Ypow = -2
        Ydiv = float(10 ** Ypow)
        Xdiv = float(10 ** Xpow)
        

        for i in range(numOfData):
            x[i] = int(float(x[i]) / Xdiv)

        for i in range(numOfAlgo):
            for j in range(numOfData):
                y[i][j] = float(y[i][j]) / Ydiv
                maxData = max(maxData, y[i][j])
                minData = min(minData, y[i][j])

        Yend = 0.74
        Ystart = 0.5
        Yinterval = 0.05

        # if maxData > 1.1:
        #     Yinterval = int(math.ceil(Yinterval))
        #     Yend = int(Yend)
        # else:
        #     Yend = 1
        #     Ystart = 20
        #     Yinterval = 0.2

        marker = ['o', 'x', 's', 'v']
        line_style = ["-.", "-", ":", "--"]
        AlgoName = ["Fidelity", "Waiting Time"]
        # color.reverse()
        # marker.reverse()
        #for i in range(numOfAlgo):
        line_plot = []
        line_plot.append(ax1.plot(x_data, y[0], color = color[0], lw = 2.5, linestyle = line_style[0], marker = marker[0], markersize = 15, markerfacecolor = "none", markeredgewidth = 2.5, zorder=0, label=AlgoName[0]))
        # plt.show()

        plt.xticks(fontsize = Xticks_fontsize)
        plt.yticks(fontsize = Yticks_fontsize)
        
        leg = plt.legend(
            ["Fidelity"],
            loc = 10,
            bbox_to_anchor = (0.15, 1.25),
            prop = {"size": fontsize, "family": "Times New Roman"},
            frameon = "False",
            labelspacing = 0.2,
            handletextpad = 0.2,
            handlelength = 1,
            columnspacing = 0.2,
            ncol = 2,
            facecolor = "None",
        )

        leg.get_frame().set_linewidth(0.0)
        
        
        #Ylabel += self.genMultiName(Ypow)
        Xlabel += self.genMultiName(Xpow)
        plt.subplots_adjust(top = 0.75)
        plt.subplots_adjust(left = 0.225)
        plt.subplots_adjust(right = 0.825)
        plt.subplots_adjust(bottom = 0.25)

        plt.yticks(np.arange(Ystart, Yend + Yinterval, step = Yinterval), fontsize = Yticks_fontsize)
        plt.xticks(x_data, x)
        plt.ylabel(Ylabel, fontsize = Ylabel_fontsize, labelpad = 15)
        plt.xlabel(Xlabel, fontsize = Xlabel_fontsize, labelpad = 10)
        ax1.yaxis.set_label_coords(-0.225, 0.5)
        ax1.xaxis.set_label_coords(0.45, -0.27)
        # plt.show()
        # plt.tight_layout()

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        
        ax2.set_ylabel('Waiting Time', fontsize = Ylabel_fontsize, labelpad = 15, rotation=-90)  # we already handled the x-label with ax1
        line_plot.append(ax2.plot(x_data, y[1], color = color[1], lw = 2.5, linestyle = line_style[1], marker = marker[1], markersize = 15, markerfacecolor = "none", markeredgewidth = 2.5, zorder=1, label=AlgoName[1]))
        ax2.yaxis.set_label_coords(1.25, 0.5)
        plt.yticks(fontsize=Ylabel_fontsize)

        leg = plt.legend(
            ["Waiting Time"],
            loc = 10,
            bbox_to_anchor = (0.75, 1.25),
            prop = {"size": fontsize, "family": "Times New Roman"},
            frameon = "False",
            labelspacing = 0.2,
            handletextpad = 0.2,
            handlelength = 1,
            columnspacing = 0.2,
            ncol = 2,
            facecolor = "None",
        )

        leg.get_frame().set_linewidth(0.0)

        pdfName = dataName[0:-4]
        plt.savefig(directory_path + '{}.eps'.format(pdfName)) 
        plt.savefig(directory_path + '{}.jpg'.format(pdfName)) 
        # Xlabel = Xlabel.replace(' (%)','')
        # Xlabel = Xlabel.replace('# ','')
        # Ylabel = Ylabel.replace('# ','')
        plt.close()

    def genMultiName(self, multiple):
        if multiple == 0:
            return str()
        else:
            return "($" + "10" + "^{" + str(multiple) + "}" + "$)"

if __name__ == "__main__":
    # data檔名 Y軸名稱 X軸名稱 Y軸要除多少(10的多少次方) Y軸起始座標 Y軸終止座標 Y軸座標間的間隔
    # ChartGenerator("numOfnodes_waitingTime.txt", "need #round", "#Request of a round", 0, 0, 25, 5)
    Xlabels = ["min_fidelity", "new_request_cnt"]
    #Ylabels = ["encode_cnt", "unencode_cnt", "encode_ratio", "divide_cnt", "finished_throughputs", "encode_use_one_path_rate"]
    Ylabels = ["encode_ratio", "divide_cnt", "encode_use_one_path_rate", "use_channel_ratio", "use_memory_ratio"]
    
    # for Xlabel in Xlabels:
    #     for Ylabel in Ylabels:
    #         dataFileName = Xlabel + '_' + Ylabel + '.ans'
    #         ChartGenerator(dataFileName, Ylabel, Xlabel)


    Xlabel = "# Encoded Process"
    Ylabel = "Fidelity"
    dataFileName = "encode_math.ans"
    ChartGenerator(dataFileName, Ylabel, Xlabel)
