import os
from matplotlib import pyplot as plot
from matplotlib import patches as mplPatches
import numpy as np
import sklearn


def readKmerFiles(directory):
    keys = ("start", "end", "mean", "std", "min", "max", "duration") #the column labels for kmer segment data
    kmerDataSeparated = list()
    kmerDataBulk = list()

    for f,filename in enumerate(os.listdir(directory)):
        kmerDataSeparated.append(list())

        with open(directory+'/'+filename,'r') as file:
            file.readline()

            for line in file:
                kmer = list(map(float,line.split('\t')))
                kmerDataBulk.append(kmer)
                kmerDataSeparated[f].append(kmer)

    return kmerDataSeparated,kmerDataBulk


def readStandardKmerMeans(filename): # file with 2 columns: kmer sequence, and expected mean. Includes column headers

    standardKmerMeans = dict()

    with open(filename,'r') as stdMeansFile:
        stdMeansFile.readline()

        for line in stdMeansFile:
            data = line.strip().split()

            try:
                standardKmerMeans[data[0]] = float(data[1])
            except:
                print("WARNING: duplicate kmer found in standard means reference list: %s"%data[0])

    return standardKmerMeans


def calculateExpectedSignal(standardKmerMeans,sequence):
    k = 6
    signal = list()

    for i in range(0,len(sequence)-k+1):
        kmer = sequence[i:i+k]

        signal.append(standardKmerMeans[kmer])

    return signal


def generateHistogram(xmin,xmax,interval,y):
    # Establish axis with defined limits
    panel = plot.axes()

    bins = np.arange(xmin,xmax,interval)

    frequencies, binsPlaceholder = np.histogram(y, bins)

    for i,frequency in enumerate(frequencies):
        left = bins[i]
        bottom = 0
        width = bins[i+1]-bins[i]
        height = frequency

        rectangle = mplPatches.Rectangle((left, bottom), width, height,
                                          linewidth=0.1,
                                          facecolor=(0.5, 0.5, 0.5),
                                          edgecolor=(0, 0, 0))
        panel.add_patch(rectangle)

    print(frequencies)

    panel.set_xlim([0, xmax])
    panel.set_ylim([0, 1.1*max(frequencies)])

    # Turn off top/right ticks
    panel.tick_params(axis='both', which='both',
                      bottom='on', labelbottom='on',
                      left='on', labelleft='on',
                      right='off', labelright='off',
                      top='off', labeltop='off')

def extractKmerMeans(kmerDataSeparated):
    signalsSeparated = list()
    for sequence in kmerDataSeparated:
        signalsSeparated.append(np.array(sequence)[:,2])

    return signalsSeparated


def histKmers(kmerMeans):
    plot.figure()

    generateHistogram(0,240,2.5,kmerMeans)

    plot.show()


def plotSegmentedSignals(signalsSeparated):
    nplots = len(signalsSeparated)

    f, axes = plot.subplots(nplots, sharex=True, sharey=True)

    for s,signal in enumerate(signalsSeparated):
        print(s)
        print(len(signal))

        colors = ["red","blue","orange","purple","green","yellow","brown","black","gray"]

        if s >= len(signalsSeparated)-2:
            color = colors[-(s+1)]
        else:
            color = "black"

        # plot.figure()
        # panel = plot.axes()
        # panel.plot(means)

        n = len(signal)
        # length = 100 / float(n)
        length = 5

        for i,kmerMean in enumerate(signal):
            x0 = (i+1)*length
            x1 = x0+length
            y = kmerMean

            if i >0:
                xprev = i*length + length
                yprev = signal[i-1]
                axes[s].plot([xprev,x0],[yprev,y],color=color)

            # print i,n,length,x0,x1,y

            axes[s].plot([x0,x1],[y,y],color=color)
            axes[s].text(x1, y, "%.3f"%y, ha="right", va="top", color="red", fontsize="5")
            axes[s].set_ylabel("%d"%s)
            # axes[s].set_ylim([0,150])
            # axes[s].set_yticks(np.linspace(0,150,4))

            axes[s].tick_params(axis='both', which='both',
                          bottom='off', labelbottom='off',
                          left='on', labelleft='on',
                          right='off', labelright='off',
                          top='off', labeltop='off')

    axes[0].set_title("Observed vs Expected Signal")

    # f.savefig("barcodeComparison2.png",dpi=600)
    plot.show()




#
directory = "allLeadKmers"
kmerDataSeparated,kmerDataBulk = readKmerFiles(directory)

# # standardKmerMeans = readStandardKmerMeans("kmerMeans")
# #
# # barcodeSignal1 = calculateExpectedSignal(standardKmerMeans,"TTGATTGAACCCTTGATTGAACCCTTGATTGAACCCTTGATTGAACCCAAAAAA")
# # print(barcodeSignal1)
# #
# # barcodeSignal2 = calculateExpectedSignal(standardKmerMeans,"GGGTTCAATCAAGGGTTCAATCAAGGGTTCAATCAAGGGTTCAATCAAT")
# # print(barcodeSignal2)
# #
# # barcodeSignalBoth = barcodeSignal2 + barcodeSignal1
# #
# # signalsBulk = np.array(kmerDataBulk)[:,2]
#
# signals = extractKmerMeans(kmerDataSeparated)
#
# plotSegmentedSignals(signals)
#
# # signalsSeparated.append(np.array(barcodeSignal1))
# # signalsSeparated.append(np.array(barcodeSignal2))
# # # signalsSeparated.append(np.array(barcodeSignalBoth))
# #
# # np.set_printoptions(suppress=True)