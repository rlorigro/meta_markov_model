from Fast5Types import *
import numpy as np
from matplotlib import pyplot as plot
import copy
import os
import yahmm
from kmerAnalysis import *
import random
import sys


class KmerCalculator:
    '''
    Some tools for estimation of a signal from a sequence (if you happen to have a table of means for all kmers... :P)
    '''
    def __init__(self, file_stdKmerMeans):
        self.k = 6
        self.standardKmerMeans = self.readStandardKmerMeans(file_stdKmerMeans)


    def readStandardKmerMeans(self, file_stdKmerMeans):
        '''
        Read a file containing the list of all kmers and their expected signal means (2 columns, with headers)
        '''

        standardKmerMeans = dict()

        with open(file_stdKmerMeans, 'r') as file:
            file.readline()

            for line in file:
                data = line.strip().split()

                if len(data[0])!=self.k:
                    sys.exit("ERROR: kmer length not equal to parameter K")

                try:
                    standardKmerMeans[data[0]] = float(data[1])
                except:
                    print("WARNING: duplicate kmer found in standard means reference list: %s" % data[0])

        return standardKmerMeans


    def calculateExpectedSignal(self, sequence):
        '''
        Given a sequence, use standard kmer signal values to estimate its expected signal (over all kmers in sequence)
        '''

        signal = list()

        for i in range(0, len(sequence) - self.k + 1):
            kmer = sequence[i:i + self.k]

            signal.append(self.standardKmerMeans[kmer])

        return signal


class metaHMMBuilder:
    def __init__(self,file_standardKmerMeans):
        # self.metaModel = Model(name="meta")
        # self.metaStartNode = self.metaModel.start
        # self.metaEndNode = self.metaModel.end
        self.metaNodes = dict()
        self.metaTransitions = dict()
        self.HMMs = dict()
        self.HMM = None
        self.model = Model(name="meta")
        self.startNode = self.model.start
        self.endNode = self.model.end
        self.kmerCalc = KmerCalculator(file_standardKmerMeans)
        self.expectedMeans = None
        self.k = self.kmerCalc.k
        self.sequence = None
        self.sequenceIsString = None
        self.unknownBridgeStdev = 7.5

        #----Preset transition probabilities----#
        #For all non-terminal nodes. Node->end probabilities wil replace node->next probabilities.

        self.transitions = {"skip": {"next": 0.01, #unused
                                     "emit": 1.0},

                            "blip": {"emit": 0.5,
                                     "next": 0.5},

                            "drop": {"next": 0.7,
                                     "emit": 0.3},

                            "emit": {"self": 0.015,
                                     "next": 0.953,
                                     "prev": 0.0165,
                                     "skip": 0.01,
                                     "blip": 0.005,
                                     "drop": 0.0005},

                            "start": {"emit": 0.982,
                                      "skip": 0.015,
                                      "drop": 0.0025,
                                      "blip": 0.0005}}

        for key1 in self.transitions:
            totalProbability = float(sum([self.transitions[key1][key2] for key2 in self.transitions[key1]]))

            if totalProbability != 1.0:
                print("WARNING: Total outgoing transition probability for %s nodes != 1.0, p = %f"%(key1,totalProbability))

    def newPhase(self,skipNode=None,blipNode=None,dropNode=None,emitNode=None,startNode=None):
        '''
        Generate a container to be used for each phase in the linear HMM
        '''
        phase =  {"skip": skipNode,
                  "blip": blipNode,
                  "drop": dropNode,
                  "emit": emitNode}

        return phase

    def buildMetaHMM(self,nodes,transitions): #nodes is a dict of names: sequences, transitions is the dict of dict w/ probabilities
        self.metaNodes = nodes
        self.metaTransitions = transitions

        for nodeName in self.metaNodes:
            if nodeName != "start":

                data = nodes[nodeName]

                if type(data) == str or type(data) == list:
                    sequence = data

                    self.buildDraftHMM(nodeName,sequence)

                if type(data) == tuple:
                    print(len(data[0]))
                    print(len(data[1]))
                    self.buildDraftHMM(nodeName,data[0],stdevs=data[1])


        for nodeName in self.metaTransitions:
            if nodeName == "start":
                for nodeName2 in self.metaTransitions["start"]:
                    self.addStartTransitions(self.HMMs[nodeName2],self.metaTransitions["start"][nodeName2])
            else:
                for nodeName2 in self.metaTransitions[nodeName]:
                    if nodeName2 == "end":
                        self.addEndTransitions(self.HMMs[nodeName],self.metaTransitions[nodeName][nodeName2])
                    else:
                        self.bridgeHMMs(self.HMMs[nodeName],self.HMMs[nodeName2],self.metaTransitions[nodeName][nodeName2])


    class linearHMM:
        def __init__(self, name, HMM, expectedMeans, sequence, sequenceIsString):
            self.name = name
            self.HMM = HMM
            self.expectedMeans = expectedMeans
            self.sequence = sequence
            self.sequenceIsString = sequenceIsString

    def storeHMM(self,name):
        self.HMMs[name] = self.linearHMM(name,self.HMM,self.expectedMeans,self.sequence,self.sequenceIsString)

        self.clearHMM()

    def clearHMM(self):
        self.HMM = None
        self.expectedMeans = None
        self.sequence = None
        self.sequenceIsString = None

    def bridgeHMMs(self,HMM1,HMM2,transitionProbability):
        # if self.HMM != None:
        #     self.storeHMM()

        bridgeName = "%s-%s"%(HMM1.name,HMM2.name)

        if HMM1.sequenceIsString and HMM2.sequenceIsString:     #both HMMs have a known sequence
            bridgeSeq = self.bridgeSequences(HMM1,HMM2)
            self.buildDraftHMM(bridgeName, bridgeSeq, reverseEmitAllowed=True)
            # print("SPLICING 2 NT SEQS")

        else:
            bridgeSeq = self.bridgeMeans(HMM1,HMM2)
            self.buildDraftHMM(bridgeName, bridgeSeq,defaultEmitStdev=float(self.unknownBridgeStdev), reverseEmitAllowed=True)
            # print("SPLICING 2 MEAN SEQS")

        self.spliceHMMs(HMM1,self.HMMs[bridgeName],transitionProbability,reverseEmitAllowed=False)
        self.spliceHMMs(self.HMMs[bridgeName],HMM2,1)

    def bridgeMeans(self,HMM1,HMM2):
        if HMM1.sequenceIsString:
            seq = HMM1.sequence[-(self.k):]
            # print(seq)
            startMean = self.kmerCalc.calculateExpectedSignal(seq)[0]
            endMean = HMM2.expectedMeans[0]

        elif HMM2.sequenceIsString:
            seq = HMM2.sequence[:self.k]
            # print(seq)
            endMean = self.kmerCalc.calculateExpectedSignal(seq)[0]
            startMean = HMM1.expectedMeans[len(HMM1.expectedMeans)-1]
        else:
            startMean = HMM1.expectedMeans[len(HMM1.expectedMeans)-1]
            endMean = HMM2.expectedMeans[0]


        interval = (endMean - startMean)/6.0
        bridgeMeanList = [startMean+interval*(i+1) for i in range(5)]

        return bridgeMeanList

    def bridgeSequences(self,HMM1,HMM2):
        seq1 = self.HMMs[HMM1.name].sequence
        seq2 = self.HMMs[HMM2.name].sequence

        bridgeSequence = seq1[-(self.k-1):] + seq2[:self.k-1]

        # print(bridgeSequence)

        return bridgeSequence

    def addStartTransitions(self,HMM,transitionProbability):
        inPhase = HMM.HMM[0]

        # # start -> skip
        self.model.add_transition(self.startNode, inPhase["skip"], transitionProbability*self.transitions["start"]["skip"])

        # start -> emit
        self.model.add_transition(self.startNode, inPhase["emit"], transitionProbability*self.transitions["start"]["emit"])

        # start -> blip
        self.model.add_transition(self.startNode, inPhase["blip"], transitionProbability*self.transitions["start"]["blip"])

        # start -> drop
        self.model.add_transition(self.startNode, inPhase["drop"], transitionProbability*self.transitions["start"]["drop"])

    def addEndTransitions(self,HMM,transitionProbability):
        outPhase = HMM.HMM[len(HMM.HMM)-1]

        # # skip -> end
        self.model.add_transition(outPhase["skip"], self.endNode, transitionProbability*self.transitions["skip"]["next"])

        # emit -> end
        self.model.add_transition(outPhase["emit"], self.endNode, transitionProbability*self.transitions["emit"]["next"])

        # drop -> end
        self.model.add_transition(outPhase["drop"], self.endNode, transitionProbability*self.transitions["skip"]["next"])

        # blip -> end
        self.model.add_transition(outPhase["blip"], self.endNode, transitionProbability*self.transitions["blip"]["next"])

    def spliceHMMs(self,startHMM,endHMM,transitionProbability, reverseEmitAllowed=True):
        phaseA = startHMM.HMM[len(startHMM.HMM)-1]
        phaseB = endHMM.HMM[0]

        # if n > 0:
        if reverseEmitAllowed:
            # emit -> emit (prev) .
            self.model.add_transition(phaseB["emit"], phaseA["emit"], transitionProbability*self.transitions["emit"]["prev"])

        # emit (prev) -> emit .
        self.model.add_transition(phaseA["emit"], phaseB["emit"], transitionProbability*self.transitions["emit"]["next"])

        # blip (prev) -> blip .
        self.model.add_transition(phaseA["blip"], phaseB["blip"], transitionProbability*self.transitions["blip"]["next"])

        # emit (prev) -> skip .
        self.model.add_transition(phaseA["emit"], phaseB["skip"], transitionProbability*self.transitions["emit"]["skip"])

        # skip (prev) -> skip .
        # self.model.add_transition(phaseA["skip"], phaseB["skip"], transitionProbability*self.transitions["skip"]["next"])

        # skip (prev) -> emit .
        self.model.add_transition(phaseA["skip"], phaseB["emit"], transitionProbability*self.transitions["skip"]["emit"])

        # drop (prev) -> drop .
        self.model.add_transition(phaseA["drop"], phaseB["drop"], transitionProbability*self.transitions["drop"]["next"])


    def buildDraftHMM(self,name,sequence,defaultEmitStdev=5.0,stdevs=None,reverseEmitAllowed = True):
        self.sequence = sequence

        if type(sequence) == str:
            self.expectedMeans = self.kmerCalc.calculateExpectedSignal(sequence)
            self.sequenceIsString = True

        elif type(sequence) == list:
            if type(sequence[0]) == float:
                self.expectedMeans = sequence
                self.sequenceIsString = False
            else:
                sys.exit("ERROR: list of means provided not in float format")
        else:
            sys.exit("ERROR: Sequence provided is not nucleotide or list of means")

        # print(self.expectedMeans)

        self.HMM = list()
        nPhases = len(self.expectedMeans)

        for n in range(nPhases):
            print(n)
            self.HMM.append(self.newPhase()) #add template to HMM to contain each state

            if stdevs != None:
                emitNodeStdev = stdevs[n]
            else:
                emitNodeStdev = defaultEmitStdev

            emission = NormalDistribution(self.expectedMeans[n],emitNodeStdev)
            blipEmission = NormalDistribution(self.expectedMeans[n],3)
            dropEmission = NormalDistribution(10,3)

            #generate all the states for current phase
            self.HMM[n]["skip"] = State(None, "skip_%d_%s" % (n,name))    #silent
            self.HMM[n]["blip"] = State(blipEmission, "blip_%d_%s" % (n,name))
            self.HMM[n]["drop"] = State(dropEmission, "drop_%d_%s" % (n,name))
            self.HMM[n]["emit"] = State(emission, "emit_%d_%s" % (n,name))

            if n > 0:
                if reverseEmitAllowed:
                    # emit -> emit (prev)
                    self.model.add_transition(self.HMM[n]["emit"],self.HMM[n-1]["emit"],self.transitions["emit"]["prev"])

                # emit (prev) -> emit
                self.model.add_transition(self.HMM[n-1]["emit"],self.HMM[n]["emit"],self.transitions["emit"]["next"])

                # blip (prev) -> blip
                self.model.add_transition(self.HMM[n-1]["blip"],self.HMM[n]["blip"],self.transitions["blip"]["next"])

                # emit (prev) -> skip
                self.model.add_transition(self.HMM[n-1]["emit"],self.HMM[n]["skip"],self.transitions["emit"]["skip"])

                # skip (prev) -> skip
                # self.model.add_transition(self.HMM[n-1]["skip"],self.HMM[n]["skip"],self.transitions["skip"]["next"])

                # skip (prev) -> emit
                self.model.add_transition(self.HMM[n-1]["skip"],self.HMM[n]["emit"],self.transitions["skip"]["emit"])

                # drop (prev) -> drop
                self.model.add_transition(self.HMM[n-1]["drop"],self.HMM[n]["drop"],self.transitions["drop"]["next"])

            # emit -> self
            self.model.add_transition(self.HMM[n]["emit"], self.HMM[n]["emit"], self.transitions["emit"]["self"])

            # emit -> blip
            self.model.add_transition(self.HMM[n]["emit"], self.HMM[n]["blip"], self.transitions["emit"]["blip"])

            # emit -> drop
            self.model.add_transition(self.HMM[n]["emit"], self.HMM[n]["drop"], self.transitions["emit"]["drop"])

            # blip -> emit
            self.model.add_transition(self.HMM[n]["blip"], self.HMM[n]["emit"], self.transitions["blip"]["emit"])

            # drop -> emit
            self.model.add_transition(self.HMM[n]["drop"], self.HMM[n]["emit"], self.transitions["drop"]["emit"])

        self.storeHMM(name)
        self.printHMM(name)

    def sample(self):
        nplots = 4
        f, axes = plot.subplots(nplots, sharex=True, sharey=True)

        # signals = [HMM.expectedMeans for HMM in self.HMMs]
        signals= list()
        metaStates = set()
        paths = list()

        for i in range(nplots):
            signal,pathOutput = self.model.sample(path=True)
            signals.append(signal)
            print
            path = [(node.name).split('_')[-1] for node in pathOutput][1:-1]
            print(" -> ".join(path))

            for item in path:
                metaStates.add(item)

            paths.append(path)

        metaStates = sorted(list(metaStates))
        colors = {"B1": "blue",
                  "B2": "purple",
                  "fLead": "black",
                  "rLead": "brown",
                  "rLeadTail": "yellow",
                  "pA": "red",
                  "pT": "orange",
                  "hairpin": "green"}

        for s, signal in enumerate(signals):

            n = len(signal)

            try:
                length = 100 / float(n)
            except ZeroDivisionError:
                length = 0

            for i, kmerMean in enumerate(signal):
                x0 = (i) * length
                x1 = x0 + length
                y = kmerMean

                try: color = colors[paths[s][i]]
                except KeyError: color = "gray"


                if i > 0:
                    xprev = i * length
                    yprev = signal[i - 1]
                    axes[s].plot([xprev, x0], [yprev, y], color=color)

                axes[s].plot([x0, x1], [y, y], color=color)
                if s == 3: axes[s].set_ylabel("Current (pA)")
                # axes[s].set_ylim([0,150])
                # axes[s].set_yticks(np.linspace(0,150,4))

                axes[s].tick_params(axis='both', which='both',
                                    bottom='off', labelbottom='off',
                                    left='on', labelleft='on',
                                    right='off', labelright='off',
                                    top='off', labeltop='off')

        axes[0].set_title("Simulated Signals")
        plot.show()

    def printHMM(self,name):
        for p,phase in enumerate(self.HMMs[name].HMM):
            print("PHASE%d - "%p),
            for slot in phase:
                try:
                    print("%s"%(phase[slot].name)),
                except:
                    print("%s"%("None")),
            print

#------------------------------------------------------------------------------------------------------


def plotAlignedSignals(signalsSeparated, paths, colors):
    nplots = len(signalsSeparated)

    f, axes = plot.subplots(nplots, sharex=True, sharey=True)

    for s, signal in enumerate(signalsSeparated):
        # print(s)
        # print(len(signal))
        # print(paths[s])

        # colors = ["red","blue","orange","purple","green","yellow","brown","black","gray"]

        # if s >= len(signalsSeparated)-2:
        #     color = colors[-(s+1)]
        # else:
        #     color = "black"

        # plot.figure()
        # panel = plot.axes()
        # panel.plot(means)

        n = len(signal)
        length = 100/float(n)

        for i, kmerMean in enumerate(signal):
            x0 = (i+1)*length
            x1 = x0+length
            y = kmerMean

            try:
                color = colors[paths[s][i]]
            except KeyError:
                color = "gray"

            if i > 0:
                xprev = i*length+length
                yprev = signal[i-1]

                axes[s].plot([xprev, x0], [yprev, y], color=color)

            # print i,n,length,x0,x1,y

            axes[s].plot([x0, x1], [y, y], color=color)
            if s == 3: axes[s].set_ylabel("Current (pA)")
            # axes[s].set_ylim([0,150])
            # axes[s].set_yticks(np.linspace(0,150,4))

            axes[s].tick_params(axis='both', which='both',
                                bottom='off', labelbottom='off',
                                left='on', labelleft='on',
                                right='off', labelright='off',
                                top='off', labeltop='off')

    axes[0].set_title("Aligned Observed Signal")

    # f.savefig("barcodeComparison2.png",dpi=600)
    plot.show()


filename = "kmerMeans"


builder = metaHMMBuilder(filename)


fLeadMeans = list(map(float,[125,73,94,76,96,80,81,98,98,65,102,74.6, 62.4, 95.33333333, 99.75, 87.8, 92.75, 77.8]))
fLeadStdevs = list(map(float,[14,4.3,3,6,2.7,5,8,5,5,5,5,1.949358869, 3.361547263, 1.577350269, 1.707825128, 3.033150178, 4.573474245, 1.788854382]))

rLeadMeans = list(map(float, [95.75, 70.66666667, 62.75, 131, 98.75, 62, 80.75, 57.25, 47, 54.75]))
rLeadStdevs = list(map(float, [5.560275773, 3.055050463, 1.258305739, 5.715476066, 2.217355783, 4.082482905, 8.057087977, 3.095695937, 2.645751311, 1.258305739]))

rLeadTailMeans = list(map(float, [62.75, 81.66666667, 93, 77.75, 80]))
rLeadTailStdevs = list(map(float, [2.217355783, 4.041451884, 7.874007874, 44.68687354, 3.559026084]))

hairpinMeans = list(map(float,[103.5, 96.66666667, 69.875, 99.375, 106.875, 84.75, 109, 73.125, 91.75, 97.5, 83.5, 106.625, 73, 88.625, 101.375, 83.75, 97.375, 82.625, 62.75, 77.33333333, 94, 85.5, 94.375, 100.2857143, 75.14285714, 61.85714286, 98.375, 75, 84.125, 88.5, 95.75, 74.28571429, 134.25, 132.625, 129.375, 121.5, 113.375, 85, 74.875, 79.625, 86.25, 87.375, 131.5, 133, 125.875, 119, 110.25, 86.375, 80.375, 75.375, 83.25, 88.4, 99.625, 87.16666667, 106.75, 91.57142857, 77.375, 83.83333333, 104, 86.66666667, 73.33333333, 65.33333333, 84.33333333]))
hairpinStdevs = list(map(float,[5.042675027, 3.076794869, 6.937218463, 6.696214282, 8.42509008, 9.42956344, 5.631543813, 2.90012315, 5.092010549, 5.732115042, 3.271085447, 3.739270364, 2.563479778, 4.405759218, 6.926914382, 3.195979617, 6.412877669, 3.925648263, 5.063877679, 1.751190072, 3.023715784, 2.672612419, 3.020761493, 5.529143566, 6.256425269, 3.287784027, 3.335416016, 3.927922024, 8.626164518, 6.718843438, 1.908627031, 4.644505202, 10.41633333, 10.74293788, 4.897156609, 4.535573676, 6.345695954, 5.126959556, 2.474873734, 5.705573716, 3.412163118, 2.924648941, 5.903993806, 7.0305456, 7.772432971, 5.451081151, 6.541078985, 4.926241685, 4.56500665, 4.068608046, 3.991061441, 2.408318916, 2.924648941, 2.786873995, 4.131758533, 2.149196971, 2.924648941, 3.656045222, 4.242640687, 4.546060566, 3.326659987, 3.265986324, 3.011090611]))

homopolymerLength = 5
pA_means = list(map(float, [83.459321]*homopolymerLength))
pA_stdevs = list(map(float, [3.5]*homopolymerLength))

pT_means = list(map(float, [87.762283]*homopolymerLength))
pT_stdevs = list(map(float, [3.5]*homopolymerLength))

# builder.buildDraftHMM("lead",leadMeans,stdevs=leadStdevs)

n = 4.0

# builder.buildDraftHMM("a",pA_means,pA_stdevs)
#
# sys.exit()

transitions = {"fLead":{"pT":0.4,
                       "pA":0.4,
                       "B1":0.1,
                       "B2":0.1},
               "rLead":{"end":1/2.0,
                        "rLeadTail":1/2.0},
               "hairpin":{"pA":1/3.0,
                          "pT":1/3.0,
                          "B1":1/6.0,
                          "B2":1/6.0},
               "B1":{"B2":0.2,
                     "hairpin":0.1,
                     "pA":0.3,
                     "pT":0.3,
                     "rLead":0.1},
               "B2":{"B1":0.2,
                     "hairpin":0.1,
                     "pA":0.3,
                     "pT":0.3,
                     "rLead":0.1},
               "start":{"fLead":1},
               "pA":{"hairpin":0.1,
                     "pA":0.4,
                     "B1":0.2,
                     "B2":0.2,
                     "rLead":0.1},
               "pT": {"hairpin": 0.1,
                      "pT":0.4,
                      "B1": 0.2,
                      "B2": 0.2,
                      "rLead": 0.1},
               "rLeadTail":{"end":1}}

nodes = {"B1": "GGGTTCAATCAAGGGTTCAATCAAGGGTTCAATCAAGGGTTCAATCAAT",
         "B2": "TTGATTGAACCCTTGATTGAACCCTTGATTGAACCCTTGATTGAACCCAAAAAA",
         "fLead": (fLeadMeans,fLeadStdevs),
         "rLead": (rLeadMeans, rLeadStdevs),
         "rLeadTail": (rLeadTailMeans,rLeadTailStdevs),
         "pA": (pA_means,pA_stdevs),
         "pT": (pT_means,pT_stdevs),
         "hairpin": (hairpinMeans,hairpinStdevs)}

builder.buildMetaHMM(nodes,transitions)
builder.model.bake(verbose=True)
# builder.sample()
