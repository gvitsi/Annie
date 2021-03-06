#This script has functions needed to make a Cumulative Distribution Plot from
#Different variables output in  PhaseIITreeMaker root file.

import glob

import sys
import uproot
import lib.ROOTProcessor as rp
import lib.EventSelection as es
import lib.ProfileLikelihoodBuilder as plb
import lib.AmBePlots as abp
import lib.BeamPlots as bp
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as scp
import numpy as np
import scipy.misc as scm
from pylab import figure, axes, pie, title, show
from sklearn.utils import shuffle

plt.rc('font', family='Times', size=12)
import pylab
pylab.rcParams['figure.figsize'] = 10, 7.6


z=-50
SIGNAL_DIR = f"../Data/Calibration_2021/Signal/{z}/"
BKG_DIR = f"../Data/Calibration_2021/BKG/{z}/"

#SIGNAL_DIR = "/Users/edrakopo/work/ANNIETools_ntuples/ANNIETools/ANNIENtupleAnalysis/Data/V3_5PE100ns/Pos0Data/"
#BKG_DIR = "/Users/edrakopo/work/ANNIETools_ntuples/ANNIETools/ANNIENtupleAnalysis/Data/V3_5PE100ns/BkgPos0Data/"

#PEPERMEV = 12.
#expoPFlat= lambda x,C1,tau,mu,B: C1*np.exp(-(x-mu)/tau) + B
#mypoisson = lambda x,mu: (mu**x)*np.exp(-mu)/scm.factorial(x)

def GetDataFrame(mytreename,mybranches,filelist):
    RProcessor = rp.ROOTProcessor(treename=mytreename)
    for f1 in filelist:
        RProcessor.addROOTFile(f1,branches_to_get=mybranches)
    data = RProcessor.getProcessedData()
    df = pd.DataFrame(data)
    return df

def PlotDemo(Sdf,Bdf,Sdf_trig,Bdf_trig): 
   Sdf['label'] = '1'
   print("----- Signal------")
   print(Sdf.head())
   print("Sdf.shape: ", Sdf.shape)
   print("All columns are: ", Sdf.columns.values.tolist())
   Sdf.to_csv("vars_DNN_Signal.csv",  index=False,float_format = '%.3f')
   #print(type(Sdf.hitDetID))

   Bdf['label'] = '0'
   Bdf = shuffle(Bdf, random_state=0)
   print("----- Bkgd------")
   print(Bdf.head())
   print("Bdf.shape: ", Bdf.shape)
   print("All columns are: ", Bdf.columns.values.tolist())
   Bdf.to_csv("vars_DNN_Bkgd.csv",  index=False,float_format = '%.3f')
   #print(type(Bdf.hitDetID))
    
   data = pd.concat((Sdf,Bdf[:27645]))
   print("----- Signal+Bkgd------")
   #data['hitDetID'].to_csv("testing.csv")

   data['hitDetID'] = [','.join(str(y) for y in x) for x in data['hitDetID']] #dropping brackets in pd.Series
   data['hitPE'] = [','.join(str(y) for y in x) for x in data['hitPE']]
   #data['hitQ'] = [','.join(str(y) for y in x) for x in data['hitQ']]
   #data['hitT'] = [','.join(str(y) for y in x) for x in data['hitT']]
   print(data.head())
   #print(data.tail())
   print("data.shape: ", data.shape)
  
   #randomly shuffle the data
   data = shuffle(data, random_state=0) 
   print("after shuffling: ", data.head())
   print("data.shape: ", data.shape)
   data.to_csv("labels_DNN_Signal_BkgdNEW.csv",  index=False,float_format = '%.3f', sep=",")   
   data.drop(['label'], axis=1).to_csv("vars_DNN_Signal_BkgdNEW.csv",header=False,index=False,float_format = '%.3f', sep=",")

   #-------- selecting only prompt events as signal: --------#
   print("Selecting only prompt events (t<2us) as signal")
   Sdf_prompt=Sdf.loc[Sdf['clusterTime']<2000].reset_index(drop=True)
   print(Sdf_prompt.head())
   print("Sdf_prompt.shape: ", Sdf_prompt.shape)
   data2 = pd.concat((Sdf_prompt,Bdf[:3570]))
   data2['hitDetID'] = [','.join(str(y) for y in x) for x in data2['hitDetID']]
   data2['hitPE'] = [','.join(str(y) for y in x) for x in data2['hitPE']]
   print("data2.shape: ", data2.shape)

   #randomly shuffle the data
   data2 = shuffle(data2, random_state=0)
   print("after shuffling: ", data2.head())
   print("data2.shape: ", data2.shape)
   data2.to_csv("labels_DNN_Signal_Bkgd_promptNEW.csv",  index=False,float_format = '%.3f', sep=",")
   data2.drop(['label'], axis=1).to_csv("vars_DNN_Signal_Bkgd_promptNEW.csv",header=False,index=False,float_format = '%.3f', sep=",")

   #-------- selecting only delayed events as signal: --------#
   print("Selecting only delayed events (t>=2us) as signal")
   Sdf_del=Sdf.loc[Sdf['clusterTime']>=2000].reset_index(drop=True)
   print(Sdf_del.head())
   print("Sdf_del.shape: ", Sdf_del.shape)
   data3 = pd.concat((Sdf_del,Bdf[:24075]))
   data3['hitDetID'] = [','.join(str(y) for y in x) for x in data3['hitDetID']]
   data3['hitPE'] = [','.join(str(y) for y in x) for x in data3['hitPE']]
   print("data3.shape: ", data3.shape)

   #randomly shuffle the data
   data3 = shuffle(data3, random_state=0)
   print("after shuffling: ", data3.head())
   print("data3.shape: ", data3.shape)
   data3.to_csv("labels_DNN_Signal_Bkgd_delNEW.csv",  index=False,float_format = '%.3f', sep=",")
   data3.drop(['label'], axis=1).to_csv("vars_DNN_Signal_Bkgd_delNEW.csv",header=False,index=False,float_format = '%.3f', sep=",")

   '''
   Sdf['label'] = '1'
   Bdf['label'] = '0'
   labels = pd.concat((Sdf,Bdf))
   assert(data.shape[0]==labels.shape[0])
   labels.to_csv("labels_DNN_Signal_Bkgd.csv",  index=False,float_format = '%.3f', sep=",")
   '''
if __name__=='__main__':
    slist = glob.glob(SIGNAL_DIR+"*.ntuple.root")
    blist = glob.glob(BKG_DIR+"*.ntuple.root")

    livetime_estimate = es.EstimateLivetime(slist)
    print("SIGNAL LIVETIME ESTIMATE IN SECONDS IS: " + str(livetime_estimate))
    livetime_estimate = es.EstimateLivetime(blist)
    print("BKG LIVETIME ESTIMATE IN SECONDS IS: " + str(livetime_estimate))

    #mybranches = ['eventNumber','eventTimeTank','clusterTime','SiPMhitT','SiPMhitQ','SiPMhitAmplitude','clusterChargeBalance','clusterPE','SiPM1NPulses','SiPM2NPulses','SiPMNum','clusterHits']
    mybranches = ['eventNumber','eventTimeTank','clusterTime','hitT','hitQ','hitPE','hitDetID','clusterChargeBalance','clusterPE','clusterMaxPE'    ,'clusterHits', 'SiPMhitT','SiPMhitQ','SiPMhitAmplitude','SiPM1NPulses','SiPM2NPulses','SiPMNum']
    #mybranches = ['clusterTime','hitT','hitQ','hitPE','hitDetID','clusterChargeBalance','clusterPE','clusterMaxPE','clusterHits']
    #mybranches = ['clusterTime','hitDetID','hitPE','clusterChargeBalance','clusterPE','clusterMaxPE','clusterHits']

    SProcessor = rp.ROOTProcessor(treename="phaseIITankClusterTree")
    for f1 in slist:
        SProcessor.addROOTFile(f1,branches_to_get=mybranches)
    Sdata = SProcessor.getProcessedData()
    Sdf = pd.DataFrame(Sdata)

    BProcessor = rp.ROOTProcessor(treename="phaseIITankClusterTree")
    for f1 in blist:
        BProcessor.addROOTFile(f1,branches_to_get=mybranches)
    Bdata = BProcessor.getProcessedData()
    Bdf = pd.DataFrame(Bdata)

    SProcessor = rp.ROOTProcessor(treename="phaseIITriggerTree")
    for f1 in slist:
        SProcessor.addROOTFile(f1,branches_to_get=mybranches)
    Sdata = SProcessor.getProcessedData()
    Sdf_trig = pd.DataFrame(Sdata)

    BProcessor = rp.ROOTProcessor(treename="phaseIITriggerTree")
    for f1 in blist:
        BProcessor.addROOTFile(f1,branches_to_get=mybranches)
    Bdata = BProcessor.getProcessedData()
    Bdf_trig = pd.DataFrame(Bdata)

    PlotDemo(Sdf,Bdf,Sdf_trig,Bdf_trig)


