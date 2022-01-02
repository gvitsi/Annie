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

plt.rc('font', family='Times', size=12)
import pylab
pylab.rcParams['figure.figsize'] = 10, 7.6


SIGNAL_DIR = "../Data/V3_5PE100ns/Pos0Data/"
BKG_DIR = "../Data/Calibration_2021/BKG/"

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
'''
def EventSelectionLosses(df,df_trig):
    print("TOTAL NUMBER OF EVENT TIME TANKS SET: " + str(len(set(df_trig['eventTimeTank']))))
    print("TOTAL NUMBER OF EVENT TIME TANKS, LIST: " + str(len(df_trig['eventTimeTank'])))
    print("TOTAL NUMBER OF ENTRIES: " + str(len(df_trig)))
    
    df_cleanSiPM = es.SingleSiPMPulses(df_trig)
    print("TOTAL NUMBER OF EVENTS W/ ONE PULSE IN EACH SIPM: " + str(len(set(df_cleanSiPM['eventTimeTank']))))
    df_SinglePulses = es.SingleSiPMPulses(df) #Clusters with a single SiPM pulse

    df_cleanSiPMDT = es.SingleSiPMPulsesDeltaT(df_trig,200) 
    print("TOTAL NUMBER OF EVENTS W/ ONE PULSE IN EACH SIPM, PEAKS WITHIN 200 NS: " + str(len(set(df_cleanSiPMDT['eventTimeTank']))))
    
    df_trig_cleanSiPM = df_trig.loc[(df_trig['SiPM1NPulses']==1) & (df_trig['SiPM2NPulses']==1)].reset_index(drop=True)
    df_cleanPrompt = es.NoPromptClusters(df_SinglePulses,2000)
    df_trig_cleanPrompt = es.NoPromptClusters_WholeFile(df_SinglePulses,df_trig_cleanSiPM,2000)
    print("TOTAL NUMBER OF EVENTS W/ ONE PULSE IN EACH SIPM AND NO PROMPT CLUSTER: " + str(len(set(df_trig_cleanPrompt['eventTimeTank']))))

    #Late burst cut
    df_trig_cleanWindow = es.NoBurst_WholeFile(df_cleanPrompt,df_trig_cleanPrompt,2000,150)
    print("TOTAL NUMBER OF EVENTS W/ CLEAN PROMPT AND NO BURST ABOVE 150 PE AND 2 MICROSECONDS: " + str(len(set(df_trig_cleanWindow['eventTimeTank']))))
'''
def PlotDemo(Bdf,Bdf_trig): 
    
    #print("EVENT SELECTION LOSSES FOR CENTRAL SOURCE RUN")
    #EventSelectionLosses(Sdf,Sdf_trig)
 
    #print("EVENT SELECTION LOSSES FOR BKG CENTRAL SOURCE RUN")
    #EventSelectionLosses(Bdf,Bdf_trig)
    Bdf['label'] = '0'
    print(Bdf.head())
    print("Bdf.shape: ", Bdf.shape)
    print("All columns are: ", Bdf.columns.values.tolist())
    #Bdf.to_csv("vars_DNN_Bkgd.csv",  index=False,float_format = '%.3f')
#    print(type(Bdf.hitDetID))

    #---- My Plots:
    Bdf_prompt=Bdf.loc[Bdf['clusterTime']<2000].reset_index(drop=True) #prompt events
    plt.hist(Bdf_prompt['clusterTime'],bins=100,range=(0,2000))
    plt.title("Prompt window Tank cluster times - no cuts")
    plt.xlabel("Cluster time [ns]")
    plt.savefig("plots_AmBe_Bkgd/time_prompt.png")
    plt.show()
#    

    Bdf_del=Bdf.loc[Bdf['clusterTime']>=2000].reset_index(drop=True) #delayed events
    plt.hist(Bdf_del['clusterTime'],bins=100,range=(2000,70000))
    plt.title("Delayed window Tank cluster times - no cuts")
    plt.xlabel("Cluster time [ns]")
    plt.savefig("plots_AmBe_Bkgd/time_del.png")
    plt.show()
#    

    #--- CB to cluster Time:   
    labels = {'title': 'Charge balance parameters in time window \n (Beam data, $t_{c}>=2 \, \mu s$)', 
             'xlabel': 'Cluster time (ns)', 'ylabel': 'Charge balance'}
    ranges = {'xbins': 58, 'ybins':50, 'xrange':[2000,60000],'yrange':[0,1]}
    abp.Make2DHist(Bdf_del,'clusterTime','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/CB_time_del.png")
    plt.show()
#    
    labels = {'title': 'Charge balance parameters in time window \n (Beam data, $t_{c}<2 \, \mu s$)',
              'xlabel': 'Cluster time (ns)', 'ylabel': 'Charge balance'}
    ranges = {'xbins': 20, 'ybins':50, 'xrange':[0,2000],'yrange':[0,1]}
    abp.Make2DHist(Bdf_prompt,'clusterTime','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/CB_time_prompt.png")
    plt.show()
#    

    #--- CB to clusterPE: 
    labels = {'title': 'Charge balance parameters in time window \n (Beam data, $t_{c}>=2 \, \mu s$)',
              'xlabel': 'Cluster PE', 'ylabel': 'Charge balance'}
    ranges = {'xbins': 58, 'ybins':50, 'xrange':[0,500],'yrange':[0,1]}
    abp.Make2DHist(Bdf_del,'clusterPE','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/CB_PE_del.png")
    plt.show()
#    
    labels = {'title': 'Charge balance parameters in time window \n (Beam data, $t_{c}<2 \, \mu s$)',
              'xlabel': 'Cluster PE', 'ylabel': 'Charge balance'}
    ranges = {'xbins': 20, 'ybins':50, 'xrange':[0,500],'yrange':[0,1]}
    abp.Make2DHist(Bdf_prompt,'clusterPE','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/CB_PE_prompt.png")
    plt.show()
#    

    #splitting to CB categories:
    #--- CB>=0.9 
    Bdf_prompt_highCB = Bdf_prompt.loc[Bdf_prompt['clusterChargeBalance']>=0.9].reset_index(drop=True) 
    Bdf_del_highCB = Bdf_del.loc[Bdf_del['clusterChargeBalance']>=0.9].reset_index(drop=True)

    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Beam data, $t_{c}<2 \, \mu s$) \n CB>=0.9 ',
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,200],'yrange':[0,200]}
    #abp.Make2DHist(Sdf_prompt_highCB,'clusterPE','clusterMaxPE',labels,ranges)
    abp.Make2DHist(Bdf_prompt_highCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/PE_maxPE_prompt_highCB.png")
    plt.show()
#    

    #PE = np.hstack(Sdf_del_highCB['hitPE'])
    #ID = np.hstack(Sdf_del_highCB['hitDetID'])
    #T = np.hstack(Sdf_del_highCB['hitT'])
    #maxPE_highCB = max(np.hstack(Sdf_prompt_highCB.hitPE))
    #print("maxPE_highCB ",maxPE_highCB," clusterMaxPE ",Sdf_prompt_highCB.clusterMaxPE)

    highCB_PE = np.hstack(Bdf_prompt_highCB.hitPE)
    highCB_DetID = np.hstack(Bdf_prompt_highCB.hitDetID)
#    highCB_PE = np.hstack(Sdf_del_highCB.hitPE)
#    highCB_DetID = np.hstack(Sdf_del_highCB.hitDetID)
    plt.hist2d(highCB_DetID,highCB_PE)
    plt.title("PE distribution for all hits in clusters, CB>=0.9)")
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.savefig("plots_AmBe_Bkgd/TubeID_PE_prompt_highCB.png")
    plt.show()
#    
    
    #--- 0.6<CB<0.9
    Bdf_prompt_upperCB = Bdf_prompt.loc[(Bdf_prompt['clusterChargeBalance']<0.9) & (Bdf_prompt['clusterChargeBalance']>=0.6)].reset_index(drop=True)
    Bdf_del_upperCB = Bdf_del.loc[(Bdf_del['clusterChargeBalance']<0.9) & (Bdf_prompt['clusterChargeBalance']>=0.6)].reset_index(drop=True)

    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Beam data, $t_{c}<2 \, \mu s$) \n 0.6<=CB<0.9',
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,200],'yrange':[0,200]}
    abp.Make2DHist(Bdf_prompt_upperCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/PE_maxPE_prompt_upperCB.png")
    plt.show()
#    

    upperCB_PE = np.hstack(Bdf_prompt_upperCB.hitPE)
    upperCB_DetID = np.hstack(Bdf_prompt_upperCB.hitDetID)
    plt.hist2d(upperCB_DetID,upperCB_PE)
    plt.title("PE distribution for all hits in clusters, 0.6=<CB<0.9)")
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.savefig("plots_AmBe_Bkgd/TubeID_PE_prompt_upperCB.png")
    plt.show()
    
#    

    #--- 0.4<CB<0.6
    Bdf_prompt_midCB = Bdf_prompt.loc[(Bdf_prompt['clusterChargeBalance']<0.6) & (Bdf_prompt['clusterChargeBalance']>=0.4)].reset_index(drop=True)
    Bdf_del_midCB = Bdf_del.loc[(Bdf_del['clusterChargeBalance']<0.6) & (Bdf_prompt['clusterChargeBalance']>=0.4)].reset_index(drop=True)
   
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Beam data, $t_{c}<2 \, \mu s$)\n 0.4<=CB<0.6',
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,200],'yrange':[0,200]}
    abp.Make2DHist(Bdf_prompt_midCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/PE_maxPE_prompt_midCB.png")
    plt.show()
#    
     
    midCB_PE = np.hstack(Bdf_prompt_midCB.hitPE)
    midCB_DetID = np.hstack(Bdf_prompt_midCB.hitDetID)
    plt.hist2d(midCB_DetID,midCB_PE)
    plt.title("PE distribution for all hits in clusters, 0.4=<CB<0.6)")
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.savefig("plots_AmBe_Bkgd/TubeID_PE_prompt_midCB.png")
    plt.show()
#    

    #--- CB<0.4
    Bdf_prompt_lowCB = Bdf_prompt.loc[Bdf_prompt['clusterChargeBalance']<0.4].reset_index(drop=True)
    Bdf_del_lowCB = Bdf_del.loc[Bdf_del['clusterChargeBalance']<0.4].reset_index(drop=True)
     
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Beam data, $t_{c}<2 \, \mu s$) \n CB<0.4',
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,200],'yrange':[0,200]}
    abp.Make2DHist(Bdf_prompt_lowCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/PE_maxPE_prompt_lowCB.png")
    plt.show()
#    
      
    lowCB_PE = np.hstack(Bdf_prompt_lowCB.hitPE)
    lowCB_DetID = np.hstack(Bdf_prompt_lowCB.hitDetID)
    plt.hist2d(lowCB_DetID,lowCB_PE)
    plt.title("PE distribution for all hits in clusters, CB<=0.4)")
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.savefig("plots_AmBe_Bkgd/TubeID_PE_prompt_lowCB.png")
    plt.show()
#    


    #splitting to CB categories - plots for delayed events:
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Beam data, $t_{c}>=2 \, \mu s$) \n CB>=0.9 ',
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 100, 'ybins':100, 'xrange':[0,100],'yrange':[0,100]}
    abp.Make2DHist(Bdf_del_highCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/TotPE_maxPEclust_high.png")
    plt.show()

    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Beam data, $t_{c}>=2 \, \mu s$) \n 0.6<=CB<0.9',
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 100, 'ybins':100, 'xrange':[0,100],'yrange':[0,100]}
    abp.Make2DHist(Bdf_del_upperCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/TotPE_maxPEclust_upper.png")
    plt.show()

    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Beam data, $t_{c}>=2 \, \mu s$)\n 0.4<=CB<0.6',
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 100, 'ybins':100, 'xrange':[0,100],'yrange':[0,100]}
    abp.Make2DHist(Bdf_del_midCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/TotPE_maxPEclust_mid.png")
    plt.show()

    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Beam data, $t_{c}>=2 \, \mu s$) \n CB<0.4',
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 100, 'ybins':100, 'xrange':[0,100],'yrange':[0,100]}
    abp.Make2DHist(Bdf_del_lowCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/TotPE_maxPEclust_low.png")
    plt.show()


if __name__=='__main__':
    #slist = glob.glob(SIGNAL_DIR+"*.ntuple.root")
    blist = glob.glob(BKG_DIR+"*.ntuple.root")

    #livetime_estimate = es.EstimateLivetime(slist)
    #print("SIGNAL LIVETIME ESTIMATE IN SECONDS IS: " + str(livetime_estimate))
    #livetime_estimate = es.EstimateLivetime(blist)
    #print("BKG LIVETIME ESTIMATE IN SECONDS IS: " + str(livetime_estimate))

    #mybranches = ['eventNumber','eventTimeTank','clusterTime','SiPMhitT','SiPMhitQ','SiPMhitAmplitude','clusterChargeBalance','clusterPE','SiPM1NPulses','SiPM2NPulses','SiPMNum','clusterHits']
    mybranches = ['clusterTime','hitT','hitQ','hitPE','hitDetID','clusterChargeBalance','clusterPE','clusterMaxPE','clusterHits'] #,'SiPMhitT','SiPMhitQ','SiPMhitAmplitude','SiPM1NPulses','SiPM2NPulses','SiPMNum']

    #SProcessor = rp.ROOTProcessor(treename="phaseIITankClusterTree")
    #for f1 in slist:
     #   SProcessor.addROOTFile(f1,branches_to_get=mybranches)
    #Sdata = SProcessor.getProcessedData()
    #Sdf = pd.DataFrame(Sdata)

    BProcessor = rp.ROOTProcessor(treename="phaseIITankClusterTree")
    for f1 in blist:
        BProcessor.addROOTFile(f1,branches_to_get=mybranches)
    Bdata = BProcessor.getProcessedData()
    Bdf = pd.DataFrame(Bdata)

    #SProcessor = rp.ROOTProcessor(treename="phaseIITriggerTree")
    #for f1 in slist:
     #   SProcessor.addROOTFile(f1,branches_to_get=mybranches)
    #Sdata = SProcessor.getProcessedData()
    #Sdf_trig = pd.DataFrame(Sdata)

    BProcessor = rp.ROOTProcessor(treename="phaseIITriggerTree")
    for f1 in blist:
        BProcessor.addROOTFile(f1,branches_to_get=mybranches)
    Bdata = BProcessor.getProcessedData()
    Bdf_trig = pd.DataFrame(Bdata)

    PlotDemo(Bdf,Bdf_trig)


