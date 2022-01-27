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
import pylab


#plt.rcParams["font.family"] = "Times"
pylab.rcParams['font.size'] = 12
pylab.rcParams['figure.figsize'] = 10, 7.6

z=0
SIGNAL_DIR = "../Data/Calibration_2021/Signal/{z}/"
BKG_DIR = f"../Data/Calibration_2021/BKG/{z}/"

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
    #df_trig_cleanWindow = es.NoBurst_WholeFile(df,df_trig,2000,150)
    print("TOTAL NUMBER OF EVENTS W/ CLEAN PROMPT AND NO BURST ABOVE 150 PE AND 2 MICROSECONDS: " + str(len(set(df_trig_cleanWindow['eventTimeTank']))))

   
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
    Bdf_del=Bdf.loc[Bdf['clusterTime']>=2000].reset_index(drop=True) #delayed events
    '''
    plt.hist(Bdf_prompt['clusterTime'],bins=100,range=(0,2000))
    plt.subplots_adjust(left=0.15, right=0.92, bottom=0.14, top=0.92, hspace=0.05, wspace=0.5)
    plt.title("Prompt window Tank cluster times - no cuts")
    plt.xlabel("Cluster time [ns]")
    plt.savefig("plots_AmBe_Bkgd/{}/time_prompt.png".format(z))
    plt.show()
#    

    
    plt.hist(Bdf_del['clusterTime'],bins=100,range=(2000,70000))
    plt.subplots_adjust(left=0.15, right=0.90, bottom=0.14, top=0.92, hspace=0.05, wspace=0.5)
    plt.title("Delayed window Tank cluster times - no cuts")
    plt.xlabel("Cluster time [ns]")
    plt.savefig("plots_AmBe_Bkgd/{}/time_del.png".format(z))
    plt.show()
#    
'''
    #--- CB to cluster Time:  
    labels = {'title': 'Charge balance parameters in time window \n (Bkgd data, $0<t_{c}<=60 \, \mu s$)\n' f"for position Z={z}cm", 'xlabel': 'Cluster time (ns)', 'ylabel': 'Charge balance'} #total events
    plt.subplots_adjust(left=0.19, right=0.94, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    ranges = {'xbins': 200, 'ybins':100, 'xrange':[0,60000],'yrange':[0,1]}
    abp.Make2DHist(Bdf,'clusterTime','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/CB_time_total.png".format(z))
    plt.show()
    
    plt.subplots_adjust(left=0.17, right=0.94, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Charge balance parameters in time window \n (Bkgd data, $t_{c}>=2 \, \mu s$)\n' f"for position Z={z}cm", 
             'xlabel': 'Cluster time (ns)', 'ylabel': 'Charge balance'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[2000,30000],'yrange':[0,1]}
    abp.Make2DHist(Bdf_del,'clusterTime','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/CB_time_del.png".format(z))
    plt.show()
#    
    plt.subplots_adjust(left=0.19, right=0.94, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Charge balance parameters in time window \n (Bkgd data, $t_{c}<2 \, \mu s$)\n' f"for position Z={z}cm",
              'xlabel': 'Cluster time (ns)', 'ylabel': 'Charge balance'}
    ranges = {'xbins': 200, 'ybins':100, 'xrange':[0,2000],'yrange':[0,1]}
    abp.Make2DHist(Bdf_prompt,'clusterTime','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/CB_time_prompt.png".format(z))
    plt.show()
#    

    #--- CB to clusterPE: 
    plt.subplots_adjust(left=0.15, right=0.98, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Charge balance parameters in time window \n (Bkgd data, $0<t_{c}<=60 \, \mu s$)\n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Charge balance'} #total events
    ranges = {'xbins': 200, 'ybins':100, 'xrange':[0,300],'yrange':[0,1]}
    abp.Make2DHist(Bdf,'clusterPE','clusterChargeBalance',labels,ranges)
    plt.savefig(f"plots_AmBe_Bkgd/{z}/CB_PE_total.png")
    plt.show()
    
    plt.subplots_adjust(left=0.15, right=0.98, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Charge balance parameters in time window \n (Bkgd data, $t_{c}>=2 \, \mu s$)\n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Charge balance'} #delayed events
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,200],'yrange':[0,1]}
    abp.Make2DHist(Bdf_del,'clusterPE','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/CB_PE_del.png".format(z))
    plt.show()
#    
    plt.subplots_adjust(left=0.15, right=0.98, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2)
    labels = {'title': 'Charge balance parameters in time window \n (Bkgd data, $t_{c}<2 \, \mu s$)\n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Charge balance'} #prompt events
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,300],'yrange':[0,1]}
    abp.Make2DHist(Bdf_prompt,'clusterPE','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/CB_PE_prompt.png".format(z))
    plt.show()
#    

    #splitting to CB categories:
    #--- CB>=0.9 
    Bdf_prompt_highCB = Bdf_prompt.loc[Bdf_prompt['clusterChargeBalance']>=0.9].reset_index(drop=True) 
    Bdf_del_highCB = Bdf_del.loc[Bdf_del['clusterChargeBalance']>=0.9].reset_index(drop=True)

    plt.subplots_adjust(left=0.15, right=1, bottom=0.15, top=0.8, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Bkgd data, $t_{c}<2 \, \mu s$) \n CB>=0.9 ' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,140],'yrange':[0,60]}
    #abp.Make2DHist(Sdf_prompt_highCB,'clusterPE','clusterMaxPE',labels,ranges)
    abp.Make2DHist(Bdf_prompt_highCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/PE_maxPE_prompt_highCB.png".format(z))
    plt.show()
    
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Bkgd data, $t_{c}>=2 \, \mu s$) CB>=0.9 \n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,50],'yrange':[0,50]}
    plt.subplots_adjust(left=0.12, right=1, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2)
    abp.Make2DHist(Bdf_del_highCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/PE_maxPEclust_del_high.png".format(z))
    plt.show()
#    

    #PE = np.hstack(Sdf_del_highCB['hitPE'])
    #ID = np.hstack(Sdf_del_highCB['hitDetID'])
    #T = np.hstack(Sdf_del_highCB['hitT'])
    #maxPE_highCB = max(np.hstack(Sdf_prompt_highCB.hitPE))
    #print("maxPE_highCB ",maxPE_highCB," clusterMaxPE ",Sdf_prompt_highCB.clusterMaxPE)

    highCB_PE_prompt = np.hstack(Bdf_prompt_highCB.hitPE)
    highCB_DetID_prompt = np.hstack(Bdf_prompt_highCB.hitDetID)
    highCB_PE_del = np.hstack(Bdf_del_highCB.hitPE)
    highCB_DetID_del = np.hstack(Bdf_del_highCB.hitDetID)
    #labels = {'title': 'PE distribution for all hits in clusters, CB>=0.9) ','xlabel': 'Tube ID', 'ylabel': 'PE'}
    #ranges = {'xbins': 200, 'ybins':100, 'xrange':[0,500],'yrange':[0,100]}
    #abp.Make2DHist(Bdf_prompt_highCB, 'hitDetID', 'hitPE', labels, ranges)
    plt.subplots_adjust(left=0.14, right=0.97, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.hist2d(highCB_DetID_prompt,highCB_PE_prompt, bins=(200,136), range=[[330,465],[0,50]],cmap=plt.cm.turbo)
    plt.colorbar()
    plt.title("Prompt PE distribution for all hits in clusters, \nCB>=0.9) " f"for position Z={z}cm")
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.savefig("plots_AmBe_Bkgd/{}/TubeID_PE_prompt_highCB.png".format(z))
    plt.show()
    
    plt.subplots_adjust(left=0.14, right=0.97, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.hist2d(highCB_DetID_del,highCB_PE_del, bins=(200,136), range=[[330,465],[0,20]],cmap=plt.cm.turbo)
    plt.colorbar()
    plt.title("Delayed PE distribution for all hits in clusters, \nCB>=0.9) " f"for position Z={z}cm")
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.savefig("plots_AmBe_Bkgd/{}/TubeID_PE_del_highCB.png".format(z))
    plt.show()
#    
    
    #--- 0.6<CB<0.9
    Bdf_prompt_upperCB = Bdf_prompt.loc[(Bdf_prompt['clusterChargeBalance']<0.9) & (Bdf_prompt['clusterChargeBalance']>=0.6)].reset_index(drop=True)
    Bdf_del_upperCB = Bdf_del.loc[(Bdf_del['clusterChargeBalance']<0.9) & (Bdf_prompt['clusterChargeBalance']>=0.6)].reset_index(drop=True)

    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Bkgd data, $t_{c}<2 \, \mu s$)\n 0.6<=CB<0.9  ' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,140],'yrange':[0,60]}
    
    plt.subplots_adjust(left=0.15, right=1, bottom=0.15, top=0.8, hspace=0.05, wspace=0.2)
    abp.Make2DHist(Bdf_prompt_upperCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/PE_maxPE_prompt_upperCB.png".format(z))
    plt.show()
    
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Bkgd data, $t_{c}>=2 \, \mu s$)  0.6<=CB<0.9\n ' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,30],'yrange':[0,50]}
    plt.subplots_adjust(left=0.12, right=1, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2)
    abp.Make2DHist(Bdf_del_upperCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/PE_maxPEclust_del_upper.png".format(z))
    plt.show()
#    

    upperCB_PE_prompt = np.hstack(Bdf_prompt_upperCB.hitPE)
    upperCB_DetID_prompt = np.hstack(Bdf_prompt_upperCB.hitDetID)
    #labels = {'title': 'PE distribution for all hits in clusters, 0.6=<CB<0.9) ','xlabel': 'Tube ID', 'ylabel': 'PE'}
    #ranges = {'xbins': 200, 'ybins':100, 'xrange':[0,500],'yrange':[0,100]}
    #abp.Make2DHist(Bdf_prompt_highCB, 'hitDetID', 'hitPE', labels, ranges)
    plt.subplots_adjust(left=0.14, right=0.97, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2)
    plt.hist2d(upperCB_DetID_prompt, upperCB_PE_prompt, bins=(200,136), range=[[330,465],[0,50]],cmap=plt.cm.turbo)
    #plt.hist2d(upperCB_DetID,upperCB_PE, labels, ranges)
    plt.title("Prompt PE distribution for all hits in clusters, \n 0.6=<CB<0.9) | " f"for position Z={z}cm")
    plt.colorbar()
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.savefig("plots_AmBe_Bkgd/{}/TubeID_PE_prompt_upperCB.png".format(z))
    plt.show()
    
    upperCB_PE_del = np.hstack(Bdf_del_upperCB.hitPE)
    upperCB_DetID_del = np.hstack(Bdf_del_upperCB.hitDetID)
    
    plt.subplots_adjust(left=0.12, right=1, bottom=0.13, top=0.80, hspace=0.05, wspace=0.2)
    plt.hist2d(upperCB_DetID_del, upperCB_PE_del,bins=(200,136), range=[[330,465],[0,20]],cmap=plt.cm.turbo)
    plt.title("Delayed PE distribution \n for all hits in clusters, \n 0.6=<CB<0.9)" f"for position Z={z}cm")
    plt.colorbar()
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.savefig("plots_AmBe_Bkgd/{}/TubeID_PE_del_upperCB.png".format(z))
    plt.show()
    
#    

    #--- 0.4<CB<0.6
    Bdf_prompt_midCB = Bdf_prompt.loc[(Bdf_prompt['clusterChargeBalance']<0.6) & (Bdf_prompt['clusterChargeBalance']>=0.4)].reset_index(drop=True)
    Bdf_del_midCB = Bdf_del.loc[(Bdf_del['clusterChargeBalance']<0.6) & (Bdf_prompt['clusterChargeBalance']>=0.4)].reset_index(drop=True)
   
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Bkgd data, $t_{c}<2 \, \mu s$)  0.4<=CB<0.6\n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,140],'yrange':[0,60]}
    plt.subplots_adjust(left=0.15, right=1, bottom=0.15, top=0.8, hspace=0.05, wspace=0.2)
    abp.Make2DHist(Bdf_prompt_midCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/PE_maxPE_prompt_midCB.png".format(z))
    plt.show()
    
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Bkgd data, $t_{c}>=2 \, \mu s$) 0.4<=CB<0.6\n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,30],'yrange':[0,50]}
    plt.subplots_adjust(left=0.12, right=1, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2)
    abp.Make2DHist(Bdf_del_midCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/PE_maxPEclust_del_mid.png".format(z))
    plt.show()
#    
     
    midCB_PE_prompt = np.hstack(Bdf_prompt_midCB.hitPE)
    midCB_DetID_prompt = np.hstack(Bdf_prompt_midCB.hitDetID)
    midCB_PE_del = np.hstack(Bdf_del_midCB.hitPE)
    midCB_DetID_del = np.hstack(Bdf_del_midCB.hitDetID)
    
    plt.subplots_adjust(left=0.14, right=0.97, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2)
    plt.hist2d(midCB_DetID_prompt, midCB_PE_prompt, bins=(200,136), range=[[330,465],[0,50]],cmap=plt.cm.turbo)
    plt.title("Prompt PE distribution for all hits in clusters, \n 0.4=<CB<0.6)" f"for position Z={z}cm")
    plt.colorbar()
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.savefig("plots_AmBe_Bkgd/{}/TubeID_PE_prompt_midCB.png".format(z))
    plt.show()
    
    
    plt.subplots_adjust(left=0.14, right=0.97, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2)
    plt.hist2d(midCB_DetID_del, midCB_PE_del ,bins=(200,136), range=[[330,465],[0,50]],cmap=plt.cm.turbo)
    plt.title("Delayed PE distribution for all hits in clusters, \n 0.4=<CB<0.6)" f"for position Z={z}cm")
    plt.colorbar()
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.savefig("plots_AmBe_Bkgd/{}/TubeID_PE_del_midCB.png".format(z))
    plt.show()
#    

    #--- CB<0.4
    Bdf_prompt_lowCB = Bdf_prompt.loc[Bdf_prompt['clusterChargeBalance']<0.4].reset_index(drop=True)
    Bdf_del_lowCB = Bdf_del.loc[Bdf_del['clusterChargeBalance']<0.4].reset_index(drop=True)
     
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Bkgd data, $t_{c}<2 \, \mu s$) CB<0.4\n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    plt.subplots_adjust(left=0.15, right=1, bottom=0.15, top=0.80, hspace=0.05, wspace=0.2)
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,140],'yrange':[0,60]}
    abp.Make2DHist(Bdf_prompt_lowCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/PE_maxPE_prompt_lowCB.png".format(z))
    plt.show()
    
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (Bkgd data, $t_{c}>=2 \, \mu s$) CB<0.4\n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,50],'yrange':[0,50]}
    plt.subplots_adjust(left=0.12, right=1, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2)
    abp.Make2DHist(Bdf_del_lowCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/PE_maxPEclust_del_low.png".format(z))
    plt.show()
#    
      
    lowCB_PE_prompt = np.hstack(Bdf_prompt_lowCB.hitPE)
    lowCB_DetID_prompt = np.hstack(Bdf_prompt_lowCB.hitDetID)
    lowCB_PE_del = np.hstack(Bdf_prompt_lowCB.hitPE)
    lowCB_DetID_del = np.hstack(Bdf_prompt_lowCB.hitDetID)
    
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2)
    plt.hist2d(lowCB_DetID_prompt,lowCB_PE_prompt, bins=(200,136), range=[[330,465],[0,50]],cmap=plt.cm.turbo)
    plt.title("Prompt PE distribution for all hits in clusters,\n CB<=0.4) | " f"for position Z={z}cm")
    plt.colorbar()
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.savefig("plots_AmBe_Bkgd/{}/TubeID_PE_prompt_lowCB.png".format(z))
    plt.show()
    
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2)
    plt.hist2d(lowCB_DetID_del,lowCB_PE_del, bins=(200,136), range=[[330,465],[0,50]],cmap=plt.cm.turbo)
    plt.title("Delayed PE distribution for all hits in clusters,\n CB<=0.4) | " f"for position Z={z}cm")
    plt.colorbar()
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.savefig("plots_AmBe_Bkgd/{}/TubeID_PE_del_lowCB.png".format(z))
    plt.show()
#    


    
    
    #----- Plots only for 0.4 < CB < 0.6 -----#:
    plt.subplots_adjust(left=0.15, right=0.98, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Charge balance parameters in time window \n (Bkgd data, $t_{c}=>2 \, \mu s$, 0.4 < CB < 0.6)\n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Charge balance'} #delayed events
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,300],'yrange':[0,1]}
    abp.Make2DHist(Bdf_del_midCB,'clusterPE','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/CB_PE_del_midCB.png".format(z))
    plt.show()
    
    plt.subplots_adjust(left=0.15, right=0.98, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Charge balance parameters in time window \n (Bkgd data, $t_{c}<2 \, \mu s$, 0.4 < CB < 0.6)\n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Charge balance'} #prompt events
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,300],'yrange':[0,1]}
    abp.Make2DHist(Bdf_prompt_midCB,'clusterPE','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/CB_PE_prompt_midCB.png".format(z))
    plt.show()
    
    plt.subplots_adjust(left=0.17, right=0.94, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Charge balance parameters in time window \n (Bkgd data, $t_{c}>=2 \, \mu s$, 0.4 < CB < 0.6)\n' f"for position Z={z}cm", 
             'xlabel': 'Cluster time (ns)', 'ylabel': 'Charge balance'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[2000,30000],'yrange':[0,1]}
    abp.Make2DHist(Bdf_del_midCB,'clusterTime','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/CB_time_del_midCB.png".format(z))
    plt.show()
#    
    plt.subplots_adjust(left=0.19, right=0.94, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Charge balance parameters in time window \n (Bkgd data, $t_{c}<2 \, \mu s$, 0.4 < CB < 0.6)\n' f"for position Z={z}cm",
              'xlabel': 'Cluster time (ns)', 'ylabel': 'Charge balance'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,2000],'yrange':[0,1]}
    abp.Make2DHist(Bdf_prompt_midCB,'clusterTime','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe_Bkgd/{}/CB_time_prompt_midCB.png".format(z))
    plt.show()
    
    
    
    
    
    #----- 1D Histogram for all CB areas -----#:
    
    #------prompt charts------#
    #--- CB>=0.9 
    '''
    Bdf_prompt_highCB['clusterPE'].plot(kind = "hist", density = False, bins=100,range=(0,500))
    plt.xlabel("Cluster PE")
    Bdf_prompt_highCB['clusterPE'].plot(kind = "kde")
    plt.show()
    '''
    plt.hist(Bdf_prompt_highCB['clusterPE'],bins=100,range=(0,500), histtype='step', density=True )
    #print(n1.shape, bins1.shape,len(bins1), max(n1), bins1.max(), n1.max(),)
    #y=plt.hist(Bdf_prompt_highCB['clusterPE'],bins=100,range=(0,500),histtype='step', weights=(1/max(n1))*(bins1.shape-1))
    plt.title("Prompt Cluster PE, CB>=0.9\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe_Bkgd/{z}/PromptHistClusterPE_highCB.png")
    plt.show()
    
    #--- 0.6<=CB<0.9
    plt.hist(Bdf_prompt_upperCB['clusterPE'],bins=100,range=(0,500),histtype='step', density=True)
    plt.title("Prompt Cluster PE, 0.6<=CB<0.9\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe_Bkgd/{z}/PromptHistClusterPE_upperCB.png")
    plt.show()
    
    #--- 0.4<=CB<0.6
    plt.hist(Bdf_prompt_midCB['clusterPE'],bins=100,range=(0,500),histtype='step', density=True)
    plt.title("Prompt Cluster PE, 0.4<=CB<0.6\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe_Bkgd/{z}/PromptHistClusterPE_midCB.png")
    plt.show()
    
    #--- CB<0.4
    plt.hist(Bdf_prompt_lowCB['clusterPE'],bins=100,range=(0,500),histtype='step', density=True)
    plt.title("Prompt Cluster PE, CB<0.4\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe_Bkgd/{z}/PromptHistClusterPE_lowCB.png")
    plt.show()
    
    
    #-----delayed charts------#
    #--- CB>=0.9 
    plt.hist(Bdf_del_highCB['clusterPE'],bins=100,range=(2000,20000),histtype='step', density=True)
    plt.title("Delayed Cluster PE, CB>=0.9\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe_Bkgd/{z}/DelHistClusterPE_highCB.png")
    plt.show()
    
    #--- 0.6<=CB<0.9
    plt.hist(Bdf_del_upperCB['clusterPE'],bins=100,range=(2000,20000),histtype='step', density=True)
    plt.title("Delayed Cluster PE, 0.6<=CB<0.9\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe_Bkgd/{z}/DelHistClusterPE_upperCB.png")
    plt.show()
    
    #--- 0.4<=CB<0.6
    plt.hist(Bdf_del_midCB['clusterPE'],bins=100,range=(2000,20000),histtype='step', density=True)
    plt.title("Delayed Cluster PE, 0.4<=CB<0.6\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe_Bkgd/{z}/DelHistClusterPE_midCB.png")
    plt.show()
    
    #--- CB<0.4
    plt.hist(Bdf_del_lowCB['clusterPE'],bins=100,range=(2000,20000),histtype='step', density=True)
    plt.title("Delayed Cluster PE, CB<0.4\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe_Bkgd/{z}/DelHistClusterPE_lowCB.png")
    plt.show()
    


if __name__=='__main__':
    #slist = glob.glob(SIGNAL_DIR+"*.ntuple.root")
    blist = glob.glob(BKG_DIR+"*.ntuple.root")

    #livetime_estimate = es.EstimateLivetime(slist)
    #print("SIGNAL LIVETIME ESTIMATE IN SECONDS IS: " + str(livetime_estimate))
    #livetime_estimate = es.EstimateLivetime(blist)
    #print("BKG LIVETIME ESTIMATE IN SECONDS IS: " + str(livetime_estimate))

    #mybranches = ['eventNumber','eventTimeTank','clusterTime','SiPMhitT','SiPMhitQ','SiPMhitAmplitude','clusterChargeBalance','clusterPE','SiPM1NPulses','SiPM2NPulses','SiPMNum','clusterHits']
    mybranches = ['eventNumber','eventTimeTank','clusterTime','hitT','hitQ','hitPE','hitDetID','clusterChargeBalance','clusterPE','clusterMaxPE','clusterHits','SiPMhitT','SiPMhitQ','SiPMhitAmplitude','SiPM1NPulses','SiPM2NPulses','SiPMNum']

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

    #print("EVENT SELECTION LOSSES FOR BKG CENTRAL SOURCE RUN")
    #EventSelectionLosses(Bdf,Bdf_trig)
    PlotDemo(Bdf,Bdf_trig)


