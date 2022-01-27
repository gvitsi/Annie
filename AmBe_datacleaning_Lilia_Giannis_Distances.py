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

#plt.rcParams['font.family'] = "Times"
pylab.rcParams['font.size'] = 12
pylab.rcParams['figure.figsize'] = (10, 7.6)

z=0
SIGNAL_DIR = "../Data/Calibration_2021/Signal/{}/".format(z)
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
    return df_trig_cleanWindow

#def PlotDemo(Sdf,Bdf,Sdf_trig,Bdf_trig):
def PlotDemo(Sdf,Sdf_trig):
    
    
 
    #print("EVENT SELECTION LOSSES FOR BKG CENTRAL SOURCE RUN")
    #EventSelectionLosses(Bdf,Bdf_trig)
    
    Sdf['label'] = '1'
    print(Sdf.head())
    print("Sdf.shape: ", Sdf.shape)
    print("All columns are: ", Sdf.columns.values.tolist())
    #Bdf.to_csv("vars_DNN_Signal.csv",  index=False,float_format = '%.3f')
    #print(type(Sdf.hitDetID))

    
    #---- My Plots:
    Sdf_prompt=Sdf.loc[Sdf['clusterTime']<2000].reset_index(drop=True)	#prompt events
    Sdf_del=Sdf.loc[Sdf['clusterTime']>=2000].reset_index(drop=True) #delayed events
    '''
    plt.subplots_adjust(left=0.15, right=0.92, bottom=0.14, top=0.92, hspace=0.05, wspace=0.5)
    plt.hist(Sdf_prompt['clusterTime'],bins=100,range=(0,2000))
    plt.title("Prompt window Tank cluster times - no cuts")
    plt.xlabel("Cluster time [ns]")
    plt.savefig(f"plots_AmBe/{z}/time_prompt.png")
    plt.show()   

    
    plt.subplots_adjust(left=0.15, right=0.92, bottom=0.14, top=0.92, hspace=0.05, wspace=0.5)
    plt.hist(Sdf_del['clusterTime'],bins=100,range=(2000,70000))
    plt.title("Delayed window Tank cluster times - no cuts")
    plt.xlabel("Cluster time [ns]")
    plt.savefig(f"plots_AmBe/{z}/time_del.png")
    plt.show()
    '''
    
    #--- CB to cluster Time:   
    labels = {'title': 'Charge balance parameters in time window \n (AmBe data, $0<t_{c}<=60 \, \mu s$)\n' f"for position Z={z}cm", 
             'xlabel': 'Cluster time (ns)', 'ylabel': 'Charge balance'} #total events
    plt.subplots_adjust(left=0.19, right=0.94, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,60000],'yrange':[0,1]}
    abp.Make2DHist(Sdf,'clusterTime','clusterChargeBalance',labels,ranges)
    plt.savefig(f"plots_AmBe/{z}/CB_time_total.png")
    plt.show()
   
    labels = {'title': 'Charge balance parameters in time window \n (AmBe data, $t_{c}>=2 \, \mu s$)\n' f"for position Z={z}cm", 
             'xlabel': 'Cluster time (ns)', 'ylabel': 'Charge balance'} #delayed events
    plt.subplots_adjust(left=0.17, right=0.94, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[2000,60000],'yrange':[0,1]}
    abp.Make2DHist(Sdf_del,'clusterTime','clusterChargeBalance',labels,ranges)
    plt.savefig(f"plots_AmBe/{z}/CB_time_del.png")
    plt.show()
    
    labels = {'title': 'Charge balance parameters in time window \n (AmBe data, $t_{c}<2 \, \mu s$)\n' f"for position Z={z}cm",
              'xlabel': 'Cluster time (ns)', 'ylabel': 'Charge balance'} #prompt events
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,2000],'yrange':[0,1]}
    plt.subplots_adjust(left=0.19, right=0.94, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    abp.Make2DHist(Sdf_prompt,'clusterTime','clusterChargeBalance',labels,ranges)
    plt.savefig(f"plots_AmBe/{z}/CB_time_prompt.png")
    plt.show()

    #--- CB to clusterPE: 
    plt.subplots_adjust(left=0.15, right=0.97, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) #total events
    labels = {'title': 'Charge balance parameters in time window \n (AmBe data, $0<t_{c}<=60 \, \mu s$)\n ' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Charge balance'}
    ranges = {'xbins': 200, 'ybins':100, 'xrange':[0,300],'yrange':[0,1]}
    abp.Make2DHist(Sdf,'clusterPE','clusterChargeBalance',labels,ranges)
    plt.savefig(f"plots_AmBe/{z}/CB_PE_total.png")
    plt.show()
    
        
    plt.subplots_adjust(left=0.15, right=0.97, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Charge balance parameters in time window \n (AmBe data, $t_{c}>=2 \, \mu s$)\n ' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Charge balance'}	#delayed events
    ranges = {'xbins': 200, 'ybins':100, 'xrange':[0,300],'yrange':[0,1]}
    abp.Make2DHist(Sdf_del,'clusterPE','clusterChargeBalance',labels,ranges)
    plt.savefig(f"plots_AmBe/{z}/CB_PE_del.png")
    plt.show()
    plt.subplots_adjust(left=0.15, right=0.97, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Charge balance parameters in time window \n (AmBe data, $t_{c}<2 \, \mu s$)\n ' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Charge balance'}	#prompt events
    ranges = {'xbins': 200, 'ybins':100, 'xrange':[0,300],'yrange':[0,1]}
    abp.Make2DHist(Sdf_prompt,'clusterPE','clusterChargeBalance',labels,ranges)
    plt.savefig(f"plots_AmBe/{z}/CB_PE_prompt.png")
    plt.show()
#--- CB>=0.9 
    #splitting to CB categories:
    
    Sdf_prompt_highCB = Sdf_prompt.loc[Sdf_prompt['clusterChargeBalance']>=0.9].reset_index(drop=True) 
    Sdf_del_highCB = Sdf_del.loc[Sdf_del['clusterChargeBalance']>=0.9].reset_index(drop=True)
    plt.subplots_adjust(left=0.15, right=1, bottom=0.15, top=0.8, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (AmBe data, $t_{c}<2 \, \mu s$) CB>=0.9\n' f"for position Z={z}cm", 'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 280, 'ybins':200, 'xrange':[0,140],'yrange':[0,60]}
    abp.Make2DHist(Sdf_prompt_highCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig(f"plots_AmBe/{z}/PE_maxPE_prompt_highCB.png")
    plt.show()
    
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (AmBe data, $t_{c}>=2 \, \mu s$) CB>=0.9 \n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,50],'yrange':[0,50]}
    plt.subplots_adjust(left=0.12, right=1, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2)
    abp.Make2DHist(Sdf_del_highCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe/{}/TotPE_maxPEclust_high.png".format(z))
    plt.show()

    #PE = np.hstack(Sdf_del_highCB['hitPE'])
    #ID = np.hstack(Sdf_del_highCB['hitDetID'])
    #T = np.hstack(Sdf_del_highCB['hitT'])
    #maxPE_highCB = max(np.hstack(Sdf_prompt_highCB.hitPE))
    #print("maxPE_highCB ",maxPE_highCB," clusterMaxPE ",Sdf_prompt_highCB.clusterMaxPE)

    highCB_PE_prompt = np.hstack(Sdf_prompt_highCB.hitPE)
    highCB_DetID_prompt = np.hstack(Sdf_prompt_highCB.hitDetID)
    highCB_PE_del = np.hstack(Sdf_del_highCB.hitPE)
    highCB_DetID_del = np.hstack(Sdf_del_highCB.hitDetID)
    #labels = {'title': 'PE distribution for all hits in clusters, CB>=0.9)\n' f"for position Z={z}cm" ,'xlabel': 'Tube ID', 'ylabel': 'PE'}
    #ranges = {'xbins': 200, 'ybins':100, 'xrange':[0,140],'yrange':[0,100]}
    plt.hist2d(highCB_DetID_prompt,highCB_PE_prompt, bins=(200,136), range=[[330,465],[0,50]],cmap=plt.cm.turbo)
    plt.title("Prompt PE distribution for all hits in clusters,\n CB>=0.9 " f"for position Z={z}cm")
    plt.colorbar()
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/TubeID_PE_prompt_highCB.png")
    plt.show()
    
    plt.hist2d(highCB_DetID_del,highCB_PE_del, bins=(200,136), range=[[330,465],[0,20]],cmap=plt.cm.turbo)
    plt.title("Delayed PE distribution for all hits in clusters,\n CB>=0.9 " f"for position Z={z}cm")
    plt.colorbar()
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/TubeID_PE_del_highCB.png")
    plt.show()
   
    #--- 0.6<=CB<0.9
    Sdf_prompt_upperCB = Sdf_prompt.loc[(Sdf_prompt['clusterChargeBalance']<0.9) & (Sdf_prompt['clusterChargeBalance']>=0.6)].reset_index(drop=True)
    Sdf_del_upperCB = Sdf_del.loc[(Sdf_del['clusterChargeBalance']<0.9) & (Sdf_prompt['clusterChargeBalance']>=0.6)].reset_index(drop=True)
    plt.subplots_adjust(left=0.15, right=1, bottom=0.15, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (AmBe data, $t_{c}<2 \, \mu s$) 0.6<=CB<0.9\n ' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 280, 'ybins':100, 'xrange':[0,140],'yrange':[0,60]}
    abp.Make2DHist(Sdf_prompt_upperCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig(f"plots_AmBe/{z}/PE_maxPE_prompt_upperCB.png")
    plt.show()
    
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (AmBe data, $t_{c}>=2 \, \mu s$)  0.6<=CB<0.9\n ' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,30],'yrange':[0,50]}
    plt.subplots_adjust(left=0.12, right=1, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2)
    abp.Make2DHist(Sdf_del_upperCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe/{}/TotPE_maxPEclust_upper.png".format(z))
    plt.show()


    upperCB_PE_prompt = np.hstack(Sdf_prompt_upperCB.hitPE)
    upperCB_DetID_prompt = np.hstack(Sdf_prompt_upperCB.hitDetID)
    upperCB_PE_del = np.hstack(Sdf_del_upperCB.hitPE)
    upperCB_DetID_del = np.hstack(Sdf_del_upperCB.hitDetID)
    plt.hist2d(upperCB_DetID_prompt,upperCB_PE_prompt, bins=(200,136), range=[[330,465],[0,50]],cmap=plt.cm.turbo)
    plt.title("Prompt PE distribution for all hits in clusters, \n 0.6=<CB<0.9 " f"for position Z={z}cm")
    plt.colorbar()
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/TubeID_PE_prompt_upperCB.png")
    plt.show()
    
    plt.hist2d(upperCB_DetID_del,upperCB_PE_del, bins=(200,136), range=[[330,465],[0,20]],cmap=plt.cm.turbo)
    plt.title("Delayed PE distribution for all hits in clusters, \n 0.6=<CB<0.9 " f"for position Z={z}cm")
    plt.colorbar()
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/TubeID_PE_del_upperCB.png")
    plt.show()


    #--- 0.4<=CB<0.6
    Sdf_prompt_midCB = Sdf_prompt.loc[(Sdf_prompt['clusterChargeBalance']<0.6) & (Sdf_prompt['clusterChargeBalance']>=0.4)].reset_index(drop=True)
    Sdf_del_midCB = Sdf_del.loc[(Sdf_del['clusterChargeBalance']<0.6) & (Sdf_prompt['clusterChargeBalance']>=0.4)].reset_index(drop=True)
    plt.subplots_adjust(left=0.15, right=1, bottom=0.15, top=0.8, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (AmBe data, $t_{c}<2 \mu s$) 0.4<=CB<0.6\n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 280, 'ybins':200, 'xrange':[0,140],'yrange':[0,60]}
    abp.Make2DHist(Sdf_prompt_midCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig(f"plots_AmBe/{z}/PE_maxPE_prompt_midCB.png")
    plt.show()
    
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (AmBe data, $t_{c}>=2 \, \mu s$) 0.4<=CB<0.6\n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,30],'yrange':[0,50]}
    plt.subplots_adjust(left=0.12, right=1, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2)
    abp.Make2DHist(Sdf_del_midCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe/{}/TotPE_maxPEclust_mid.png".format(z))
    plt.show()

     
    midCB_PE_prompt = np.hstack(Sdf_prompt_midCB.hitPE)
    midCB_DetID_prompt = np.hstack(Sdf_prompt_midCB.hitDetID)
    midCB_PE_del = np.hstack(Sdf_prompt_midCB.hitPE)
    midCB_DetID_del = np.hstack(Sdf_prompt_midCB.hitDetID)
    plt.hist2d(midCB_DetID_prompt,midCB_PE_prompt, bins=(200,136), range=[[330,465],[0,50]],cmap=plt.cm.turbo)
    plt.title("Promot PE distribution for all hits in clusters,\n 0.4=<CB<0.6 "f"for position Z={z}cm")
    plt.colorbar()
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/TubeID_PE_prompt_midCB.png")
    plt.show()
    
    plt.hist2d(midCB_DetID_del,midCB_PE_del, bins=(200,136), range=[[330,465],[0,50]],cmap=plt.cm.turbo)
    plt.title("Delayed PE distribution for all hits in clusters,\n 0.4=<CB<0.6 "f"for position Z={z}cm")
    plt.colorbar()
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/TubeID_PE_del_midCB.png")
    plt.show()

    #--- CB<0.4
    Sdf_prompt_lowCB = Sdf_prompt.loc[Sdf_prompt['clusterChargeBalance']<0.4].reset_index(drop=True)
    Sdf_del_lowCB = Sdf_del.loc[Sdf_del['clusterChargeBalance']<0.4].reset_index(drop=True)
    plt.subplots_adjust(left=0.15, right=1, bottom=0.15, top=0.8, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (AmBe data, $t_{c}<2 \, \mu s$) CB<0.4 \n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,140],'yrange':[0,60]}
    abp.Make2DHist(Sdf_prompt_lowCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig(f"plots_AmBe/{z}/PE_maxPE_prompt_lowCB.png")
    plt.show()
    
    labels = {'title': 'Total PE vs Maximum PE in Cluster for \n (AmBe data, $t_{c}>=2 \, \mu s$) CB<0.4\n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Maximum PE in Cluster'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,50],'yrange':[0,50]}
    plt.subplots_adjust(left=0.12, right=1, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2)
    abp.Make2DHist(Sdf_del_lowCB,'clusterPE','clusterMaxPE',labels,ranges)
    plt.savefig("plots_AmBe/{}/TotPE_maxPEclust_low.png".format(z))
    plt.show()
      
    lowCB_PE_prompt = np.hstack(Sdf_prompt_lowCB.hitPE)
    lowCB_DetID_prompt = np.hstack(Sdf_prompt_lowCB.hitDetID)
    lowCB_PE_del = np.hstack(Sdf_prompt_lowCB.hitPE)
    lowCB_DetID_del = np.hstack(Sdf_prompt_lowCB.hitDetID)
    
    plt.hist2d(lowCB_DetID_prompt,lowCB_PE_prompt, bins=(200,136), range=[[330,465],[0,50]],cmap=plt.cm.turbo)
    plt.title("Prompt PE distribution for all hits in clusters,\n CB<=0.4 " f"for position Z={z}cm")
    plt.colorbar()
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/TubeID_PE_prompt_lowCB.png")
    plt.show()
    
    plt.hist2d(lowCB_DetID_del,lowCB_PE_del, bins=(200,136), range=[[330,465],[0,50]],cmap=plt.cm.turbo)
    plt.title("Delayed PE distribution for all hits in clusters,\n CB<=0.4 " f"for position Z={z}cm")
    plt.colorbar()
    plt.xlabel("Tube ID")
    plt.ylabel("PE")
    plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/TubeID_PE_del_lowCB.png")
    plt.show()


    #----- 1D Histogram for all CB areas -----#:
    
    #------prompt charts------#
    #--- CB>=0.9 
    plt.hist(Sdf_prompt_highCB['clusterPE'],bins=100,range=(0,500), histtype='step', density=True )
    plt.title("Prompt Cluster PE, CB>=0.9\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/PromptHistClusterPE_highCB.png")
    plt.show()
    
    #--- 0.6<=CB<0.9
    plt.hist(Sdf_prompt_upperCB['clusterPE'],bins=100,range=(0,500), histtype='step', density=True )
    plt.title("Prompt Cluster PE, 0.6<=CB<0.9\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/PromptHistClusterPE_upperCB.png")
    plt.show()
    
    #--- 0.4<=CB<0.6
    plt.hist(Sdf_prompt_midCB['clusterPE'],bins=100,range=(0,500), histtype='step', density=True )
    plt.title("Prompt Cluster PE, 0.4<=CB<0.6\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/PromptHistClusterPE_midCB.png")
    plt.show()
    
    #--- CB<0.4
    plt.hist(Sdf_prompt_lowCB['clusterPE'],bins=100,range=(0,500), histtype='step', density=True )
    plt.title("Prompt Cluster PE, CB<0.4\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/PromptHistClusterPE_lowCB.png")
    plt.show()
    
    
    #-----delayed charts------#
    #--- CB>=0.9 
    plt.hist(Sdf_del_highCB['clusterPE'],bins=100,range=(2000,20000), histtype='step', density=True )
    plt.title("Delayed Cluster PE, CB>=0.9\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/DelHistClusterPE_highCB.png")
    plt.show()
    
    #--- 0.6<=CB<0.9
    plt.hist(Sdf_del_upperCB['clusterPE'],bins=100,range=(2000,20000), histtype='step', density=True )
    plt.title("Delayed Cluster PE, 0.6<=CB<0.9\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/DelHistClusterPE_upperCB.png")
    plt.show()
    
    #--- 0.4<=CB<0.6
    plt.hist(Sdf_del_midCB['clusterPE'],bins=100,range=(2000,20000), histtype='step', density=True )
    plt.title("Delayed Cluster PE, 0.4<=CB<0.6\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/DelHistClusterPE_midCB.png")
    plt.show()
    
    #--- CB<0.4
    plt.hist(Sdf_del_lowCB['clusterPE'],bins=100,range=(2000,20000), histtype='step', density=True )
    plt.title("Delayed Cluster PE, CB<0.4\n" f"for position Z={z}cm")
    plt.xlabel("Cluster PE")
    #plt.subplots_adjust(left=0.14, right=0.94, bottom=0.14, top=0.86, hspace=0.05, wspace=0.2) 
    plt.savefig(f"plots_AmBe/{z}/DelHistClusterPE_lowCB.png")
    plt.show()
    
    
    #----- Plots only for 0.4 < CB < 0.6 -----#:
    plt.subplots_adjust(left=0.15, right=0.98, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Charge balance parameters in time window \n (Bkgd data, $t_{c}=>2 \, \mu s$, 0.4 < CB < 0.6)\n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Charge balance'} #delayed events
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,300],'yrange':[0,1]}
    abp.Make2DHist(Sdf_del_midCB,'clusterPE','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe/{}/CB_PE_del_midCB.png".format(z))
    plt.show()
    
    plt.subplots_adjust(left=0.15, right=0.98, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Charge balance parameters in time window \n (Bkgd data, $t_{c}<2 \, \mu s$, 0.4 < CB < 0.6)\n' f"for position Z={z}cm",
              'xlabel': 'Cluster PE', 'ylabel': 'Charge balance'} #prompt events
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,300],'yrange':[0,1]}
    abp.Make2DHist(Sdf_prompt_midCB,'clusterPE','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe/{}/CB_PE_prompt_midCB.png".format(z))
    plt.show()
    
    plt.subplots_adjust(left=0.17, right=0.94, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Charge balance parameters in time window \n (Bkgd data, $t_{c}>=2 \, \mu s$, 0.4 < CB < 0.6)\n' f"for position Z={z}cm", 
             'xlabel': 'Cluster time (ns)', 'ylabel': 'Charge balance'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[2000,30000],'yrange':[0,1]}
    abp.Make2DHist(Sdf_del_midCB,'clusterTime','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe/{}/CB_time_del_midCB.png".format(z))
    plt.show()
#    
    plt.subplots_adjust(left=0.19, right=0.94, bottom=0.14, top=0.80, hspace=0.05, wspace=0.2) 
    labels = {'title': 'Charge balance parameters in time window \n (Bkgd data, $t_{c}<2 \, \mu s$, 0.4 < CB < 0.6)\n' f"for position Z={z}cm",
              'xlabel': 'Cluster time (ns)', 'ylabel': 'Charge balance'}
    ranges = {'xbins': 200, 'ybins':200, 'xrange':[0,2000],'yrange':[0,1]}
    abp.Make2DHist(Sdf_prompt_midCB,'clusterTime','clusterChargeBalance',labels,ranges)
    plt.savefig("plots_AmBe/{}/CB_time_prompt_midCB.png".format(z))
    plt.show()
    
    
    
    
    
    
    
'''
#########
    Sdf_prompt_noCB = Sdf.loc[Sdf['clusterTime']<2000].reset_index(drop=True)
    Sdf_prompt = Sdf_prompt_noCB.loc[Sdf_prompt_noCB['clusterChargeBalance']<0.9].reset_index(drop=True) ##for flasherscheck clusters with CB>0.9
#    plt.hist(Sdf_prompt['clusterTime'],bins=100,range=(0,2000))
#    plt.title("Prompt window Tank cluster times")
#    plt.xlabel("Cluster time [ns]")
#    plt.show()
    print("TOTAL PROMPT TANK CLUSTERS, NO CB: " + str(len(Sdf_prompt_noCB)))
    print("TOTAL PROMPT TANK CLUSTERS: " + str(len(Sdf_prompt)))
    print("TOTAL PROMPT MRD CLUSTERS: " + str(len(Sdf_mrd)))
   
    labels = {'title': 'Charge balance parameters in time window \n (Beam data, $t_{c}<2 \, \mu s$)', 
            'xlabel': 'Cluster time (ns)', 'ylabel': 'Charge balance'}
    ranges = {'xbins': 40, 'ybins':40, 'xrange':[0,2000],'yrange':[0,1]}
    #abp.MakeHexJointPlot(Sdf,'clusterPE','clusterChargeBalance',labels,ranges)
    abp.Make2DHist(Sdf_prompt_noCB,'clusterTime','clusterChargeBalance',labels,ranges)
    plt.show()

    labels = {'title': 'Tank PMT hit cluster count as a function of time \n (Beam data, $t_{c}<2 \, \mu s$)', 
            'xlabel': 'Cluster time (ns)', 'ylabel': 'Cluster PE'}
    ranges = {'xbins': 40, 'ybins':25, 'xrange':[0,2000],'yrange':[0,5000]}
    #abp.MakeHexJointPlot(Sdf,'clusterPE','clusterChargeBalance',labels,ranges)
    abp.Make2DHist(Sdf_prompt_noCB,'clusterTime','clusterPE',labels,ranges)
    plt.show()


    #Get largest cluster in each acquisition in prompt window
    Sdf_maxPE = es.MaxPEClusters(Sdf_prompt)
    print("TOTAL HIGHEST PE PROMPT CLUSTERS: " + str(len(Sdf_maxPE)))
    Sdf_mrd_maxhit = es.MaxHitClusters(Sdf_mrd)
    print("TOTAL MOST PADDLE MRD CLUSTERS: " + str(len(Sdf_mrd_maxhit)))
'''
'''
    #Now, get the index number for clusterTime pairs in the same triggers 
    TankIndices, MRDIndices = es.MatchingEventTimes(Sdf_maxPE,Sdf_mrd_maxhit)
    TankTimes = Sdf_maxPE["clusterTime"].values[TankIndices]
    MRDTimes = Sdf_mrd_maxhit["clusterTime"].values[MRDIndices]
    Pairs_HaveVeto = Sdf_mrd_maxhit.loc[(Sdf_mrd_maxhit["vetoHit"].values[MRDIndices]==1)]
    print("NUM OF MRD CLUSTERS IN TRIG WITH A TANK CLUSTER: " + str(len(MRDTimes)))
    print("NUM OF MRD CLUSTERS WITH VETO IN SUBSET: " + str(len(Pairs_HaveVeto)))
    plt.scatter(TankTimes,MRDTimes,marker='o',s=15,color='blue',alpha=0.7)
    plt.title("Tank and MRD cluster times in prompt window \n (Largest PE tank clusters, largest paddle count MRD clusters)")
    plt.xlabel("Tank Cluster time [ns]")
    plt.ylabel("MRD Cluster time [ns]")
    plt.show()

    plt.hist(MRDTimes - TankTimes, bins = 160, color='blue', alpha=0.7)
    plt.axvline(x=700,color='black',linewidth=6)
    plt.axvline(x=800,color='black',linewidth=6)
    plt.title("Difference in MRD and Tank cluster times in acquisitions \n (Largest PE tank clusters, largest paddle count MRD clusters)")
    plt.xlabel("MRD cluster time - Tank cluster time [ns]")
    plt.show()


    print("CLUSTER COUNT IN EVENTS BEFORE 2 US: " + str(len(Sdf_ClustersInPromptCandidates.loc[Sdf_ClustersInPromptCandidates["clusterTime"]<2000].values)))
    Sdf_ValidDelayedClusters = Sdf_ClustersInPromptCandidates.loc[Sdf_ClustersInPromptCandidates['clusterTime']>12000].reset_index(drop=True)
    Sdf_ValidDelayedClustersCB = Sdf_ClustersInPromptCandidates.loc[Sdf_ClustersInPromptCandidates['clusterChargeBalance']<0.4].reset_index(drop=True)
    print("CLUSTER COUNT IN EVENTS WITH PMT/MRD ACTIVITY PAST 12 US: " + str(len(Sdf_ValidDelayedClusters)))

    plt.hist(Sdf.loc[Sdf["clusterTime"]>12000,"clusterTime"],bins=20,range=(12000,65000),label='No PMT/MRD pairing in prompt',alpha=0.8)
    plt.hist(Sdf_ValidDelayedClusters["clusterTime"], bins=20, range=(12000,65000),label='PMT/MRD pair required in prompt',alpha=0.8)
    plt.hist(Sdf_ValidDelayedClustersCB["clusterTime"], bins=20, range=(12000,65000),label=' + CB < 0.4',alpha=0.8)
    plt.title("Delayed cluster times in beam runs")
    plt.ylabel("Number of clusters")
    plt.xlabel("Cluster time [ns]")
    leg = plt.legend(loc=1,fontsize=24)
    leg.set_frame_on(True)
    leg.draw_frame(True)
    plt.show()

'''
if __name__=='__main__':
    slist = glob.glob(SIGNAL_DIR+"*.ntuple.root")
    #blist = glob.glob(BKG_DIR+"*.ntuple.root")

    livetime_estimate = es.EstimateLivetime(slist)
    print("SIGNAL LIVETIME ESTIMATE IN SECONDS IS: " + str(livetime_estimate))
    #livetime_estimate = es.EstimateLivetime(blist)
    #print("BKG LIVETIME ESTIMATE IN SECONDS IS: " + str(livetime_estimate))

    #mybranches = ['eventNumber','eventTimeTank','clusterTime','SiPMhitT','SiPMhitQ','SiPMhitAmplitude','clusterChargeBalance','clusterPE','SiPM1NPulses','SiPM2NPulses','SiPMNum','clusterHits']
    mybranches = ['eventNumber','eventTimeTank','clusterTime','hitT','hitQ','hitPE','hitDetID','clusterChargeBalance','clusterPE','clusterMaxPE' ,'clusterHits', 'SiPMhitT','SiPMhitQ','SiPMhitAmplitude','SiPM1NPulses','SiPM2NPulses','SiPMNum']
    #mybranches = ['clusterTime','hitT','hitQ','hitPE','hitDetID','clusterChargeBalance','clusterPE','clusterMaxPE','clusterHits']

    SProcessor = rp.ROOTProcessor(treename="phaseIITankClusterTree")
    for f1 in slist:
        SProcessor.addROOTFile(f1,branches_to_get=mybranches)
    Sdata = SProcessor.getProcessedData()
    Sdf = pd.DataFrame(Sdata)

    #BProcessor = rp.ROOTProcessor(treename="phaseIITankClusterTree")
    #for f1 in blist:
        #BProcessor.addROOTFile(f1,branches_to_get=mybranches)
    #Bdata = BProcessor.getProcessedData()
    #Bdf = pd.DataFrame(Bdata)

    SProcessor = rp.ROOTProcessor(treename="phaseIITriggerTree")
    for f1 in slist:
        SProcessor.addROOTFile(f1,branches_to_get=mybranches)
    Sdata = SProcessor.getProcessedData()
    Sdf_trig = pd.DataFrame(Sdata)

    #BProcessor = rp.ROOTProcessor(treename="phaseIITriggerTree")
    #for f1 in blist:
        #BProcessor.addROOTFile(f1,branches_to_get=mybranches)
    #Bdata = BProcessor.getProcessedData()
    #Bdf_trig = pd.DataFrame(Bdata)
    
    #print("EVENT SELECTION LOSSES FOR CENTRAL SOURCE RUN")
    #sdf_trig_cleanWindow = EventSelectionLosses(Sdf,Sdf_trig)
    
    PlotDemo(Sdf,Sdf_trig)


