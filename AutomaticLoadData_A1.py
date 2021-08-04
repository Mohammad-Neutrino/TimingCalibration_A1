from __future__ import print_function
#from ROOT import TCanvas, TGraph
#from ROOT import gROOT
import ROOT
import os
import matplotlib.pyplot as plt
from scipy import signal,stats
from scipy.stats import norm
from scipy.fftpack import rfft, irfft, fftfreq
from scipy.signal import iirfilter,lfilter,butter,hilbert
from scipy import signal
from scipy.fftpack import fft
from scipy import optimize
from scipy.misc import derivative
from scipy.interpolate import interp1d, Akima1DInterpolator
import numpy as np
import sys
#sys.path.append('../')
#import global_constants
import math
from math import sin
from array import array
from pynverse import inversefunc
#from TimingCalibration import partial_derivative
#from VoltageCalibration import BlockCorrector, MeanEachBlock
import matplotlib
import scipy
import itertools
#from CalPulser import loadvalues
import warnings
#import requests
#from bs4 import BeautifulSoup
from ctypes import c_int as ll
import urllib.request as request
#import urllib.error as error
#import urllib2 #for python2
ROOT.gSystem.Load("$ARA_UTIL_INSTALL_DIR/lib/libAraEvent.so")
warnings.simplefilter(action='ignore', category=FutureWarning)



def MeanEachBlock(v,channel):
    new_v = np.zeros(len(v))
    if(int(channel) not in [0,1,8,9,16,17,24,25]):
        size = 32
    else:
        size = 64
    for k in range(0,int(len(new_v)/size)):
        temp_v = v[k*size:k*size+size]
        temp_v[1::2]=temp_v[1::2]-np.mean(temp_v[1::2])
        temp_v[0::2]=temp_v[0::2]-np.mean(temp_v[0::2])
        new_v[k*size:k*size+size] = temp_v
        #new_v[k*size:k*size+size] = v[k*size:k*size+size]-np.mean(v[k*size:k*size+size])
    return(new_v)

def RemoveFirstBlock(time,volt,b_number):
    blocks_to_cut=1 #at least 1
    if(b_number%2==0):
        volt=volt[blocks_to_cut*64:]#1280
        time=time[blocks_to_cut*64:]
        b_number = (b_number + blocks_to_cut)%512
    else:
        volt=volt[(blocks_to_cut+0)*64:]#1280
        time=time[(blocks_to_cut+0)*64:]
        b_number = (b_number + blocks_to_cut+0)%512
    return(time,volt,b_number)

def SameLengthWaveforms(time,volt,length):
    n = len(volt)
    if(n>length):
        #print('too long')
        volt=volt[:length]
        time=time[:length]
    if(n<length):
        #print('too short')

        volt=np.append(volt,np.zeros(length-n))

        time=np.append(time,np.linspace(time[-1],time[-1]+(time[1]-time[0])*(-1+length-n),(length-n)))


    return(time,volt)

def fixlength(t,length):#calibrated time will always start on even block

    spacing = t[1:129]-t[:128]
    t_new = np.zeros(length)

    for b in range(1,length):
        t_new[b]=t_new[b-1]+spacing[(b-1)%128]

    return(t_new)

def PedestalFix(v,channel,block_num,ped_values): # this is for using an ARA generated ped file

    #print('ped values are', ped_values)
    my_ped_values = np.zeros(len(v))
    counter = 0
    for i in range(0,int(len(v)/64)):
        v0 =int(int(channel)/8)
        v1 =int(block_num)
        v2 =int(int(channel)%8)

        ind_test = np.where((ped_values[:,0]==v0))
        my_ped_row = np.where((ped_values[:,0]==v0) & (ped_values[:,1]==v1) &(ped_values[:,2]==v2))
        my_ped_row_vals = ped_values[my_ped_row[0],:]
        my_ped_values[counter:counter+64]=my_ped_row_vals[0,3:]
        block_num=(block_num+1)%512
        counter = counter +64
    return(v-my_ped_values)


def VoltageCorrector(v,block,p_pos,p_neg,total_samples,chi2p,chi2n):

    #test =Cubic(50.0,p_pos[50])
    new_v = np.zeros(len(v))
    new_b = block

    for t in range(0,total_samples):
        if t%2==0 :
            t_odd = t+1
        else:
            t_odd = t

        best_index= (int(new_b)*64+t_odd%64)%32768
        #print(chi2p)
        while(chi2p[best_index]>1.0):
            best_index=(best_index+2)%32768

        if(t%64==0 and t>0):
            #print('hello')
            new_b=int(new_b+1)%512
        if(v[t]>=0.0):
            val=Cubic(v[t],p_pos[best_index,:])
        if(v[t]<0.0):
            val=Cubic(v[t],p_neg[best_index,:])
        if(np.abs(val)>800):
            new_v[t]=v[t]
        else:
            new_v[t]=val

    return(new_v)

def Cubic(a,p):

    deg = len(p)
    p_p=0.0
    for pi in range(0,deg):
        p_p = p_p+p[deg-pi-1]*a**pi

    return(p_p)


def Calibrator(station,time,volt,block_number,channel,length,ped_values,kPed,kTime,kVolt):
    #MyFile = ROOT.TFile.Open('data/processed/calibration_data_full_elChan'+channel+'_run'+files+'.root')
    #block_nums0 = np.loadtxt("data/processed/block_data_elChan"+channel+"_run"+files+'.txt')
    #print(kPed,kTime,kVolt)


    if(int(channel) not in [0,1,2,3,8,9,10,11,16,17,18,19,24,25,26,27]):
        channel='18'
    #print('calibrating')
    cal_file = 'data/ARA'+station+'_cal_files'
    #cal_file = 'cal_files8'

    if(kPed==1):
        volt = PedestalFix(volt,channel,block_number,ped_values)
        volt = MeanEachBlock(volt,channel)
    """
    if(int(channel)==24 and kVolt==1):
        volt = Ch24VoltageCorrector(volt,block_number)
        volt = MeanEachBlock(volt,channel)
        return(time,volt,block_number)
    """
    #plt.figure(3)
    #plt.plot(time,volt)
    if(kTime==1):
        #print('hello! calibrating the time')
        #t_cal = np.load('/home/kahughes/ARA/'+cal_file+'/t_cal_'+str(channel)+'.npy')
        t_cal = np.load(cal_file+'/t_cal_'+str(channel)+'.npy')
        t_cal = fixlength(t_cal,length+64)
        time,volt = SameLengthWaveforms(time,volt,length+64)
        #plt.figure()
        #plt.plot(time,volt)
        if(block_number%2==0):
            #print('should be here')
            """
            if(int(channel) not in [0,1,8,9,16,17,24,25]):
                time=t_cal[:length/2]
                volt=volt[:length]
                volt=volt[1::2]
                length=length/2
            else:
            """

            time=t_cal[:length]+20.0
            volt=volt[:length]

            #plt.figure()
            #plt.plot(time,volt)
        else:
            #print('should not be here')
            """
            if(int(channel) not in [0,1,8,9,16,17,24,25]):

                time=t_cal[:length/2]
                volt=volt[32:length+32]
                volt=volt[1::2]
                length=length/2

            else:
            """

            time=t_cal[64:length+64]
            volt=volt[:length]
            block_number = block_number+1
    #else:
    #time=time+20.0
    #print(time)

    #plt.plot(time,volt)

    """
    else:
        time,volt = SameLengthWaveforms(time,volt,length+64)
        if(block_number%2==0):
            time=time[:length]
            volt=volt[:length]
        else:
            time=time[:length]
            volt=volt[64:length+64]
            block_number = block_number+1
    """
    #plt.figure()
    #plt.plot(time,volt)
    #plt.show()
    volt = volt[:length]
    time = time[:length]
    if(kVolt==1):
        if(int(channel)==24):
            volt = Ch24VoltageCorrector(volt,block_number)
            volt = MeanEachBlock(volt,channel)
        else:
            if(int(channel)%4==3):
                channel = str(int(channel)-2)
            chi2_pos = np.load(cal_file+'/chi2_pos_'+str(channel)+'.npy')
            chi2_neg = np.load(cal_file+'/chi2_neg_'+str(channel)+'.npy')
            #pedestals0 = np.load('best_pedestals/ch_'+channel+'_ped.npy')
            p_pos = np.load(cal_file+'/p_pos_'+str(channel)+'.npy')
            p_neg = np.load(cal_file+'/p_neg_'+str(channel)+'.npy')

            volt = VoltageCorrector(volt,block_number,p_pos,p_neg,length,chi2_pos,chi2_neg)
            volt = MeanEachBlock(volt,channel)


    #plt.plot(time,volt)
    #plt.show()
    #print('done')
    return(time,volt,block_number)

def RemoveBackwardsSamples(tcal,v):
    if(tcal.ndim>1):
        diffs = tcal[0,1:]-tcal[0,:-1]
        backwards_args = np.where(diffs<0)
        t_new = np.zeros([len(tcal[:,0]),len(tcal[0,:])-len(backwards_args[0])])
        #print(np.shape(t_new))
        v_new = np.zeros([len(tcal[:,0]),len(tcal[0,:])-len(backwards_args[0])])
        for z in range(0,len(tcal[:,0])):
            #print(tcal[z,1:])
            #print(tcal[z,:-1])
            diffs = tcal[z,1:]-tcal[z,:-1]
            #print(len(tcal),len(v[0]))
            #print(diffs)
            backwards_args = np.where(diffs<0)
            #print(backwards_args)
            #print('backwards samples are:',backwards_args[0])
            t_new[z,:] = np.delete(tcal[z,:],backwards_args[0])
            v_new[z,:] = np.delete(v[z,:],backwards_args[0])
            #v = np.delete(v,backwards_args[0],axis=1)
            #print(len(tcal),len(v[0]))
            #print(tcal,v)
        return(t_new,v_new)
    else:
        #print('single array!')
        diffs = tcal[1:]-tcal[:-1]
        backwards_args = np.where(diffs<0)
        #print(backwards_args[0])
        if(len(backwards_args[0])==0):
            return(tcal,v)
        else:
            #print('uh oh! Backwards samples. Removing them?')
            t_new = np.delete(tcal,backwards_args[0])
            v_new = np.delete(v,backwards_args[0])

        return(t_new,v_new)

def LoadSineWaveData(station,run,pedFile,channel,kPed,kTime,kVolt):
    
    #data_directory = "/project2/avieregg/ARA_cal_data/data/root/"
    # data_directory = "/data/user/khughes/ARA5_cal_data/"
    data_directory = "/home/mhossain/ARA_Calibration/predeployment_data/"

    if(os.path.isfile("/data/user/khughes/ARA5_cal_data/SavedCalibData/time_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy")):
        all_t = np.load("/data/user/khughes/ARA5_cal_data/SavedCalibData/time_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy")
        all_volts = np.load("/data/user/khughes/ARA5_cal_data/SavedCalibData/volts_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy")
        all_blocks = np.load("/data/user/khughes/ARA5_cal_data/SavedCalibData/blocks_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy")
        print('Testing Code August 2')
        return(all_t,all_volts,all_blocks)
    try:
        ped_values = np.loadtxt(pedFile)
    except:
        print('unknown pedestal file location. Are you pointing to the right place?')

    test = ROOT.TFile.Open(data_directory+"run000"+str(run)+"/event000"+str(run)+".root")
    calibrator = ROOT.AraEventCalibrator.Instance()
    eventTree = test.Get("eventTree")

    rawEvent = ROOT.RawAtriStationEvent()
    eventTree.SetBranchAddress("event",ROOT.AddressOf(rawEvent))
    totalEvents = eventTree.GetEntries()
    print('Total Events:', totalEvents)
    length = 384
    all_volts = np.zeros([totalEvents,length])
    all_t=np.zeros([totalEvents,length])
    all_blocks=np.zeros([totalEvents])+701
    print('all blocks from definition', all_blocks, np.shape(all_blocks))

    for i in range(0,totalEvents):#totalEvents):
        if(i%1000==0):
            print(i)
        eventTree.GetEntry(i)

        usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kNoCalib)
        gr1 = usefulEvent.getGraphFromElecChan(int(channel))
        block_number = rawEvent.blockVec[0].getBlock()
        ##print('block_number from rawEvent', block_number, np.shape(block_number))
        if(block_number%2==0): #BECAUSE WE REMOVE THE FIRST BLOCK, WANT TO IGNORE ALL THAT HAVE FIRST EVEN BLOCK
            continue
        ##print('block_number from rawEvent after ignoring even blocks', block_number, np.shape(block_number))
        t_buff = gr1.GetX()
        v_buff = gr1.GetY()
        n = gr1.GetN()
        #t_buff.SetSize(n)
        #v_buff.SetSize(n)
        t = np.frombuffer(t_buff,dtype=float,count=-1)
        v = np.frombuffer(v_buff,dtype=float,count=-1)

        #Remove first block which is corrupted
        t,v,block_number = RemoveFirstBlock(t,v,block_number)
        ##print('block_number from RemoveFirstBlock', block_number, np.shape(block_number))
        #tt,vv,bblock_number = Calibrator(station,t,v,block_number,str(channel),length,ped_values,kPed,0,0)
        #plt.scatter(tt-20,vv,color='red')
        t,v,block_number = Calibrator(station,t,v,block_number,str(channel),length,ped_values,kPed,kTime,kVolt)
        ##print('block_number from Calibrator', block_number, np.shape(block_number))
        ##print('t after calibrator', t, np.shape(t))
        ##print('v after calibrator', v, np.shape(v))
        t=t-20
        #print(block_number)
        ##print(t[0])
        #print(len(v))
        """
        plt.scatter(t,v)
        params = SineFit(t,v,0.218)
        t_up = np.linspace(t[0],t[-1],5000)
        plt.plot(t_up,SineFunc(t_up,params[0],params[1],params[2]))
        plt.show()
        """

        all_volts[i,:]=v
        all_t[i,:]=t
        all_blocks[i]=block_number
        ##print('all_volts inside loop', all_volts, np.shape(all_volts))
        ##print('all_t inside loop', all_t, np.shape(all_t))
        ##print('all blocks inside loop', all_blocks, np.shape(all_blocks))
        #print('check 0')
        gr1.Delete()
        #usefulEvent.Delete()
        #print('check 1')

    all_t = all_t[~np.all(all_volts==0,axis=1)]
    all_volts = all_volts[~np.all(all_volts == 0, axis=1)]
    all_blocks = all_blocks[all_blocks !=701]
    ##print('all_volts outside loop', all_volts, np.shape(all_volts))
    ##print('all_t outside loop', all_t, np.shape(all_t))
    ##print('all blocks outside the loop', all_blocks, np.shape(all_blocks))
    #print(all_t)
    if(kVolt==1):
        all_t,all_volts= RemoveBackwardsSamples(all_t,all_volts)

    #print(t,all_volts)
    np.save("/home/mhossain/ARA_Calibration/data/SavedCalibData/time_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_t)
    np.save("/home/mhossain/ARA_Calibration/data/SavedCalibData/volts_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_volts)
    np.save("/home/mhossain/ARA_Calibration/data/SavedCalibData/blocks_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_blocks)

    return(all_t,all_volts,all_blocks)

def LoadARACalPulsers(run,channel,kCalPulser,kPed,kTime,kVolt):
    #t_cal,volt=LoadFilteredData('5','3416','0602','2018',int(channel),total_samples,1,1,0,0)
    length = 1792
    station= '5'
    #run = '3416'
    date_list=np.load('data/date_list.npy',allow_pickle='TRUE').item()
    #print(date_list)
    t_file = "/data/user/khughes/ARA5_cal_data/SavedCalibData/time_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy"
    v_file = "/data/user/khughes/ARA5_cal_data/SavedCalibData/volts_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy"
    #b_file = "SavedCalibData/blocks_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy"
    if(os.path.isfile(t_file) and os.path.isfile(v_file)):
        t=np.load(t_file)
        v=np.load(v_file)
        #b = np.load(b_file)
        return(t,v)
    else:
        ped_values,address_ped = LoadNewPed(run)
        #test = ROOT.TWebFile.Open("http://icecube:skua@sundog.uchicago.edu/convey/data/exp/ARA/"+year+"/filtered/L0/ARA05/"+date+"/run"+run+"/event"+run+".root")
        test = ROOT.TWebFile.Open(date_list[int(run)])
        print(date_list[int(run)])
        calibrator = ROOT.AraEventCalibrator.Instance()
        eventTree = test.Get("eventTree")

        rawEvent = ROOT.RawAtriStationEvent()
        eventTree.SetBranchAddress("event",ROOT.AddressOf(rawEvent))
        totalEvents = eventTree.GetEntries()
        print('total events:', totalEvents)


        t_univ = np.arange(21.0,length*0.3125-1.0,0.03125)

        all_volts = np.zeros([totalEvents,len(t_univ)])
        all_t=np.zeros([totalEvents,len(t_univ)])
        all_blocks=np.zeros([totalEvents])+701

        print('here is t_univ:', t_univ)
        for i in range(0,totalEvents):#totalEvents):
            #print(i)
            eventTree.GetEntry(i)
            #print('success')
            if(rawEvent.isCalpulserEvent()==0 and kCalPulser==1): #if not a cal pulser and we want cal pulsers, go to next event
                continue
            if(rawEvent.isCalpulserEvent()==1 and kCalPulser==0):
                continue
            #print(i)
            usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kNoCalib)
            gr1 = usefulEvent.getGraphFromElecChan(channel)

            t_buff = gr1.GetX()
            v_buff = gr1.GetY()
            n = gr1.GetN()
            #t_buff.SetSize(n)
            #v_buff.SetSize(n)
            t = np.frombuffer(t_buff,dtype=float,count=-1)
            v = np.frombuffer(v_buff,dtype=float,count=-1)
            #print(t,v)
            #print(len(t),len(v))
            plt.plot(t,v)
            block_number = rawEvent.blockVec[0].getBlock()

            #Remove first block which is corrupted
            t,v,block_number = RemoveFirstBlock(t,v,block_number)


            t,v,block_number = Calibrator(station,t,v,block_number,str(channel),length,ped_values,kPed,kTime,kVolt)
            plt.plot(t,v)
            plt.savefig('test_'+str(i)+'.pdf')
            plt.clf()
            #print(t[0])
            #print(t[0])
            if(int(channel)==24):
                t=t-0.3125

            #print('test',t,v)
            t,v= RemoveBackwardsSamples(t,v)
            f = Akima1DInterpolator(t,v)
            v = f(t_univ)
            #print(v)
            all_volts[i,:]=v
            #all_t[i,:]=t
            all_blocks[i]=block_number
            #print('test2')
            gr1.Delete()
            #usefulEvent.Delete()
            #print('test3')
            #print(v)
            #print(t_univ)
            #plt.plot(t_univ,v)
            #plt.show()

        #all_t = all_t[~np.all(all_volts==0,axis=1)]
        all_volts = all_volts[~np.all(all_volts == 0, axis=1)]
        all_blocks = all_blocks[all_blocks !=701]
        t_univ = t_univ-21.0

        np.save("data/SavedCalibData/time_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy",t_univ)
        np.save("data/SavedCalibData/volts_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_volts)
        np.save("data/SavedCalibData/blocks_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_blocks)


        return(t_univ,all_volts)

def LoadFilteredData(station,run,date,year,channel,length,kCalPulser,kPed,kTime,kVolt):
    t_file = "SavedCalibData/time_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy"
    v_file = "SavedCalibData/volts_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy"
    b_file = "SavedCalibData/blocks_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy"
    if(os.path.isfile(t_file) and os.path.isfile(v_file)):
        t=np.load(t_file)
        v=np.load(v_file)
        b = np.load(b_file)
        return(t,v)
    else:
        ped_values,address_ped = LoadNewPed(run)
        test = ROOT.TWebFile.Open("http://icecube:skua@sundog.uchicago.edu/convey/data/exp/ARA/"+year+"/filtered/L0/ARA05/"+date+"/run"+run+"/event"+run+".root")
        calibrator = ROOT.AraEventCalibrator.Instance()
        eventTree = test.Get("eventTree")

        rawEvent = ROOT.RawAtriStationEvent()
        eventTree.SetBranchAddress("event",ROOT.AddressOf(rawEvent))
        totalEvents = eventTree.GetEntries()
        print('total events:', totalEvents)


        t_univ = np.arange(21.0,length*0.3125-1.0,0.03125)

        all_volts = np.zeros([totalEvents,len(t_univ)])
        all_t=np.zeros([totalEvents,len(t_univ)])
        all_blocks=np.zeros([totalEvents])+701

        print('here is t_univ:', t_univ)
        for i in range(0,totalEvents):#totalEvents):

            eventTree.GetEntry(i)
            if(rawEvent.isCalpulserEvent()==0 and kCalPulser==1): #if not a cal pulser and we want cal pulsers, go to next event
                continue
            if(rawEvent.isCalpulserEvent()==1 and kCalPulser==0):
                continue
            print(i)
            usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kNoCalib)
            gr1 = usefulEvent.getGraphFromElecChan(channel)

            t_buff = gr1.GetX()
            v_buff = gr1.GetY()
            n = gr1.GetN()
            #t_buff.SetSize(n)
            #v_buff.SetSize(n)
            t = np.frombuffer(t_buff,dtype=float,count=-1)
            v = np.frombuffer(v_buff,dtype=float,count=-1)


            block_number = rawEvent.blockVec[0].getBlock()
            #print('block_number from rawEvent', block_number, np.shape(block_number))
            #Remove first block which is corrupted
            t,v,block_number = RemoveFirstBlock(t,v,block_number)
            #print('block_number from RemoveFirstBlock', block_number, np.shape(block_number))

            t,v,block_number = Calibrator(station,t,v,block_number,str(channel),length,ped_values,kPed,kTime,kVolt)
            #print(t[0])
            #print('block_number from Calibrator', block_number, np.shape(block_number))

            t,v= RemoveBackwardsSamples(t,v)
            f = Akima1DInterpolator(t,v)
            v = f(t_univ)
            #print(v)
            all_volts[i,:]=v
            #all_t[i,:]=t
            all_blocks[i]=block_number

            gr1.Delete()
            usefulEvent.Delete()
            print(v)
            print(t_univ)
            #plt.plot(t_univ,v)
            #plt.show()

        #all_t = all_t[~np.all(all_volts==0,axis=1)]
        all_volts = all_volts[~np.all(all_volts == 0, axis=1)]
        all_blocks = all_blocks[all_blocks !=701]
        t_univ = t_univ-21.0

        np.save("SavedCalibData/time_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy",t_univ)
        np.save("SavedCalibData/volts_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_volts)
        np.save("SavedCalibData/blocks_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_blocks)


        return(t_univ,all_volts)

def LoadNewPed(ARA_run):
    print('load new ped')
    ped_dir = 'data/pedFiles/'
    ped_names = os.listdir(ped_dir)
    run_nums = []
    for i in ped_names:
        run_nums.append(int(i[-8:-4]))
        #print(i[-8:-4])
    best_ped = min(run_nums,key=lambda x:abs(x-int(ARA_run)))
    print('best ped is:', best_ped)
    for i in ped_names:
        if(int(i[-8:-4])==best_ped):
            best_ped_vals = np.loadtxt(ped_dir+i)
            print(np.shape(best_ped_vals))
            return(best_ped_vals,ped_dir)
    return(0)


def LoadSPICElist():

    for i in range(5124,totalEvents):#5124
        #print(i)
        test2 = eventTree.GetEntry(i)


        usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kNoCalib)
        gr1 = usefulEvent.getGraphFromElecChan(1)
        t_buff = gr1.GetX()
        v_buff = gr1.GetY()
        n = gr1.GetN()
        t_buff.SetSize(n)
        v_buff.SetSize(n)
        v1 = np.array(v_buff,copy=True)
        t1 = np.array(t_buff,copy=True)

        gr1 = usefulEvent.getGraphFromElecChan(9)
        t_buff = gr1.GetX()
        v_buff = gr1.GetY()
        n = gr1.GetN()
        t_buff.SetSize(n)
        v_buff.SetSize(n)
        v2 = np.array(v_buff,copy=True)
        t2 = np.array(t_buff,copy=True)

        if(np.max(v1)-np.min(v1)>500 and np.max(v2)-np.min(v2)>500):
            event_list.append(i)
            print(i)
    np.save('SPICE_eventlist.npy',np.asarray(event_list))

def LoadSpiceData(channel,kPed,kTime,kVolt):
    station = '5'
    run = '3989'
    year='2018'
    date='1225'
    length = 2048

    t_file = "SavedCalibData/time_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy"
    v_file = "SavedCalibData/volts_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy"
    #b_file = "/home/kahughes/PA_Analysis/data/SavedCalibData/blocks_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy"

    if(os.path.isfile(t_file) and os.path.isfile(v_file)):
        t=np.load(t_file)
        v=np.load(v_file)
        #b = np.load(b_file)
        return(t,v)

    else:

        #ped_values = DownloadPedFile(year,run,'')
        ped_values, ped_dir = LoadNewPed(run)
        ROOT.gSystem.Load("$ARA_UTIL_INSTALL_DIR/lib/libAraEvent.so")

        
        test = ROOT.TFile.Open("/project2/avieregg/ARA_cal_data/data/root/run00"+str(run)+"/event00"+str(run)+".root")
        calibrator = ROOT.AraEventCalibrator.Instance()
        eventTree = test.Get("eventTree")
        """
        except ReferenceError:
            test = ROOT.TWebFile.Open("http://icecube:skua@convey.icecube.wisc.edu/data/exp/ARA/"+year+"/filtered/L0/ARA05/"+date+"/run"+run+"/event"+run+".root")
            calibrator = ROOT.AraEventCalibrator.Instance()
            eventTree = test.Get("eventTree")
        """
        rawEvent = ROOT.RawAtriStationEvent()

        #unixTime = HK_data.unixTime
        #print('unix time is', unixTime)
        eventTree.SetBranchAddress("event",ROOT.AddressOf(rawEvent))
        totalEvents = eventTree.GetEntries()
        print('total events:', totalEvents)
        #length = 1792

        all_volts = np.zeros([totalEvents,length])
        all_t=np.zeros([totalEvents,length])
        all_blocks=np.zeros([totalEvents])+701

        event_list = []

        event_list = np.load('SPICE_eventlist.npy')
        event_times = []
        print(channel)

        for i in event_list:#totalEvents):
            if(i%10==0):
                print(i/len(event_list))
            test2 = eventTree.GetEntry(i)
            usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kNoCalib)
            gr1 = usefulEvent.getGraphFromElecChan(channel)
            HK_data = rawEvent.unixTime
            event_times.append(HK_data)
            #t_buff = gr1.GetX()
            #v_buff = gr1.GetY()
            #n = gr1.GetN()
            #t_buff.SetSize(n)
            #v_buff.SetSize(n)
            #v = np.array(v_buff,copy=True)
            #t = np.array(t_buff,copy=True)


            t_buff = gr1.GetX()
            v_buff = gr1.GetY()
            n = gr1.GetN()
            t_buff.SetSize(n)
            v_buff.SetSize(n)
            t = np.frombuffer(t_buff,dtype=float,count=-1)
            v = np.frombuffer(v_buff,dtype=float,count=-1)


            block_number = rawEvent.blockVec[0].getBlock()

            #Remove first block which is corrupted
            t,v,block_number = RemoveFirstBlock(t,v,block_number)
            #print(kTime,kPed,kVolt)
            t,v,block_number = Calibrator(station,t,v,block_number,str(channel),length,ped_values,kPed,kTime,kVolt)
            #print(len(v))
            if(len(v)<length):
                #print('here')
                v=np.append(v,np.zeros(length-len(v)))
            if(len(t)<length):
                t=np.append(t,np.zeros(length-len(t)))
            all_volts[i,:]=v

            all_t[i,:]=t
            all_blocks[i]=block_number

            gr1.Delete()
            usefulEvent.Delete()
            #print(t)

        all_t = all_t[~np.all(all_volts==0,axis=1)]
        all_volts = all_volts[~np.all(all_volts == 0, axis=1)]
        all_blocks = all_blocks[all_blocks !=701]

        if(kTime==1):
            all_t,all_volts= RemoveBackwardsSamples(all_t,all_volts)

        #print(t,all_volts)
        np.save("SavedCalibData/time_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_t)
        np.save("SavedCalibData/volts_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_volts)
        np.save("SavedCalibData/blocks_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_blocks)
        np.save("SPICE_event_times.npy",event_times)
        if(wantBlocks==1):
            return(all_t,all_volts,all_blocks)
        else:
            return(all_t,all_volts)

def main():
    '''channel = str(sys.argv[1])#'0'
    station = str(sys.argv[2])

    N1 = [0, 3, 8, 11, 16, 19, 24, 27]
    N2 = [1, 2, 9, 10, 17, 18, 25, 26]
    N_special = [9, 16, 24, 25]

    # for ARA01
    if(station=='1'):
        if(int(channel) in N1):
            rootfile = '975'
            rootfiles = ['975', '964', '967', '968']
        if(int(channel) in N2):
            rootfile = '933'
            rootfiles = ['933', '972', '974', '971']

    pedFile = '/home/mhossain/ARA_Calibration/data/pedFiles/pedMean/pedFile_'+str(rootfile)+'.dat'
    
    #TestFunction("5","5337","0529","2019",0,1,1,1,1,0,0)#run,date,year,channel)
    #LoadDataFromWeb("5","5337","0529","2019",0,1,1,1,1,0,0)#run,date,year,channel)
    #LoadSineWaveData(station, rootfile, pedFile, channel, kPed = 1, kTime = 0, kVolt = 0)
    '''
    print('hello! nothing to see here.')


if __name__=="__main__":
   main()
