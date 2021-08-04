from __future__ import print_function
from ROOT import TCanvas, TGraph
from ROOT import gROOT
import ROOT
import matplotlib.pyplot as plt
from scipy.fftpack import rfft, irfft, fftfreq,fft
from scipy import optimize
from scipy.misc import derivative
import numpy as np
import sys
import math
#sys.path.append('../')
#import global_constants
from math import sin
from pynverse import inversefunc
from AutomaticLoadData_A1 import LoadSineWaveData
import matplotlib

matplotlib.rcParams.update({'font.size': 26})
font = {'weight' : 'bold', 'size' : 18}
matplotlib.rc('font', **font)






# rejecting data outside 2\sigma error
def reject_outliers(data, m = 2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]


def invertedFit(params, tval, vval):
    k = params[0] #frequency
    phi = params[1] #phase
    A = params[2] #amplitude
    T = 1/(k) #time period
    
    '''
    print('Time period, Phase, Freq, Amplitude:', T,phi,k,A)
    print('p0, p1, p2:', params)
    print('tval', tval, np.shape(tval))
    print('vval', vval, np.shape(vval))
    '''
    

    sine = (lambda t: A*np.sin(2*np.pi*k*t - phi))

    t_list = np.linspace(tval - T/2, tval + T/2, 100)
    sine_vals = sine(t_list)
    a_max = np.argmax(sine_vals)
    a_min = np.argmin(sine_vals)
    
    '''
    print('tlist', t_list, np.shape(t_list))
    print('sine_vals', sine_vals, np.shape(sine_vals))
    print('a_max', a_max, np.shape(a_max))
    print('a_min', a_min, np.shape(a_min))
    #print('tlist, sine_vals, a_max, a_min:', t_list, sine_vals, a_max, a_min)
    ''' 
   
    t_max = t_list[a_max]
    t_min = t_list[a_min]
    
    print('t_max, t_min, vval, tval', t_max, t_min, vval, tval)

    if(t_min>t_max):
        minval = t_max
        t_max = t_min
        t_min = minval
        print('t_max, t_min, vval, tval', t_max, t_min, vval, tval)

    #while(t_min>tval):
    #    t_min=tval-0.05
    #while(t_max<tval):
    #    #   print('hello!')
    #    t_max=tval+0.025
    #    #   print(t_max)

    try:
        t_close = inversefunc(sine, y_values = vval, domain = [t_min, t_max])
        print('t_close', t_close, np.shape(t_close))
    except ValueError:
        print('t_max, t_min, vval, tval', t_max, t_min, vval, tval)
        print('params', params)

        #plt.figure(1)
        #plt.plot(np.linspace(0,40,150), sine(np.linspace(0,40,150)))
        #plt.show()

    #jitter= t_close-tval

    print('t_close - tval', t_close-tval)
    return(t_close - tval) #used to be t_close-tval which is basically jitter


def SineFunc(t, k, phi, A): #time, freq, offset, amplitude
    return A*np.sin(2*np.pi*k*t - phi)

def SineFit(t, v, freq):
    #Initial guess for the parameters (length N). If None, then the initial values will all be 1.
    params, params_covariance = optimize.curve_fit(SineFunc, t, v, p0 = [freq, np.pi/2.0, 350]) #,bounds=([-np.inf,-np.inf,200],[np.inf,np.inf,np.inf]))#freq,offset,amplitude,voff
    if(params[2]<0):
        params[2] = np.abs(params[2]) # invert the amplitude
        params[1] = params[1] + np.pi # adding pi phase because the amplitude is inverted
    params[1] = params[1]%(np.pi*2)
    while(params[1]<0):
        params[1] = params[1] + np.pi*2.0
    return(params)


def HistPlotter2D(sample, jitter, loop, channel):
    sample_even = []
    jitter_even = []
    sample_odd = []
    jitter_odd = []

    sample_1d = 4
    jitter_1d = []

    for j in range(0, len(sample)):
        if(sample[j] == sample_1d):
            jitter_1d.append(jitter[j])

        if(sample[j]%2 == 0): #if even
            sample_even.append(sample[j])
            jitter_even.append(jitter[j])
        else:
            sample_odd.append(sample[j])
            jitter_odd.append(jitter[j])
    '''
    print('jitter_1d', jitter_1d, np.shape(jitter_1d))
    print('sample_even', sample_even, np.shape(sample_even))
    print('jitter_even', jitter_even, np.shape(jitter_even))
    print('sample_odd', sample_odd, np.shape(sample_odd))
    print('jitter_odd', jitter_odd, np.shape(jitter_odd))
    '''

    fig, ax1 = plt.subplots(figsize = (12, 8))
    plt.hist2d(sample_even, jitter_even, bins = (128, 128), cmap = 'PRGn', range = np.array([(0.0, 128.0), (-1.0, 1.0)]))
    plt.title('Even Samples: Channel '+str(channel)+' Loop '+str(loop))
    plt.ylabel('Jitter (ns)')
    plt.xlabel('Sample Number')
    plt.colorbar()
    plt.savefig('plots/EvenSamples_Chan'+str(channel)+'_'+str(loop)+'.png')

    fig, ax2 = plt.subplots(figsize = (12, 8))
    plt.hist2d(sample_odd, jitter_odd, bins = (128, 128), cmap = 'RdBu', range = np.array([(0.0, 128.0), (-1.0, 1.0)]))
    plt.title('Odd Samples: Channel '+str(channel)+' Loop '+str(loop))
    plt.ylabel('Jitter (ns)')
    plt.xlabel('Sample Number')
    plt.colorbar()
    plt.savefig('plots/OddSamples_Chan'+str(channel)+'_'+str(loop)+'.png')

    fig, ax3 = plt.subplots(figsize = (12, 8))
    plt.hist(jitter_1d, bins = 25, color = 'c')
    plt.grid(color = 'm', alpha = 0.65, linestyle = 'dotted', linewidth = 0.85)
    plt.title('Timing Offsets of Odd Samples: Channel '+str(channel)+' Loop '+str(loop))
    plt.xlabel('Jitter (ns)')
    plt.ylabel('Counts')
    plt.xlim([-1.0, 1.0])
    plt.savefig('plots/Jitter_1d_Chan'+str(channel)+'_'+str(loop)+'.png')
    return()

def AddOffsets2(t, v, freq, odds):
    #Goals:
    #Fit even blocks to sine waves
    #Fit odd blocks to sine waves
    #compare phase differences between two transitions
    print(odds)
    e_freq = []
    o_freq = []
    #loop over all events:
    for j in range(0, len(v[:, 0])):
        blocks_per_event = int(len(v[0, :])/64)
        #print('total blocks:', blocks_per_event)
        for i in range(0, blocks_per_event):
            #print(i*64,i*64+64,15*64)
            params = SineFit(t[1:64:2], v[j, i*64:i*64+64:2], freq)

            if i%2 == 0:#even block
                e_freq.append(params[0])
            else:
                o_freq.append(params[0])
            #print(params)
    histogram([e_freq, o_freq], 'Frequency')


    t_scaled_e1 = t[:64]*np.mean(e_freq)/freq
    t_scaled_o = t[64:128]*np.mean(o_freq)/freq
    t_scaled_e2 = t[128:192]*np.mean(e_freq)/freq

def AddOffsets(t, v, freq, odds, old_mean_o2e, old_mean_e2o):
    print(len(v[0, :]))

    e_freq = []
    o_freq = []
    for j in range(0, len(v[:, 0])):
        val = 0
        for i in range(0, 14):
            params = SineFit(t[odds[:32]], v[j, odds[val:val + 32]], freq)
            if i%2 == 0: #even block
                e_freq.append(params[0])
            else:
                o_freq.append(params[0])
            val = val + 32

    t_scaled_e1 = t[:64]*np.mean(e_freq)/freq
    t_scaled_o = t[64:128]*np.mean(o_freq)/freq
    t_scaled_e2 = t[128:192]*np.mean(e_freq)/freq

    e2o_diff = []
    o2e_diff = []
    #plt.figure(0)
    #plt.plot(even_time[odds[:32]],v[0,odds[0:32]])
    #plt.show()

    for j in range(0, len(v[:, 0])):
        val = 0
        for i in range(0, 13):
            #print(i)
            if(i%2 == 0):
                params = SineFit(t_scaled_e1[odds[:32]], v[j, odds[val:val+32]], freq)
                #print(params[0])
                e_offset = (params[1])
                params = SineFit(t_scaled_o[odds[:32]], v[j, odds[val + 32:val + 64]], freq)
                o_offset = (params[1])
                diff = e_offset - o_offset
                while(diff>np.pi):
                    diff = diff - 2*np.pi
                while(diff<-1*np.pi):
                    diff = diff + 2*np.pi
                #print(e_offset,o_offset)
                e2o_diff.append(((diff))/(2*np.pi*freq))
            else:
                #print(oe_time[:32])
                #print(len(oe_time[:32]))
                params = SineFit(t_scaled_o[odds[:32]], v[j, odds[val:val + 32]], freq)
                o_offset = (params[1])
                params = SineFit(t_scaled_e2[odds[:32]], v[j, odds[val + 32:val + 64]], freq)
                e_offset = (params[1])
                diff = e_offset - o_offset
                while(diff>np.pi):
                    diff = diff - 2*np.pi
                while(diff<-1*np.pi):
                    diff = diff + 2*np.pi
                o2e_diff.append(((diff))/(2*np.pi*freq))
            val = val + 32

    #o2e_diff = reject_outliers(np.asarray(o2e_diff))
    #histogram([e2o_diff,o2e_diff],'Wrap Around Time (ns)')
    #Scale time so that offset is taken into account
    e2o_mean = np.abs(np.mean(e2o_diff) + old_mean_e2o)
    o2e_mean = np.abs(np.mean(o2e_diff) + old_mean_o2e)
    print(e2o_mean, o2e_mean, old_mean_e2o, old_mean_o2e)
    if(e2o_mean>1.0):
        e2o_mean = 0.31
    if(o2e_mean>1.0):
        o2e_mean = 0.31
    print('the offsets are:', e2o_mean, o2e_mean)
    
    #histogram([e2o_diff,o2e_diff],'time offset')
    t_updated = np.zeros(wf_len)
    val = 0
    for i in range(0, 7):
        if i==0:
            t_updated[val:val + 64] = t_scaled_e1
            t_updated[val + 64:val + 128] = t_scaled_o + e2o_mean
        else:
            t_updated[val:val + 64] = t_scaled_e1 + o2e_mean + t_updated[val - 1]
            t_updated[val + 64:val + 128] = t_scaled_o + e2o_mean + t_updated[val - 1]
        #else:
        #    t_updated[val:val+64]=even_time+(o2e_mean+t_updated[val-1])
        #t_updated[val+64:val+128]=odd_time+(e2o_mean+t_updated[val+64-1])
        val = val + 128
    print(t_updated)
    return(t_updated, o2e_mean, e2o_mean)

    #print(t)
def histogram(vals, string):
    fig, axs = plt.subplots(2, 1, facecolor = 'w')
    counter = 0
    for ax in axs.reshape(-1):
        print('success')
        ax.hist(vals[counter], color = 'green', edgecolor = 'none', bins = 20)
        ax.axvline(x = np.mean(vals[counter]), color = 'red', ls = '-', linewidth = 2.0)
        #ax.text(230,250,"mean (MHz):  "+str(round(np.mean(freq_array[counter]*1000),2)))
        #ax.set_xlim(200,250)
        #ax.set_ylim(0,300)
        ax.set_xlabel(string)
        ax.set_ylabel('Counts')
        counter = counter + 1
    plt.savefig('plots/Hist.png')
    plt.show()

def SinePlotter(t, v, params, sample, loop, channel):
    fig, ax = plt.subplots(figsize = (12, 8))
    plt.scatter(t, v[sample, :], label = 'Data Points', color = 'b')
    print('sample from sineplotter', sample, np.shape(sample))
    print('t from sineplotter', t, np.shape(t))
    print('v from sineplotter', v[sample,:], np.shape(v[sample,:]))
    station = 1
    rootfiles = 975
    wf_title = 'Pre-Deployment Waveform: ARA'+str(station)+' run'+str(rootfiles)+' ElecChan('+str(channel)+')'
    t_up = np.linspace(t[0], t[-1], len(t)*50)
    print('t_up from sineplotter', t_up, np.shape(t_up))
    plt.plot(t_up, SineFunc(t_up, params[sample, 0], params[sample, 1], params[sample, 2]), '-', color = 'g', linewidth = 1, label = 'Best Fit: '+str(np.round(params[sample, 0]*1000, 3))+' MHz')
    plt.grid(color = 'm', alpha = 0.65, linestyle = 'dotted', linewidth = 0.85)
    plt.legend()
    plt.title(wf_title)
    plt.xlabel('Time (ns)')
    plt.ylabel('ADC')
    #plt.show()
    plt.savefig('plots/SineWaveFit_Chan'+str(channel)+'_'+str(loop)+'.png')



def SlopeFinder(t, v, sample):
    if sample==0:
        slope = (v[sample + 2] - v[sample])/(t[sample + 2] - t[sample])
    else:
        slope = (v[sample] - v[sample - 2])/(t[sample] - t[sample - 2])

    return slope



def CorrectTimingSample(rootfile, channel, freq, t_cal, station):
    wf_len = 384  # this length should be a multiple of 128 as 128 samples are calibrated in one go
    print('length of the waveform is', wf_len)
    
    pedFile = '/home/mhossain/ARA_Calibration/data/pedFiles/pedMean/pedFile_'+str(rootfile)+'.dat' #pedestal file directory

    print('pedestal file is:', pedFile)

    all_times, volt, blocks = LoadSineWaveData(station, rootfile, pedFile, channel, kPed = 1, kTime = 0, kVolt = 0)

    time = all_times[0] - all_times[0][0] #uncalibrated time
    print('all_times', all_times, np.shape(all_times))
    print('volt', volt, np.shape(volt))
    print('blocks', blocks, np.shape(blocks))
    print('time = all_times[0] - all_times[0][0]', time, np.shape(time))
    print('number of events is', np.shape(volt))

    num_blocks = len(volt[:, 0])
    print('num_blocks', num_blocks, np.shape(num_blocks))
    best_params = np.zeros([num_blocks, 4])

    odds = np.linspace(1, wf_len - 1, int(wf_len/2), dtype = int)
    evens = np.linspace(0, wf_len - 2, int(wf_len/2), dtype = int)

    odd_params = np.zeros([num_blocks, 3])
    even_params = np.zeros([num_blocks, 3])

    print('odds', odds, np.shape(odds))
    print('evens', evens, np.shape(evens))

    t_cal = np.zeros(128)


    if(t_cal[5]==0.0):
        print('clearing out old t_cal')
        t_cal = time[:128]
        print(t_cal)

    #set calibrated time array to be equal to uncalibrated time array, so length is the same

    t_cal_full = time 
    #print('t_cal before is', t_cal_full)
    #print(np.shape(t_cal_full),np.shape(volt))
    odd_mean = 0.0
    even_mean = 0.0
    for l in range(0, 4):
        print('loop number ', l)
        #print(np.shape(t_cal_full),np.shape(volt))
        #First fix offsets between blocks, as that can be larger than the channel to channel fixes.
        #if l==0:
        #    t_cal_full,odd_mean,even_mean=AddOffsets(t_cal_full,volt,freq,odds,odd_mean,even_mean)
            #t_cal_full,odd_mean,even_mean=AddOffsets2(t_cal_full,volt,freq,odds)
        #Fit each waveform to a sine wave
        volt = volt[:, :len(t_cal_full)]
        print('t_cal_full', t_cal_full, np.shape(t_cal_full))


        #STEP ONE: For each event, calculate the best sine wave fit. Then remove outlier frequencies and calculate the mean frequency.
        #Since we know the expected frequency from the lab data, we can then correct the overall time step.

        for i in range(0, num_blocks):
            odd_params[i, :] = SineFit(t_cal_full[odds], volt[i, odds], freq)
            #SinePlotter(t_cal_full,volt,odd_params,i)

        #plt.scatter(t_cal_full[odds],volt[i,odds])
        #t_up = np.linspace(t_cal_full[0],t_cal_full[-1],500)
        #plt.plot(t_up,SineFunc(t_up,odd_params[i,0],odd_params[i,1],odd_params[i,2]))
        #plt.show()

        
        #if(l>-1):
        fig, ax = plt.subplots(figsize = (12, 8))
        plt.hist(odd_params[:, 0]*1000, bins = 20, color = 'green', alpha = 0.65)
        plt.grid(color = 'm', alpha = 0.65, linestyle = 'dotted', linewidth = 0.85)
        plt.axvline(np.mean(odd_params[:, 0])*1000, color = 'r', linestyle = 'dashed', label = 'Mean Frequency = '+str(np.round(np.mean(odd_params[:, 0])*1000, 3))+' MHz')
        plt.axvline(214.0, color = 'b', linestyle = 'dashed', label = 'Input Frequency = '+str(freq*1000)+' MHz')
        plt.xlabel('Best Fit Frequency (MHz)')
        plt.legend()
        #plt.show()
        plt.savefig('plots/Frequency_hist_Chan'+str(channel)+'_'+str(l)+'.png')
        
        SinePlotter(t_cal_full, volt, odd_params, 4, l, channel)
        
        freq_no_outliers = reject_outliers(np.asarray(odd_params[:, 0]))
        print('freq_no_outliers', freq_no_outliers, np.shape(freq_no_outliers))
        mean_freq = np.mean(freq_no_outliers)
        print('mean frequency is', mean_freq)
        #histogram([odd_params[:,1],even_params[:,1]],'')

        #Scale timing to reflect true frequency
        t_cal_full=t_cal_full*mean_freq/freq
        print('normalized t_cal_full', t_cal_full, np.shape(t_cal_full))
        print('Fitting to sine:')
        #Re-fit using new time
        for i in range(0, num_blocks):

            odd_params[i,:] = SineFit(t_cal_full[odds], volt[i, odds], freq)
            #SinePlotter(t_cal_full,volt,odd_params,i)
            #print('here')
            #print(l)
            #if(l>0):
                #print("here!")
                #SinePlotter(t_cal_full,volt,odd_params,i)
                #plt.show()

        t_cal = t_cal_full[:128]

        jitter_array = []
        sample_array = []
        slope_array = []
        jitter_slope = []
        new_spacing = np.zeros(128) #spacing between 0 and 1, 1 and 2, etc.
        #Go through each sample and correct based on fit
        '''
        if(int(channel)==24):
            cutval = 5.0
        else:
            cutval = 30.0'''
        cutval = 30.0
        print('here is the slow part')
        for k in range(0, wf_len):
            counter = 0
            if(k%50==0):
                print(k/wf_len*100,'percent done.')
            for i in range(0, num_blocks):
                if(np.abs(volt[i, k])<cutval and (freq - odd_params[i, 0])<0.002):# and np.abs(odd_params[i,2]>200)):
                #if(np.abs(volt[i, k])<cutval and (freq - odd_params[i, 0])<30.0):
                    try:
                        print('odd_params[i, :]', odd_params[i, :])
                        print('tval', t_cal_full[k])
                        print('vval', volt[i, k])
                        invert_fit = invertedFit(odd_params[i, :], t_cal_full[k], volt[i, k])
                        print("Inverted Fit", invert_fit)
                        jitter_array.append(invert_fit)
                        sample_array.append(k%128)
                        #print('jitter_array', jitter_array)
                        #print("sample_array", sample_array)
                        print('loop', l)
                        counter = counter + 1 
                    except:
                        #print('odd_params[i, :]', odd_params[i, :], np.shape(odd_params[i, :]))
                        #print('t_cal_full[k]', t_cal_full[k], np.shape(t_cal_full[k]))
                        #print('volt[i, k]', volt[i, k], np.shape(volt[i, k]))
                        print("error in finding inverse!!")
                        print('loop', l)

            t_cal_full[k] = t_cal_full[k] + np.mean(jitter_array[-counter:])
            print("t_cal_full with jitter", t_cal_full[k], np.shape(t_cal_full[k]))
            if(k>0):
                new_spacing[k%128] = new_spacing[k%128] + t_cal_full[k] - t_cal_full[k - 1]


        new_spacing[1:] = new_spacing[1:]/3.0 # This is divided by 3 as wf_len is 128*3
        new_spacing[0] = new_spacing[0]/2.0 # As wf_len is 128*3
        #print('spacing is', new_spacing)



        for i in range(0, wf_len):
            if(i==0):
                t_cal_full[i] = 0.0
            else:
                t_cal_full[i] = t_cal_full[i - 1] + new_spacing[(i)%128]


        #print('final t_cal is',t_cal_full)

        t_cal = t_cal_full[:128]
        print('final t_cal is', t_cal, np.shape(t_cal))

        #if(l>0):
        #    SinePlotter(time[odds],volt[:,odds],odd_params,5)


        """
        plt.figure(6,facecolor='w')
        plt.hist2d(slope_array,jitter_slope,bins=(250,128),cmap=plt.cm.jet,range=np.array([(-1100.0,1100.0),(-1.0,1.0)]))
        plt.title('Even Samples vs Slope')
        plt.show()
        """


        if(l<1):
            np.save('/home/mhossain/ARA_Calibration/data/ARA'+str(station)+'_cal_files/samples_'+rootfile+'_'+channel+'first.npy', np.asarray(sample_array))
            np.save('/home/mhossain/ARA_Calibration/data/ARA'+str(station)+'_cal_files/jitter_'+rootfile+'_'+channel+'first.npy', np.asarray(jitter_array))
        #if(l>-1):
        HistPlotter2D(sample_array, jitter_array, l, channel)
    #print('final t_cal is', t_cal_full)
    np.save('/home/mhossain/ARA_Calibration/data/ARA'+str(station)+'_cal_files/t_cal_'+channel+'.npy', t_cal_full)
    np.save('/home/mhossain/ARA_Calibration/data/ARA'+str(station)+'_cal_files/samples_'+channel+'final.npy', np.asarray(sample_array))
    np.save('/home/mhossain/ARA_Calibration/data/ARA'+str(station)+'_cal_files/jitter_'+channel+'final.npy', np.asarray(jitter_array))
    #HistPlotter2D(sample_array,jitter_array)
    print('t_cal after all loops is', t_cal)
    return(t_cal)

def main():
    channel = str(sys.argv[1])#'0'
    station = str(sys.argv[2])


    N1 = [0, 1, 2, 3, 8, 9, 10, 11]
    N2 = [16, 17, 18, 19, 24, 25, 26, 27]
    N_special = [9, 16, 24, 25]
    
    # for ARA01
    if(station=='1'):
        if(int(channel) in N2):
            rootfile = '975' # Pedestal File Number
            rootfiles = ['975', '964', '967', '968'] # CW Wave Number
        if(int(channel) in N1):
            rootfile = '933'
            rootfiles = ['933', '972', '974', '971']


    # for ARA04 and ARA05
    if(station=='5'):
        if(int(channel) in N1):
            rootfile = '1402'
            rootfiles = ['1402', '1403', '1404', '1405']
        if(int(channel) in N2):
            rootfile = '1411'
            rootfiles = ['1411', '1412', '1413', '1414']
    if(station=='4'):
        if(int(channel) in N1 and int(channel) not in N_special):
            rootfile = '2829'
            rootfiles = ['2829', '2830', '2831', '2832']
        if(int(channel) in N2 and int(channel) not in N_special):
            rootfile = '2840'
            rootfiles = ['2840', '2841', '2842', '2843']
        if(int(channel)in N_special):
            rootfiles = ['2855', '2856']


    #sample_final=np.load('ARA'+str(station)+'_cal_files/samples_'+rootfile+'_'+channel+'final.npy')
    #jitter_final = np.load('ARA'+str(station)+'_cal_files/jitter_'+rootfile+'_'+channel+'final.npy')

    #sample_first=np.load('ARA'+str(station)+'_cal_files/samples_'+rootfile+'_'+channel+'first.npy')
    #jitter_first = np.load('ARA'+str(station)+'_cal_files/jitter_'+rootfile+'_'+channel+'first.npy')
    #HistPlotter2D(sample_first,jitter_first)
    #HistPlotter2D(sample_final,jitter_final)


    freqs = [0.214, 0.218, 0.521, 0.702] # frequency used for ARA01 is 214 MHz

    cal_t = np.zeros(128)

    #jitter = np.load('jitter_1404_3final.npy')
    #samples = np.load('samples_1404_3final.npy')
    #HistPlotter2D(samples,jitter)
    #CorrectTimingSample(rootfiles[0],channel,freqs[0],cal_t,station)

    for a in range(0, 1):
        #exists = os.path.isfile('ARA'+str(station)+'_cal_files/jitter_'+rootfiles[a]+'_'+channel+'final.npy')
        #if(exists):
        #    print('file exists!')
        #else:
        #try:
        CorrectTimingSample(rootfiles[a], channel, freqs[a], cal_t, station)
        #except:
        #    print('Error')

    #average all results together
    #average_tcals('ARA'+str(station)+'_cal_files/',channel,rootfiles)

    #account for wrap around time
    #FindWrapAround(rootfiles[1],channel,freqs[1])

if __name__=="__main__":
   main()

