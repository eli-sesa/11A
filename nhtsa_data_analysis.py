# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 13:58:12 2020

@author: eli
"""

import pandas as pd
import numpy as np
import scipy.signal
import scipy.integrate
import matplotlib.pyplot as plt


def cfc(data, cfc_class, T):
    #Apply the 2 pole butterworth backward and forward to make a cfc filter
    # Math from http://zone.ni.com/reference/en-XX/help/370858P-01/crash/misc_cfc/
    
    def J211(data, cfc_class, T):
        X = data
        Y = np.zeros_like(X)
        
        wd = 2*np.pi * cfc_class * 2.0775
        wa = (np.sin(wd * T/2)) / (np.cos(wd * T/2))
        a0 = wa**2 / (1 + np.sqrt(2)*wa + wa**2)
        a1 = 2*a0
        a2 = a0
        b1 = -2 * (wa**2 - 1) / (1 + np.sqrt(2)*wa + wa**2)
        b2 = (-1 + np.sqrt(2)*wa - wa**2) / (1 + np.sqrt(2)*wa + wa**2)

        for i in range(3, len(X)):
            Y[i] = a0*X[i] + a1*X[i-1] + a2*X[i-2] + b1*Y[i-1] + b2*Y[i-2]
            
        return Y

    Y2 = np.flip(
        J211(
            np.flip(
                J211(
                    data, cfc_class, T)),
            cfc_class, T)
        )
    return Y2



def finite_difference(x, timeStep):  
    dx = np.zeros(len(x))
    for i in range(len(x)-1):
        dx[i] = (x[i+1] - x[i]) / timeStep
    dx[-1] = dx[-2]
    return dx




timeStep = 1/20000

cfcClass = 20

filenames = ['AX.csv', 'AY.csv']
fig, axs = plt.subplots(4, 2, sharey='row', sharex=True)

fig.suptitle("CFC %d"%cfcClass)

plotcol = 0

dataOut = pd.DataFrame()
# Linear Accelerations
for filename in filenames:
    
    data = pd.read_csv(filename, usecols=range(1,11))
    time = np.arange(-.05, .3+timeStep, timeStep)
    
    for headers in data.columns:  
        
        # Filter and Integrate
        accFilt = cfc(data[headers], cfcClass, timeStep)
        vel = scipy.integrate.cumtrapz(9.81 * accFilt, time, initial=0)
        pos = scipy.integrate.cumtrapz(vel, time, initial=0)
        
        # compute raw position without all of the extra steps.  Messy, but saves variables
        posraw = scipy.integrate.cumtrapz(
            scipy.integrate.cumtrapz(
                9.81 * data[headers],
                time, 
                initial=0)
            ,time,
            initial=0)
            
        poserror = 1000 * (posraw - pos) # position error due to filtering in mm

        # Plot everything 
        axs[0, plotcol].set_title(filename.replace('.csv',''))
        
        axs[0, plotcol].plot(time, accFilt)
        axs[1, plotcol].plot(time, vel)
        axs[2, plotcol].plot(time, pos)
        axs[3, plotcol].plot(time, poserror)
        
        
        # Compute things for  Data table
        
        # imax = np.argmax(poserror) # index of the max position error
        # print('Final Position Error = %f mm' %poserror[-1])
        # print("Max Position error of %f mm at %f seconds" %(poserror[imax], time[imax]))
        # print('\n\n')
        
        # scipy.signal.findpeaks
        
        # Output filtered data
        dataOut[headers] = accFilt
        
        
    # At the end, increase the column to plot into
    plotcol +=1


    dataOut.to_csv(filename[0:2]+'Filtered'+str(cfcClass)+'.csv')

# With the figure, add ylabels to the left column        
ylabels = ('Acc (g)', 'Vel (m/s)', 'Pos (m)', ' $\Delta$ Pos (mm)')
for i in range(4):
    axs[i, 0].set_ylabel(ylabels[i])

plt.legend(data.columns)
        
#%%
# Yaw Data


    
    
    
    
    
    
    
    
    
    
#     ############################# Old WOrk 
#     wRaw = np.array(data[headers])
#     np.fft.fft(wRaw)
    
#     wFilt = cfc(wRaw, cfcClass, timeStep)
#     alphaFilt = finite_difference(wFilt, timeStep)

#     # plt.plot(time, alphaRaw)
#     plt.plot(time, alphaFilt)
    
#     # wFilt = cfc(data[headers], cfcClass, timeStep)
    
    
#     # alphaRaw = scipy.integrate.cumtrapz(wRaw, time, initial=0)
#     # alphaFilt = scipy.integrate.cumtrapz(wFilt, time, initial=0)
    
     
    
#     # axs[0].plot(time, cfc(alphaRaw, cfcClass, timeStep))
#     # axs[0].set_ylabel('Integrate Raw, then Filter')
#     # axs[1].plot(time, alphaFilt)
#     # axs[1].set_ylabel('Integrate Filtered')
#     # axs[2].plot(time, wRaw)

# plt.legend(data.columns)
#%% 
# Derivitive Sandbox
# https://www.youtube.com/watch?v=y8SqkjoKV4k

# wRaw = data['9586']
# wFilt = cfc(wRaw, cfcClass, timeStep)
dataOut = pd.DataFrame()


def filt_dif(x, filterClass, timeStep):
    # Filter the original signal using cfc, then take finite difference
    xFilt = cfc(x, filterClass, timeStep)
    dx = finite_difference(xFilt, timeStep)
    return dx

def dif_filt(x, filterClass, timeStep):
    # Take finite difference and then filter the result
    dxRaw = finite_difference(x, timeStep)
    dx = cfc(dxRaw, filterClass, timeStep)
    return dx

def spectral_dif(x, timeStep):
    # Take derivitive using fft per https://www.youtube.com/watch?v=y8SqkjoKV4k
    # Note: Does not seem to work with un-filtered signal
    xHat = np.fft.fft(x)
    omega = np.fft.fftfreq (x.size, timeStep/(2*np.pi))
    dxHat = omega * xHat * [1j] # Uses the fact that we know  fft(dx/dt) = i*omega*fft(x)
    dx = np.real(np.fft.ifft(dxHat))
    return dx

filename = 'WZ.csv'

data = pd.read_csv(filename, usecols=range(1, 11))
fig, axs = plt.subplots(4, 1, sharex=True)

for headers in data.columns:
      wRaw = np.array(data[headers])
      wFilt = cfc(wRaw, cfcClass, timeStep)
      theta = scipy.integrate.cumtrapz(wFilt, time, initial=0)
      
      axs[0].plot(time, filt_dif(wRaw, cfcClass, timeStep))
      axs[0].set_ylabel('Acc: Filter then Diff')
      
      axs[1].plot(time, dif_filt(wRaw, cfcClass, timeStep))
      axs[1].set_ylabel('Acc: Diff then Filter')

      axs[2].plot(time, wFilt)
      axs[2].set_ylabel('$\omega_z$ (deg/s)')

      axs[3].plot(time, theta)
      axs[3].set_ylabel('$\Theta$ (deg)')

      dataOut[headers] = filt_dif(wRaw, cfcClass, timeStep)

plt.legend(data.columns)
    
dataOut.to_csv(filename[0:2]+'Filtered'+str(cfcClass)+'.csv')

# plt.plot(time, filt_dif(wRaw, cfcClass, timeStep),'.')
# plt.plot(time, dif_filt(wRaw, cfcClass, timeStep),'.')
# plt.plot(time[50:-50], spectral_dif(wFilt[50:-50], timeStep),'.')





#%%
######WORKSPACE#######
L = wRaw.size * timeStep
# Compute FD of filtered data
alphaRaw = finite_difference(wRaw, timeStep)
alphaFilt = finite_difference(wFilt, timeStep)

# Compute Derivitive of data w FFT.  Let "hat" denote fft'd signal
wRawHat = np.fft.fft(wFilt)
omega = np.fft.fftfreq (wFilt.size, timeStep/(2*np.pi))
alphaHat = omega * wRawHat * [1j] # Uses the fact that we know  fft(dx/dt) = i*omega*fft(x)
alphaFFT = np.real(np.fft.ifft(alphaHat))

plt.figure()

plt.plot(time, alphaFilt,'.')
plt.plot(time, cfc(alphaRaw, cfcClass, timeStep),'.')
plt.plot(time, alphaFFT,'.')


plt.legend(('Filter then FD','FD then Filter', 'Spectral Derivitive'))


# Wrap all in functions if needed





