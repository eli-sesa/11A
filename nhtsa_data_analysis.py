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

timeStep = 1/20000

cfcClass = 20
SAVE_CSV = False

maxStroke = 2
maxTime = .150
Xmax = 2
Ymax = .600

filenames = ['AX.csv', 'AY.csv']

fig, axs = plt.subplots(4, 2, sharey='row', sharex=True)
fig.suptitle("CFC %d"%cfcClass)

dataOut = pd.DataFrame()
output = pd.DataFrame(
    columns=(
        'TestNo',
        'Direction',
        'MaxAccX', 
        'MaxAccY', 
        'XFinal',
        'YFinal',
        'TwhereX',
        'TwhereY',
        'MaxXError',
        'FinalXError',
        'MaxYError',
        'FinalYError',
        ))

plotcol = 0
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
            
        posError = 1000 * (posraw - pos) # position error due to filtering in mm

        # Plot everything 
        axs[0, plotcol].set_title(filename.replace('.csv',''))
        
        axs[0, plotcol].plot(time, accFilt)
        axs[1, plotcol].plot(time, vel)
        axs[2, plotcol].plot(time, pos)
        axs[3, plotcol].plot(time, posError)
        
        
        # Compute things for  Data table
    
        # imax = np.argmax(posError) # index of the max position error
        # print('Final Position Error = %f mm' %posError[-1])
        # print("Max Position error of %f mm at %f seconds" %(posError[imax], time[imax]))
        # print('\n\n')
        
        output['TestNo'] = headers
        output['Direction'] = filename[1]
        # This doesnt work but leaving here for readability while I figure it out 
        # if filename[1] == 'X':
        #     output['MaxAccX'].append(np.max(accFilt))
        #     output['XFinal'] = time[(time - maxTime).argmin()]
        #     output['TwhereX'] = time[(pos - Xmax).argmin()]
        #     output['MaxXError'] = np.max(posError)
        #     output['FinalXError'] = posError[-1]
               
        outputStats = np.array([headers,
                                filename[1],
                                np.max(accFilt),
                                pos[np.abs(time - maxTime).argmin()], 
                                time[np.abs(pos - Xmax).argmin()],
                                np.max(posError),
                                posError[-1]])
        
        print(outputStats)
        
        
        if filename[1] == 'Y':
            output['MaxAccY'] = np.max(accFilt)
            output['YFinal'] = time[np.abs(time - maxTime).argmin()]
            output['TwhereY'] = time[np.abs(pos - Xmax).argmin()]
            output['MaxYError'] = np.max(posError)
            output['FinalYError'] = posError[-1] 
        
        
        # Output filtered data to csv
        dataOut[headers] = accFilt       
       # scipy.signal.findpeaks
    # At the end, increase the column to plot into
    plotcol +=1

    if SAVE_CSV:
        dataOut.to_csv(filename[0:2]+'Filtered'+str(cfcClass)+'.csv')

# With the figure, add ylabels to the left column        
ylabels = ('Acc (g)', 'Vel (m/s)', 'Pos (m)', ' $\Delta$ Pos (mm)')
for i in range(4):
    axs[i, 0].set_ylabel(ylabels[i])

plt.legend(data.columns)
axs[2, 0].hlines(2.0, time[0], .150, linestyle='dashed', zorder=1)
axs[2, 0].vlines(.150, 0, 2.0, linestyle='dashed',  zorder=1)

axs[2, 1].hlines(.6, time[0], .150,  linestyle='dashed', zorder=1)  
axs[2, 1].vlines(.150, 0, .6, linestyle='dashed', zorder=1)      

#%% 
# Yaw Data

# wRaw = data['9586']
# wFilt = cfc(wRaw, cfcClass, timeStep)
dataOut = pd.DataFrame()


filename = 'WZ.csv'

data = pd.read_csv(filename, usecols=range(1, 11))
fig, axs = plt.subplots(4, 1, sharex=True)
fig.suptitle('$\omega_z$')
for headers in data.columns:
      wRaw = np.array(data[headers])
      wFilt = cfc(wRaw, cfcClass, timeStep)
      theta = scipy.integrate.cumtrapz(wFilt, time, initial=0)
      
      thetaRaw = scipy.integrate.cumtrapz(wRaw, time, initial=0)
      angleError = (thetaRaw - theta)
      
      axs[0].plot(time, filt_dif(wRaw, cfcClass, timeStep))
      axs[0].set_ylabel('$ \\alpha$ (deg/s^2)')

      axs[1].plot(time, wFilt)
      axs[1].set_ylabel('$\omega_z$ (deg/s)')

      axs[2].plot(time, theta)
      axs[2].set_ylabel('$\Theta$ (deg)')

      axs[3].plot(time, angleError)
      axs[3].set_ylabel('Angle Error (deg)')
      
      
      dataOut[headers] = filt_dif(wRaw, cfcClass, timeStep)

plt.legend(data.columns)

if SAVE_CSV:
    dataOut.to_csv(filename[0:2]+'Filtered'+str(cfcClass)+'.csv')


#%%
plt.figure()
plt.plot(filt_dif(wRaw, cfcClass, timeStep),'x')
plt.plot(dif_filt(wRaw, cfcClass, timeStep),'o')
plt.plot(spectral_dif(cfc(wRaw, cfcClass, timeStep), timeStep))

plt.legend(('filter then diff','diff then filter','spectral derivitave'))
