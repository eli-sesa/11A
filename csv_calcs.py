# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 16:50:05 2020

@author: eli
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


def rows_to_skip(path):
    rowsToSkip = 24
    with open(path) as f: 
       for ind, line in enumerate(f):
           if line.startswith('AccTgt'):
               rowsToSkip = ind
       return rowsToSkip

def calcs(data, scaleLoadCell=False):
    
    ### Constants
    mLoadCell = 156.005024; # lb, cell AB18
    dPiston = 5.125; # in
    dRod = 2.125; # in
    areaRatio = dPiston**2 / (dPiston**2 - dRod**2);
    K = 1.09; # left right correction factor
     
    # Initilaize constant pre-clamp value.  This could bemore sophisticated if desired
    iPreClamp = 301 #index of when we pull stable preclamp pressure
    p0PreClamp = data['P_PreclampYaw'].iloc[iPreClamp]
    
    # If we find the pulse in the list of ones we know to scale, do it
    if scaleLoadCell:
        data['YawLoadcellLeft'] /= 5
        data['YawLoadcellRight'] /= 5

    # Big Block 'O Computation.  Easy to add new channels as desired
    dataFull = data.copy()
    dataFull["fTarget"] = 51000/5000 * data['PnetAcqFire_Yaw']
    dataFull['fTotalAcq'] = data['YawLoadcellLeft'] + data['YawLoadcellRight']
    dataFull['fLeft'] = data['YawLoadcellLeft'] - mLoadCell * data['Acc'] 
    dataFull['fRight'] = data['YawLoadcellRight'] - mLoadCell * data['Acc']
    dataFull['fTotal'] = dataFull['fLeft'] + dataFull['fRight']
    dataFull['pLeftEstimate'] = p0PreClamp * areaRatio + data['PnetAcqFire_Yaw']/2 * (1-(K-1))
    dataFull['pRightEstimate'] = p0PreClamp * areaRatio - data['PnetAcqFire_Yaw']/2 * K
    dataFull['cofLeft'] =  dataFull['fLeft'] / (2*dataFull['pLeftEstimate'] * np.pi/4 * (dPiston**2 - dRod**2))
    dataFull['cofRight'] = dataFull['fRight'] / (2*dataFull['pRightEstimate'] * np.pi/4 * (dPiston**2 - dRod**2))
    dataFull['deltaF'] = dataFull['fRight'] - dataFull['fLeft'] 
    dataFull['pvleft'] =  dataFull['pLeftEstimate'] * data['VelAcq']
    dataFull['pvRight'] = dataFull['pRightEstimate'] * data['VelAcq']
    
    return dataFull


# Manually enter the files that we need to scale the loadCell datas
pulsesToScale = ["YB Brake In X 5kph_testy_2020-07-16-1428_2020-07-16-1438_eng.csv", 
                 "YB Brake In X 5kph_testy_2020-07-16-1504_2020-07-16-1510_eng.csv", 
                 "YB Brake In X 5kph_testy_2020-07-16-1526_2020-07-16-1534_eng.csv", 
                 "YB Brake In X 5kph_testy_2020-07-17-1047_2020-07-17-1055_eng.csv", 
                 "YB Brake In X_2020-07-17-1136_2020-07-17-1142_eng.csv", 
                 "YB Brake In X_2020-07-17-1201_2020-07-17-1210_eng.csv", 
                 "YB Brake In X_2020-07-17-1227_2020-07-17-1233_eng.csv", 
                 "YB Brake In X_2020-07-17-1248_2020-07-17-1253_eng.csv", 
                 "YB Brake In X_2020-07-17-1510_2020-07-17-1516_eng.csv", 
                 "10 kph Step X_2020-07-17-1603_2020-07-17-1609_eng.csv", 
                 "20 kph Step X_2020-07-20-0934_2020-07-20-0943_eng.csv", 
                 "20 kph Step X_2020-07-20-0956_2020-07-20-1005_eng.csv"]

# Specify what prefix to add to our new files
prefix = 'NEW'

for filename in os.listdir(os.getcwd()):
    #avoid re-doing for already created files
    if filename.endswith('csv') and not filename.startswith(prefix): 
        
        print(filename)
        rowsToSkip = rows_to_skip(filename)
        # Is in the load cell scaling?
        editLoadCell = filename in pulsesToScale
        # Load in the data
        data = pd.read_csv(filename, skiprows=rowsToSkip)
        # Add our new columns
        dataFull = calcs(data, scaleLoadCell=editLoadCell)
        # And save the csv file with a prefix
        dataFull.to_csv(prefix+filename, index=False)






















