# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%

# Specify input and output file locations
path = os.getcwd()
dataPath = path + '\\' + 'destination'


### Constants ####

timeStep = 1/10000
tFinal = 0.2
rowsToSkip = 23

samplePath = os.listdir(dataPath)[1]
#loop through all of the eng files that I have
data_store = []
cof_store = []


csv_list = []
for file in os.listdir('destination'):

    if file.endswith('csv'):
        csv_list.append(file)
         


for filename in csv_list:
    absPath = 'C:/E/0Design Work/Dynamic Yaw/COF Data Investigation/destination/'
    destination = absPath+filename

    with open(destination) as f:
        for ind, line in enumerate(f):
            if line.startswith('AccTgt'):
                rowsToSkip = ind
                
        
        # Read in csv to dataframe
    dataFull = pd.read_csv(absPath+filename, skiprows=rowsToSkip)
    dataFull['Time'] = np.arange(0, dataFull.shape[0] * timeStep, step=timeStep)
        # add test identifier column.
        #append to giant dataframe
    # print(filename)
    # print(type(dataFull.iloc[0,0]))
    # usefulData  = ['Time','MasterValveCmdCorrected', 'MstrPnet','Acc','VelAcq']
    # data = dataFull[usefulData]
    
    # cof = dataFull['Fire_UsedCOF']
    # data = data[data['Time'] < tFinal]
    # cof = cof.iloc[0:data.shape[0]]
    
    data_store.append(dataFull)
    # cof_store.append(cof)
result = pd.concat(data_store)


usefulData  = ['Time','MasterValveCmdCorrected', 'MstrPnet','Acc','VelAcq']
usefulDataFull = usefulData + ['Fire_UsedCOF']
data = result[usefulData]
cof = result['Fire_UsedCOF']
#%%

 
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
# from sklearn.metrics import accuracy_score

data_train, data_test, cof_train, cof_test = train_test_split(data, cof)
reg = LinearRegression().fit(data_train, cof_train)
cof_pred = reg.predict(data_test)
# plt.gca().set_aspect('equal')

# plt.plot(cof_pred, cof_test, 'o')
# plt.xlabel('true COF')
# plt.ylabel('Predectied COF')
# plt.xlim([0,.5])
# plt.ylim([0,.5])
#%%
# accuracy_score(np.array(cof_test), cof_pred)

#Perhaps save this big boi to be able to come back to? Look into efficiency 
# List important headers and only take those for analysis



### Some ideas ###

#svd on all of the data? What kind of energy spectrum do we get?

u, s, v = np.linalg.svd(data - np.mean(data), full_matrices=False)
# Plot matrix of the relevant params?  

# pd.plotting.scatter_matrix(data, diagonal='kde')
# Linear model to predict COF
# Ml to predict COF
# LSTM or other RNN?  Am I just dreaming at this point?

### How might we use this analysis, code or other for future projects.  11A? ###
# Is there some OOP hiding?  Test object?  THis may be beyond my abilities
