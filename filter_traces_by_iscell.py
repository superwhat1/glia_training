# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 13:09:22 2022

@author: BioCraze
"""

import numpy as np
import os
import glob

#Set folder path and then use it to open the iscell with shape: two columns, first contains 0 or 1 depending on if ROI is a cell and second the probability of that ROI of being a cell according to the suite2p classifier. 
#And to open F.npy containing an array of fluorescent trace arrays
rootdir = 'C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia_training\\data\\training_19may23'
files =[f.path[:f.path.rfind('\\')+1] for i in glob.glob(f'{rootdir}/*/**/***/',recursive=True) for f in os.scandir(i) if f.path.endswith('iscell.npy')]

for f in files:
    is_cell = np.load(f + 'iscell.npy', allow_pickle=True)
    F = np.load(f + 'F.npy',allow_pickle=True)

    #create array of traces containing only the ROIs that are cells according to iscell 
    F_fltrd = np.array([F[i] for i in range(len(is_cell)) if is_cell[i][0] == True])


    #write data to csv files
    #np.savetxt(f + f[f.find("A"):-1] +'_iscell_downsampled_traces.csv', downsample(F_fltrd, True), delimiter = ", ", fmt = "% 1.4f")

    np.savetxt(f[:f.find("a\\")+1] +f[f.find("\A"):f.find("\\s")] +'_neuron_traces.csv', F_fltrd, delimiter=", ", fmt="% 1.4f")

'''
for i in range(len(is_cell)):
    name = "ROI_" + str(i)
    if is_cell[i][0]==True:
        is_cell[i]=F[i]
'''