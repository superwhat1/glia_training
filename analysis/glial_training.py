# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 13:09:22 2022

@author: BioCraze
"""

import numpy as np
import os
import glob


def downsample(data: np.ndarray, taper_edge: bool = True) -> np.ndarray:
    """
    Returns a downsampled trace from 'data', tapering the edges of 'taper_edge' is True.

    Parameters
    ----------
    data: cells, f_intensity
        The fluorescence intensity traces of the recorded active cells to downsample
    taper_edge: bool
        Whether to taper the edges

    Returns
    -------
    F_down:
        The downsampled fluorenscence traces
   """
    cells, f_intensity = data.shape

    # bin along Y
    F_down = np.zeros((cells, int(np.ceil(f_intensity / 2))), 'float32')
    F_down[:, :f_intensity//2] = np.mean([data[:, 0:-1:2], data[:, 1::2]], axis=0)
    
    if f_intensity % 2 == 1:
        F_down[:, -1] = data[:, -1] / 2 if taper_edge else data[:, -1]


    return F_down

#Set folder path and then use it to open the iscell with shape: two columns, first 0 or 1 depending on if ROI is a cell and second the probability of that ROI of being a cell according to the suite2p classifier. 
#And to open F.npy containing an array of fluorescent trace arrays
rootdir = 'C:\\Users\\BioCraze\\Documents\\Ruthazer lab\\glia_training\\data'
files =[f.path[:f.path.rfind('\\')+1] for i in glob.glob(f'{rootdir}/*/**/***/****/*****/',recursive=True) for f in os.scandir(i) if f.path.endswith('tectal_neuron_iscell.npy')]

for f in files:
    is_cell = np.load(f + 'iscell.npy', allow_pickle=True)
    F = np.load(f + 'F.npy',allow_pickle=True)

    #create array of traces containing only the ROIs that are cells according to iscell 
    F_fltrd = np.array([F[i] for i in range(len(is_cell)) if is_cell[i][0] == True])


    #write data to csv files
    #np.savetxt(f + f[f.find("A"):-1] +'_iscell_downsampled_traces.csv', downsample(F_fltrd, True), delimiter = ", ", fmt = "% 1.4f")

    np.savetxt(f[:f.find("\m")] +'_glia_traces.csv', F_fltrd, delimiter=", ", fmt="% 1.4f")

'''
for i in range(len(is_cell)):
    name = "ROI_" + str(i)
    if is_cell[i][0]==True:
        is_cell[i]=F[i]
'''