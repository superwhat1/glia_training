# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 02:28:31 2023

@author: BioCraze
"""


import numpy as np
import csv
import scipy.ndimage as sn

def downsample(data):
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
    F_down = sn.gaussian_filter1d(data, 5)


    return F_down

file = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/glia activity/A1_cap_activation_21sep22_glia_traces.csv"
with open(file, newline='') as fn:
    data = []

    reader = csv.reader(fn)

    for row in reader:

        float_row = [float(numeric_string) for numeric_string in row[1:]]

        data.append(float_row)
        
nd = np.array(data)

np.savetxt(file[:file.find("A")] +'_downsampled_traces_gaussian_s5.csv', downsample(nd), delimiter = ", ", fmt = "% 1.4f")