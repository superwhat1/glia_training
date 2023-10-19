# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 14:51:44 2023

@author: David
"""

from scipy.fft import rfft, rfftfreq, irfft
import numpy as np
import matplotlib.pyplot as plt

data  = np.loadtxt("C:/Users/David/Documents/Ruthazer lab/glia_training/analysis/neuropil_raw_traces.csv", delimiter = ",")

SR = 15 #Sampling Rate of the microscope
D = 300#duration of imaging session in seconds
SF = 1/60 #Frequency of stimuli presentation

x = rfftfreq(SR*D, 1/SR)
y = rfft(data)

fft_plot = plt.plot(x[:300], np.abs(y[3][:300]))

target_idx = np.argwhere(x == SF)[0,0]

y[:, target_idx+30:]=0
    

fltrd_signal = irfft(y)
    

raw_plot = plt.plot(data[3])
fltrd_plot = plt.plot(fltrd_signal[3])

