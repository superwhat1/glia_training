# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 22:28:37 2023

@author: BioCraze
"""

import numpy as np
from scipy.fft import rfft, rfftfreq, irfft
import matplotlib.pyplot as plt

data = np.loadtxt("C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/neuropil activity/cumulated_neuropil.csv", delimiter=",")

SR = 300 #sampling rate in Hz
D = 15 #Duration in frames

x = rfftfreq(SR*D, 1/SR) #create x axis of length SR*D with values spaced at 1/SR
y = rfft(data) #run the fast fourier transform for data containing real values

fft = plt.plot(x[:300],np.abs(y[3][:300])) #plot the power freqency data

#set all frequencies other than the stimuli(highest amplitude) freqencies to zero 

for i, h in enumerate(y):
    target_idx = np.where(y[i]==y[i].max())[0][0]
    y[i][target_idx+30:] = 0

fltrd_signal = irfft(y) #calculate inverse fft to create a filtered trace

#plot the filtered signals
fltrd_plot = plt.plot(fltrd_signal[4])

#plot the raw signal traces
raw = plt.plot(data[4])
