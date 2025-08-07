# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 20:00:18 2023

@author: David
"""

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt


data  = np.loadtxt("C:/Users/David/Documents/Ruthazer lab/glia_training/analysis/neuropil_raw_traces.csv", delimiter = ",")

SR = 15 #Sampling Rate of the microscope
D = 300#duration of imaging session in seconds
SF = 1/60 #Frequency of stimuli presentation

x = range(len(data[0]) )
widths = np.arange(1,12,2)

peakind = signal.find_peaks_cwt(data[3], widths = widths, noise_perc=10, min_snr=1)

yp= data[3][peakind]

plt.plot(x,data[3])
plt.plot(peakind, yp)


'''
for i,h in enumerate(data):
    plt.plot(x, data[i])
'''

t, dt = np.linspace(0, 1, 200, retstep=True)

fs = 1/dt

w = 6

freq = np.linspace(1, fs/2, 20)

widths = w*fs / (2*freq*np.pi)

for i in widths:
    wavelet = signal.morlet2(100, i, w = w)
    plt.plot(abs(wavelet), label = int(i))
plt.legend(prop={'size': 6})
plt.show()
