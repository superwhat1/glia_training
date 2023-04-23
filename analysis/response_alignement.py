# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 14:59:31 2023

@author: David
"""

import numpy as np
import csv
import matplotlib.pyplot as plt

data = []

reader = csv.reader(open("C:/Users/David/Documents/Ruthazer lab/glia_training/analysis/neuropil_raw_traces.csv"))

for row in reader:

    float_row = [float(numeric_string) for numeric_string in row]

    data.append(float_row)
    
            

responses = [[] for i in range(len(data))]

first_stim = 380


for row_index, row in enumerate(data):
    for i in range(5):
        window = first_stim + i*915
        if window - 250 < 0:
            left = 0
            right =window + 250
        elif window + 250 > 4500:
            left = window - 250
            right = 4500
        else:
            left = window - 250
            right = window + 250
            
        responses[row_index].append(data[row_index][left:right])

a_responses = np.array(responses)

target_idx = 250


for sessions, responses in enumerate(a_responses[:3]):
    for response, time in enumerate(responses):
        plt.plot(np.roll(time, target_idx - np.argmax(time)))      
plt.show()


