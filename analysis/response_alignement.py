# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 14:59:31 2023

@author: David
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import os

def data_loader(exp_file): #response arrays should have a shape of cells*stims*time
    files =[f.path for i in glob.glob(f'{exp_file}/*/', recursive=True) for f in os.scandir(i) if f.path.endswith('.npy')]

    db = {}
    for file in files:
        print(file[file.rfind("A"):-4])
        db[file[file.rfind("A"):-4]] = np.load(file, allow_pickle=True)

    return db

def roll_and_average(data): 

    target_idx = int(np.floor(data.shape[-1]/2))

    rolled = np.zeros(np.shape(data))

    for session, responses in enumerate(data):
        for response, time in enumerate(responses):
            rolled[session,response] = np.roll(time, target_idx - np.argmax(time))   
            
    meaned = np.mean(rolled,axis=0)
    return rolled, meaned

def plot_averaged(meaned):

    plt.plot(meaned)

    plt.show()


file = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/neuropil activity/"
db = data_loader(file)
data = db['A2_min100_27apr22_df_RESPONSES']
rolled, meaned = roll_and_average(data)
write()
#plot_averaged(meaned)
