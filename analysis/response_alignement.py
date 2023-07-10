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

def data_loader(exp_file, data_loc): #response arrays should have a shape of cells*stims*time
    responses = [f.path for f in os.scandir(exp_file) if f.path.endswith('.npy')]
    thresholds = [f.path for f in os.scandir(exp_file) if f.path.endswith('THRESHOLD.csv')]
    
    resp_db = {}
    thresh_db = {}
    for response in responses:
        
        resp_db[response[response.rfind("A"):response.rfind("df")-1]] = np.load(response, allow_pickle=True)
    
    for threshold in thresholds:
        try:
            df = pd.read_csv(threshold, header = None)
        except pd.errors.EmptyDataError:
            pass
        
        if df.shape == (1,5):
            data = df.iloc[0,:].to_numpy()
            thresh_db[threshold[threshold.rfind("A"):threshold.rfind("df")-1]] = data
            
        else:
            data = df.iloc[:,1:].to_numpy()
            thresh_db[threshold[threshold.rfind("A"):threshold.rfind("df")-1]] = data
            
    #remove cell responses where the max is less than the threshold        
    if data_loc == "neuropil":
        for k, v in resp_db.items():
            for i in range(len(v)):
                if v[i].max() < thresh_db[k][i]:
                    resp_db[k][i] *= np.nan
    elif data_loc == 'tectum':
        for k, v in resp_db.items():
            for w, x in enumerate(v):
                for y, z in enumerate(x):
                    if z.max() < thresh_db[k][w][y]:
                        resp_db[k][w][y] *= np.nan

    return responses, thresholds, resp_db, thresh_db

def mean_stack(stacked): 
            
    meaned = np.nanmean(stacked, axis=0, dtype = np.float64)
    return  meaned

def plot_averaged(meaned, title):
    
    plt.figure(figsize=(20,6))
    position = 150
    for i in meaned:
        position+=1
        plt.subplot(position)
        plt.suptitle(title)
        plt.ylim(0,1.5)
        plt.plot(i)
    plt.show()

def group_by_treatment(treatments, resp_db):
    temp ={}
    for key, values in treatments.items():
        temp[key] ={}
        for value in values:
            kl=[]
            vl=[]
            for k, v in resp_db.items():
                if value in k:
                    kl.append(k)
                    vl.append(v)
            kv = dict(zip(kl,vl))
            temp[key][value] = kv
    return temp  


#Set parameters and run code
exp_file = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/neuropil activity/data-peaks/"
data_loc = "neuropil"
treatments = {"cap_and_train": ["21sep22", "22jun22", "26may22", "13aug22"], "cap_notrain": ["01apr23", "02feb23", "30mar23"],
"nocap_train": ["07sep22", "15jun22", "20jul22", "27apr22"]}
time_list=['min15', 'min45', 'min60', 'min80', 'min100']

responses, thresholds, resp_db, thresh_db = data_loader(exp_file, data_loc)
grouped_by_treatment = group_by_treatment(treatments, resp_db)

if data_loc == "neuropil":
    base_min = {}
    base_max = {}
    ready_to_plot ={}
    for i in time_list:
        ready_to_plot[i] ={}
        for treatment, animals in grouped_by_treatment.items():
            if treatment != "cap_notrain":
                to_stack = []
                for animal, times in animals.items():    
                    for time, trace in times.items():    
                        if i in time:
                            to_stack.append(trace)
                stacked = np.stack(to_stack)
                meaned = mean_stack(stacked)
                meaned_again = mean_stack(meaned)
                if i == 'min15':
                    base_min[treatment] = meaned_again.min()
                    base_max[treatment] = meaned_again.max()
                normalized = (meaned_again - base_min[treatment])/ (base_max[treatment] - base_min[treatment])
                title = treatment + i
                plt.title(title)
                plt.plot(normalized)
                plt.ylim(0,1.1)
        plt.show()
            
            #plot_averaged(meaned_again, title)

elif data_loc == "tectum":
    base_min = {}
    base_max = {}
    ready_to_plot = {}
    for i in time_list:
        ready_to_plot[i]={}
        for treatment, animals in grouped_by_treatment.items():
            if treatment != "cap_notrain":
                to_stack = []
                for animal, times in animals.items():    
                    for time, cells in times.items():
                        for traces in cells:
                            if i in time:
                                to_stack.append(traces)
                                
                stacked = np.stack(to_stack)
                meaned = mean_stack(stacked)
                meaned_again = mean_stack(meaned)
                if i == 'min15':
                    base_min[treatment] = meaned_again.min()
                    base_max[treatment] = meaned_again.max()
                normalized = (meaned_again - base_min[treatment])/ (base_max[treatment] - base_min[treatment])
                ready_to_plot[i][treatment]=normalized
                
                title = treatment + i
                plt.title(title)
                plt.plot(normalized)
                plt.ylim(0,1.5)
        
        plt.show()
                #plot_averaged(meaned, title)
    
    
    
    
    
    