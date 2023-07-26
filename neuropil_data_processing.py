# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 15:15:09 2023

@author: BioCraze
"""

import numpy as np
import pandas as pd
import os, glob, statistics, csv

def read_in_neuropil_csv(file):
    raw = pd.read_csv(file)
    Ftrace = raw["Mean"]
    data = Ftrace.to_list()
    return data

def calculate_neuropil_percentiles(data, percentile=7.5):

    window_size = 150
    
    baselines = []

    for frames, value in enumerate(data): # One column at a time

    # Get the index for window/2 before and after column_index

        before = frames - (window_size / 2) if (frames > window_size / 2) else 0

        after = frames + (window_size / 2) if (frames + (window_size / 2)) < (len(data) - 1) else len(data) - 1

    #create the sliding window from the values of column with before and after indexes

        window = data[int(before):int(after)]
            
        baselines.append(np.percentile(window, percentile))

    return baselines


def neuropil_deltaF(data, percentiles):

    deltaf = []

    count=0

    for frames, value in enumerate(data):
        count+=1
        ratioed = []

        try:
            ratioed.append((float(value) - percentiles[frames]) / float(percentiles[frames]))

        except IndexError:
            ratioed.append((float(value) - percentiles[frames-1]) / float(percentiles[frames-1]))
         
        deltaf.append(ratioed)

    Fo_stdev =[ f"trace_{count}",(1.5 * statistics.stdev(percentiles)) / statistics.mean(percentiles)] #calculate the delta F/Fo of the 1.5*STDEV of Fo 
    return Fo_stdev, deltaf


def output_neuropil_csv(Fo, data, file):
      
    print("writing")
    if not os.path.exists(file[:file.rfind("/")+1]+'deltaF/'):

        os.makedirs(file[:file.rfind("/")+1]+'deltaF/')



    with open(file[:file.rfind("/")+1] + 'deltaF/' + file[file.rfind("y/")+2:-4] + '_df.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in data:
            writer.writerows([row])
            
            
    with open(file[:file.rfind("/")+1] + 'deltaF/' + file[file.rfind("y/")+2:-4] + '_Fo.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)
        
        writer.writerows(Fo)
    
    
def main():

    files = [i.path for i in os.scandir("C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/neuropil activity/") if i.path.endswith('.csv')]

    for file in files:
        
        data = read_in_neuropil_csv(file)
        
        percentiles = calculate_neuropil_percentiles(data)
        
        Fo_temp, deltaf = neuropil_deltaF(data, percentiles)
        
        output_neuropil_csv(Fo_temp, deltaf, file)
    
