import os

import csv

from sklearn.metrics import auc

import numpy as np

import statistics

import scipy.ndimage as sn

import pandas as pd

def find_peaks_in_data(data, blurred, threshold_multiplier):
    
    peaks = []

    times = []

    area = []
    
    threshold = []
    
    responses = []
    
    #search_extension = 50
    
    start_of_window = -1
    end_of_window = -1
    
    ''' scan over gaussian blurred data to find response windows then pull responses from raw data'''
    threshold_to_start = (threshold_multiplier * statistics.stdev(blurred)) + statistics.mean(blurred)
    threshold_to_end = (threshold_multiplier * statistics.stdev(blurred)) + statistics.mean(blurred)
    print(threshold_to_start)
    
    for value_index, value in enumerate(blurred):

        if start_of_window == -1: # If we don't have a starting point yet
    
            if value >= threshold_to_start:
                
                start_of_window = value_index
                #start_of_window_d = row.index(min(row[value_index - search_extension:value_index])) if value_index > search_extension else row.index(min(row[:value_index+1]))
    
        elif end_of_window == -1: # We have a starting point, now we're looking for the end of our window
    
            if value <= threshold_to_end:
                
                end_of_window = value_index
                #end_of_window_d = row.index(min(row[value_index:value_index + search_extension])) if value_index < len(row) - search_extension else row.index(min(row[value_index-1:]))
    
        else: # We have a start and an end - find the max in this range
            
            if len(data[start_of_window:end_of_window]) > 1:
    
                peaks.append(max(data[start_of_window:end_of_window]))
            
                area.append(auc(list(range(start_of_window, end_of_window)), list(data[start_of_window:end_of_window])))
            
                times.append(data.index(max(data[start_of_window:end_of_window])) + 1)
                
                responses.append(data[start_of_window:end_of_window])
                try:
                    threshold.append(max(data[start_of_window-100:start_of_window]))
                except:
                    threshold.append(max(data[end_of_window:end_of_window+100]))
                    
    
            start_of_window = -1
    
            end_of_window = -1
    return area, peaks, times, threshold, responses



def read_in_csv(file):

    df = pd.read_csv(file, header = None)
    data = df.iloc[:,0].tolist()

    return data


def blur(data: list, sigma: int) -> list:
    """
    Returns a blurred trace from 'data'.

    Parameters
    ----------
    data: cells, f_intensity
        The fluorescence intensity traces of the recorded active cells to blur
        
    sigma: scalar or sequence of scalars
        Standard deviation for Gaussian kernel. The standard deviations of the Gaussian filter are given for each axis as a sequence, or as a single number, in which case it is equal for all axes.

    Returns
    -------
    F_blur:
        The blurred fluorenscence traces
   """
    nd_data = np.array(data)

    #create array that has m = cell rows and n = time/2 columns.
    F_blur = sn.gaussian_filter1d(nd_data, sigma=sigma)
        
    return F_blur.tolist()


def output_new_csv(area, peaks, times, threshold, responses, file):
        

    if not os.path.exists('data-peaks/'):

        os.makedirs('data-peaks/')
        
        
    with open('data-peaks/' + file[file.find("A"):-4] + '_AREA.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)
        writer.writerow(area)


    with open('data-peaks/' + file[file.find("A"):-4] + '_PEAKS.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        writer.writerow(peaks)


    with open('data-peaks/' + file[file.find("A"):-4] + '_TIMES.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        writer.writerow(times)
            
    with open('data-peaks/' + file[file.find("A"):-4] + '_THRESHOLD.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        writer.writerow(threshold)
    
    with open('data-peaks/' + file[file.find("A"):-4] + '_RESPONSES.npy', 'wb') as fn:
        np.save(fn, np.array(responses), allow_pickle=True)


def main(threshold_multiplier):

    files = [i for i in os.listdir() if i.endswith('df.csv')]
    for file in files:
    
        data = read_in_csv(file)
    
        blurred =blur(data, sigma)

        area, peaks, times, threshold, responses = find_peaks_in_data(data, blurred, threshold_multiplier)

        output_new_csv(area, peaks, times, threshold, responses, file)

#Set change in fluorescence threshold to use for identifying responses.
threshold_multiplier = 2
#Set amount to blur
sigma = 5

main(threshold_multiplier)