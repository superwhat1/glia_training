import os

import csv

from sklearn.metrics import auc

import numpy as np

import statistics

import scipy.ndimage as sn

def find_peaks_in_data(data, blurred, threshold_multiplier):
    
    peaks = [[] for i in range(len(data))]

    times = [[] for i in range(len(data))]

    area = [[] for i in range(len(data))]
    
    responses = [[] for i in range(len(data))]
    
    thresholds = [[] for i in range(len(data))]

    #search_extension = 50
    
    start_of_window = -1
    end_of_window = -1
    
    ''' scan over gaussian blurred data to find response windows then pull responses from raw data'''
    
    for row_index, row in enumerate(blurred):

        threshold_to_start = (threshold_multiplier * statistics.stdev(row)) + statistics.mean(row)
        threshold_to_end = (threshold_multiplier * statistics.stdev(row)) + statistics.mean(row)
        print(threshold_to_start)
        for value_index, value in enumerate(row):

            if start_of_window == -1: # If we don't have a starting point yet

                if value >= threshold_to_start:
                    
                    start_of_window = value_index
                    #start_of_window_d = row.index(min(row[value_index - search_extension:value_index])) if value_index > search_extension else row.index(min(row[:value_index+1]))

            elif end_of_window == -1: # We have a starting point, now we're looking for the end of our window

                if value <= threshold_to_end:
                    
                    end_of_window = value_index
                    #end_of_window_d = row.index(min(row[value_index:value_index + search_extension])) if value_index < len(row) - search_extension else row.index(min(row[value_index-1:]))

            else: # We have a start and an end - find the max in this range
                
                if len(data[row_index][start_of_window:end_of_window]) > 1:

                    peaks[row_index].append(max(data[row_index][start_of_window:end_of_window]))
                
                    area[row_index].append(auc(list(range(start_of_window, end_of_window)), list(data[row_index][start_of_window:end_of_window])))
                
                    times[row_index].append(data[row_index].index(max(data[row_index][start_of_window:end_of_window])) + 1)
                    
                    responses[row_index].append(data[row_index][start_of_window:end_of_window])
                    
                    try:
                        thresholds.append(max(data[row_index][start_of_window-100:start_of_window]))
                    except:
                        thresholds.append(max(data[row_index][end_of_window:end_of_window+100]))

                start_of_window = -1

                end_of_window = -1
    return area, peaks, times, responses, thresholds



def read_in_csv(file):

    with open(file) as fn:

        data = []

        reader = csv.reader(fn)

        for row in reader:

            float_row = [float(numeric_string) for numeric_string in row[1:]]

            data.append(float_row)

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


def output_new_csv(area, peaks, times, responses, thresholds, file):


    for index, row in enumerate(area):

        row.insert(0, f"trace_{index+1}")
        

    for index, row in enumerate(peaks):

        row.insert(0, f"trace_{index+1}")


    for index, row in enumerate(times):

        row.insert(0, f"trace_{index+1}")
        

    for index, row in enumerate(thresholds):
    
        row.insert(0, f"trace_{index+1}")
        

    if not os.path.exists('data-peaks/'):

        os.makedirs('data-peaks/')
    

    with open('data-peaks/' + file[file.find("F/")+2:-4] + '_AREA.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in area:

            writer.writerows([row])


    with open('data-peaks/' + file[file.find("F/")+2:-4] + '_PEAKS.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in peaks:

            writer.writerows([row])


    with open('data-peaks/' + file[file.find("F/")+2:-4] + '_TIMES.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in times:

            writer.writerows([row])


    with open('data-peaks/' + file[file.find("A"):-4] + '_THRESHOLD.csv', 'w', newline='') as fn:
        
        writer = csv.writer(fn)

        for row in thresholds:

            writer.writerows([row])
            
            
    with open('data-peaks/' + file[file.find("A"):-4] + '_RESPONSES.npy', 'wb') as fn:
        
        np.save(fn, np.array(responses), allow_pickle=True)


def main(threshold_multiplier):

    files = [i for i in os.listdir() if i.endswith('df.csv')]
    for file in files:
    
        data = read_in_csv(file)
    
        blurred =blur(data, sigma)

        area, peaks, times, responses, thresholds = find_peaks_in_data(data, blurred, threshold_multiplier)

        output_new_csv(area, peaks, times,responses, thresholds, file)

#Set change in fluorescence threshold to use for identifying responses.
threshold_multiplier = 2
#Set amount to blur
sigma = 5

main(threshold_multiplier)