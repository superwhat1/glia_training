import os

import csv

from sklearn.metrics import auc

import numpy as np

from statistics import stdev


def find_peaks_in_data(data, downsampled, threshold_start_multiplier, threshold_end_multiplier, percentile):
    
    peaks = [[] for i in range(len(data))]

    times = [[] for i in range(len(data))]

    area = [[] for i in range(len(data))]

    search_extension = 50
    
    start_of_window_d = -1
    end_of_window_d = -1
    
    ''' convert function to scan over downsampled data to find response windows then pull responses from raw'''
    
    for row_index, row in enumerate(downsampled):
        
        baseline = np.percentile(row, percentile)
        standard_dev = stdev(row)
        threshold_to_start = baseline + standard_dev * threshold_start_multiplier
        threshold_to_end = baseline + standard_dev * threshold_end_multiplier
        
        for value_index, value in enumerate(row):

            if start_of_window_d == -1: # If we don't have a starting point yet

                if value >= threshold_to_start:
                    print(value_index)
                    print(threshold_to_start)
                    print(value)
                    start_of_window_d = row.index(min(row[value_index - search_extension:value_index])) if value_index > search_extension else row.index(min(row[:value_index+1]))

            elif end_of_window_d == -1: # We have a starting point, now we're looking for the end of our window

                if value <= threshold_to_end:

                    end_of_window_d = row.index(min(row[value_index:value_index + search_extension])) if value_index < len(row) - search_extension else row.index(min(row[value_index-1:]))

            else: # We have a start and an end - find the max in this range
                
                start_of_window = start_of_window_d * 2
                end_of_window = end_of_window_d * 2
                
                if len(data[row_index][start_of_window:end_of_window]) > 1:

                    peaks[row_index].append(max(data[row_index][start_of_window:end_of_window]))
                
                    area[row_index].append(auc(list(range(start_of_window, end_of_window)), list(data[row_index][start_of_window:end_of_window])))
                
                    times[row_index].append(data[row_index].index(max(data[row_index][start_of_window:end_of_window])) + 1)

                start_of_window_d = -1

                end_of_window_d = -1

    return area, peaks, times



def read_in_csv(file):

    with open(file) as fn:

        data = []

        reader = csv.reader(fn)

        for row in reader:

            float_row = [float(numeric_string) for numeric_string in row[1:]]

            data.append(float_row)

    return data


def downsample(data: list, taper_edge: bool = True) -> list:
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

    #create array that has m = cell rows and n = time/2 columns.
    F_down = np.zeros((cells, int(np.ceil(f_intensity / 2))), 'float32')
    
    #create the downsampled array by taking the mean of a list made of the first and the second number in the row then repeat until the end of the row is reached.
    F_down[:, :f_intensity//2] = np.mean([data[:, 0:-1:2], data[:, 1::2]], axis=0)
    
    #check to see if the length of the rows are even, if not devide the last value in the row by 2 if taper_edge = True
    if f_intensity % 2 == 1:
        F_down[:, -1] = data[:, -1] / 2 if taper_edge else data[:, -1]


    return F_down()


def output_new_csv(area, peaks, times, file):


    for index, row in enumerate(area):

        row.insert(0, f"trace_{index+1}")
        

    for index, row in enumerate(peaks):

        row.insert(0, f"trace_{index+1}")



    for index, row in enumerate(times):

        row.insert(0, f"trace_{index+1}")



    if not os.path.exists('data-peaks/'):

        os.makedirs('data-peaks/')
        

    with open('data-peaks/' + file[file.find("A"):-7] + '_AREA.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in area:

            writer.writerows([row])


    with open('data-peaks/' + file[file.find("A"):-7] + '_PEAKS.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in peaks:

            writer.writerows([row])



    with open('data-peaks/' + file[file.find("A"):-7] + '_TIMES.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in times:

            writer.writerows([row])


def main(threshold_start_multiplier, threshold_end_multiplier, percnetile):

#    files = [i for i in os.listdir() if i.endswith('.csv')]
#    for file in files:

    file = "C:/Users/BioCraze/Documents/Ruthazer lab/glial training/cumulated_neuropil.csv"
    
    data = read_in_csv(file)
    
    downsampled =downsample(data)
    
    area, peaks, times = find_peaks_in_data(data,downsampled, threshold_start_multiplier, threshold_end_multiplier, percentile)

    output_new_csv(area, peaks, times, file)


#Set change in fluorescence threshold to use for identifying responses.

threshold_start_multiplier = 2; threshold_end_multiplier = 2; percentile = 7.5

main(threshold_start_multiplier, threshold_end_multiplier, percentile)