import os

import csv

from sklearn.metrics import auc

import pandas as pd

import numpy as np

from statistics import stdev


def find_peaks_in_data(data, threshold_start_multiplier, threshold_end_multiplier, percentile):
    
    peaks = [[] for i in range(len(data))]

    times = [[] for i in range(len(data))]

    area = [[] for i in range(len(data))]

    search_extension = 50
    
    start_of_window = -1
    end_of_window = -1
    
    for row_index, row in enumerate(data):
        
        baseline = np.percentile(row, percentile)
        standard_dev = stdev(row)
        threshold_to_start = baseline + standard_dev * threshold_start_multiplier
        threshold_to_end = baseline + standard_dev * threshold_end_multiplier
        
        for value_index, value in enumerate(row):

            if start_of_window == -1: # If we don't have a starting point yet

                if value >= threshold_to_start:
                    print(value_index)
                    print(threshold_to_start)
                    print(value)
                    start_of_window = row.index(min(row[value_index - search_extension:value_index])) if value_index > search_extension else row.index(min(row[:value_index+1]))

            elif end_of_window == -1: # We have a starting point, now we're looking for the end of our window

                if value <= threshold_to_end:

                    end_of_window = row.index(min(row[value_index:value_index + search_extension])) if value_index < len(row) - search_extension else row.index(min(row[value_index-1:]))

            else: # We have a start and an end - find the max in this range

                if len(data[row_index][start_of_window:end_of_window]) > 1:

                    peaks[row_index].append(max(data[row_index][start_of_window:end_of_window]))
                
                    area[row_index].append(auc(list(range(start_of_window, end_of_window)), list(data[row_index][start_of_window:end_of_window])))
                
                    times[row_index].append(data[row_index].index(max(data[row_index][start_of_window:end_of_window])) + 1)

                start_of_window = -1

                end_of_window = -1

    return area, peaks, times



def read_in_csv(file):

    with open(file) as fn:

        data = []

        reader = csv.reader(fn)

        for row in reader:

            float_row = [float(numeric_string) for numeric_string in row[1:]]

            data.append(float_row)

    return data



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
    
    area, peaks, times = find_peaks_in_data(data, threshold_start_multiplier, threshold_end_multiplier, percentile)

    output_new_csv(area, peaks, times, file)


#Set change in fluorescence threshold to use for identifying responses.

threshold_start_multiplier = 2; threshold_end_multiplier = 2; percentile = 7.5

main(threshold_start_multiplier, threshold_end_multiplier, percentile)
