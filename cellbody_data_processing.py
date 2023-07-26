# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 15:15:43 2023

@author: BioCraze
"""
import numpy as np
import pandas as pd
import os, glob, statistics, csv
from sklearn.metrics import auc


def is_cell_responsive(f):
    
    #Set folder path and then use it to open the iscell with shape: two columns, first contains 0 or 1 depending on if ROI is a cell and second the probability of that ROI of being a cell according to the suite2p classifier. 
    #And to open F.npy containing an array of fluorescent trace arrays

    is_cell = np.load(f + 'iscell.npy', allow_pickle=True)
    F = np.load(f + 'F.npy',allow_pickle=True)

    #create array of traces containing only the ROIs that are cells according to iscell 
    F_fltrd = np.array([F[i] for i in range(len(is_cell)) if is_cell[i][0] == True])


    #write data to csv files
    #np.savetxt(f + f[f.find("A"):-1] +'_iscell_downsampled_traces.csv', downsample(F_fltrd, True), delimiter = ", ", fmt = "% 1.4f")

    np.savetxt(f[:f.find("a\\")+1] +f[f.find("\A"):f.find("\\s")] +'_neuron_traces.csv', F_fltrd, delimiter=", ", fmt="% 1.4f")
    
    return F_fltrd


def calculate_percentiles(data, percentile=7.5):

    window_size = 150
    
    percentiles = [[] for i in range(len(data[0:]))]

    for cell, frames in enumerate(data): # One cell at a time
    

        for frame, value in enumerate(frames): # One frame at a time
        

            # Get the index for window/2 before and after frame

            before = frame - (window_size / 2) if (frame > window_size / 2) else 0

            after = frame + (window_size / 2) if (frame + (window_size / 2)) < (len(frames) - 1) else len(frames) - 1

            #create the sliding window from the values of column with before and after indexes

            window = frames[int(before):int(after)]
            
            #Band filter the window to pull the 10th(exclusive) to 5th percentile(exclusive) values
            #btm_prcntl = list(filter(lambda a: a < np.percentile(window,percentile),window))
            #baseline_band = list(filter(lambda a: a > np.median(btm_prcntl),btm_prcntl))
            
            percentiles[cell].append(np.percentile(window, percentile))

    return percentiles


def deltaF(data, percentiles):

    Fo_stdev = []
    deltaf = []

    count=0
    
    for cell, frames in enumerate(data):
        count+=1
        new_column = []
        for frame, value in enumerate(frames):
            try:
                new_column.append((float(value) - percentiles[cell][frame]) / float(percentiles[cell][frame]))
    
            except IndexError:
                new_column.append((float(value) - percentiles[cell][frame-1]) / float(percentiles[cell][frame-1]))
            
        Fo_stdev.append([ f"trace_{count}",(1.5 * statistics.stdev(percentiles[cell])) / statistics.mean(percentiles[cell])]) #append the delta F/Fo of the 1.5*STDEV of Fo 
        deltaf.append(new_column)

    return Fo_stdev, deltaf
    
    
def find_peaks_in_data(data, stims, animal):


    peaks = [[] for i in range(len(data))]

    times = [[] for i in range(len(data))]

    area = [[] for i in range(len(data))]
    
    responses = [[] for i in range(len(data))]
    
    thresholds = [[] for i in range(len(data))]
    
    a=0
    
    first_stim = stims.at[animal, "first stim"]

  
    for row_index, row in enumerate(data):
        a+=1
        for i in range(5):
            window = first_stim + i*915
            if window - 100 < 0:
                left = 0
                right =window + 350 - first_stim
            elif window + 250 > 4500:
                left = window - 350 + (4500 - window)
                right = 4500
            else:
                left = window - 100
                right = window + 250
            
            area[row_index].append(auc(list(range(left, right)), list(data[row_index][left:right])))
            if left == 0:
                rolling_trgt = 100 - first_stim
                responses[row_index].append(np.roll(np.array(data[row_index][left:right]), rolling_trgt))
            
            elif right == 4500:
                rolling_trgt = 4500 - window - 100
                responses[row_index].append(np.roll(np.array(data[row_index][left:right]), rolling_trgt))
            else:
                responses[row_index].append(data[row_index][left:right])
                
            peaks[row_index].append(max(data[row_index][left:right]))
            
            times[row_index].append(data[row_index].index(max(data[row_index][left:right])) + 1)
            try:
                thresholds[row_index].append(max(data[row_index][left-100:left]))
            except:
                thresholds[row_index].append(max(data[row_index][right:right+100]))

    return area, peaks, times, thresholds, responses
    
    
    
def output_csv(Fo, data, area, peaks, times, thresholds, responses, file):
      
    print("writing")
    if not os.path.exists(file[:file.rfind("/")+1]+'deltaF/'):

        os.makedirs(file[:file.rfind("/")+1]+'deltaF/')



    with open(file[:file.rfind("/")+1] + 'deltaF/' + file[file.rfind("/A")+1:-4] + '_df.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in data:
            writer.writerows([row])
            
            
    with open(file[:file.rfind("/")+1] + 'deltaF/' + file[file.rfind("/A")+1:-4] + '_Fo.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)
        
        writer.writerows(Fo)
        
    if len(area) > 0: 
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
            
            
        with open('data-peaks/' + file[file.find("A"):-4] + '_AREA.csv', 'w', newline='') as fn:
    
            writer = csv.writer(fn)
    
            for row in area:
    
                writer.writerows([row])
    
    
        with open('data-peaks/' + file[file.find("A"):-4] + '_PEAKS.csv', 'w', newline='') as fn:
    
            writer = csv.writer(fn)
    
            for row in peaks:
    
                writer.writerows([row])
    
    
        with open('data-peaks/' + file[file.find("A"):-4] + '_TIMES.csv', 'w', newline='') as fn:
    
            writer = csv.writer(fn)
    
            for row in times:
    
                writer.writerows([row])
                
        with open('data-peaks/' + file[file.find("A"):-4] + '_THRESHOLD.csv', 'w', newline='') as fn:
            
            writer = csv.writer(fn)
    
            for row in thresholds:
    
                writer.writerows([row])
                
        with open('data-peaks/' + file[file.find("A"):-4] + '_RESPONSES.npy', 'wb') as fn:
            
            np.save(fn, np.array(responses), allow_pickle=True)
                



def main():
    rootdir = 'C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia_training\\data\\training_19may23'
    files =[f.path[:f.path.rfind('\\')+1] for i in glob.glob(f'{rootdir}/*/**/***/',recursive=True) for f in os.scandir(i) if f.path.endswith('iscell.npy')]
    
    for file in files:
        #perform deltaF operation
        raw_trace = is_cell_responsive(file)
        percentiles = calculate_percentiles(raw_trace)
        Fo_temp, deltaf = deltaF(raw_trace, percentiles)
            
        #perform response metric operations
        stims_timings = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/summaries/stim_timings.csv"
        stims = pd.read_csv(stims_timings)
        stims = stims.set_index("file")        
        
        animal = file[file.rfind("A"):file.rfind("neu")-1]
        try:
            area, peaks, times, thresholds, responses = find_peaks_in_data(deltaF, stims, animal)
            
        except KeyError:
            print(file)
            pass
        
        output_csv(Fo_temp, deltaf, area, peaks, times, thresholds, responses, file)