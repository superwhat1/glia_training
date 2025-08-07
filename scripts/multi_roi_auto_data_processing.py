# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 15:15:43 2023

@author: BioCraze
"""
#Import packages needed to install packages in case they are not installed.
import sys, subprocess

#Attempt to import packages and if they are not installed, install them.
try:
    import numpy as np
    import pandas as pd
    import os, glob, statistics, csv, re
    from sklearn.metrics import auc
    import matplotlib.pyplot as plt
    import pickle
    from datetime import date
    from copy import deepcopy
    import scipy.ndimage as sn
    import math
    
except ImportError as e:
    package = re.search("\'.+\'", str(e)).group()[1:-1]
    print(package + " not installed")
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])


def is_cell_responsive(file, recording_name, output_dir):
    
    #Set folder path and then use it to open the iscell with shape: two columns, first contains 0 or 1 depending on if ROI is a cell and second the probability of that ROI of being a cell according to the suite2p classifier. 
    #And to open F.npy containing an array of fluorescent trace arrays
    is_cell = np.load(file + 'glia_iscell.npy', allow_pickle=True)
    F = np.load(file + 'F.npy',allow_pickle=True)

    #create array of traces containing only the ROIs that are cells according to iscell 
    F_fltrd = np.array([F[i] for i in range(len(is_cell)) if is_cell[i][0] == True])

    #write data to csv files
    np.savetxt(output_dir + recording_name +'_glia_traces.csv', F_fltrd, delimiter=", ", fmt="% 1.4f")
    
    return F_fltrd


def calculate_baselines(data, percentile=7.5):

    window_size = 150
    baselines = [[] for i in range(len(data[0:]))]

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
            
            baselines[cell].append(np.percentile(window, percentile))

    return baselines

def deltaF(data, baselines):

    deltaf = []
    count=0
    
    for cell, frames in enumerate(data):
        
        count+=1
        new_trace = []
        
        for frame, value in enumerate(frames):
            try:
                new_trace.append((float(value) - baselines[cell][frame]) / float(baselines[cell][frame]))
    
            except IndexError:
                new_trace.append((float(value) - baselines[cell][frame-1]) / float(baselines[cell][frame-1]))
            
        deltaf.append(new_trace)

    return  deltaf
    
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


def find_peaks_in_data(data, blurred, recording_name, threshold_multiplier):

    peaks = [[] for i in range(len(data))]
    times = [[] for i in range(len(data))]
    area = [[] for i in range(len(data))]
    responses = [[] for i in range(len(data))]
    thresholds = [[] for i in range(len(data))]
    start_of_window = -1
    end_of_window = -1
    
    ''' scan over gaussian blurred data to find response windows then pull responses from raw data'''
    for row_index, row in enumerate(blurred):
    
        threshold_to_start = (threshold_multiplier * statistics.stdev(row)) + statistics.mean(row)
        threshold_to_end = (threshold_multiplier * statistics.stdev(row)) + statistics.mean(row)

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
                
                if len(data[row_index][start_of_window:end_of_window]) > 15:
    
                    peaks[row_index].append(max(data[row_index][start_of_window:end_of_window]))
                
                    area[row_index].append(auc(list(range(start_of_window, end_of_window)), list(data[row_index][start_of_window:end_of_window])))
                
                    times[row_index].append(data[row_index].index(max(data[row_index][start_of_window:end_of_window])) + 1)
                    
                    if start_of_window-15<=0:
                        responses[row_index].append(data[row_index][0 : end_of_window+15 + (15-start_of_window) ])
                        
                    elif end_of_window+15>=(len(data[row_index])-1):
                        responses[row_index].append(data[row_index][start_of_window-15 - (15-len(data[row_index])-end_of_window):len(data[row_index])])
                        
                    else:
                        responses[row_index].append(data[row_index][start_of_window-15:end_of_window+15])
                
                    try:
                        thresholds[row_index].append(max(data[row_index][start_of_window-30:start_of_window-15]))
                    except:
                        thresholds[row_index].append(max(data[row_index][end_of_window+15:end_of_window+30]))

                start_of_window = -1
                end_of_window = -1 

    return area, peaks, times, thresholds, responses

    
def output_csv(baselines, deltaf, area, peaks, times, thresholds, responses, output_dir, recording_name):
      
    print("writing")
    if not os.path.exists(output_dir + 'deltaF/'):
        os.makedirs(output_dir + 'deltaF/')

    with open(output_dir + 'deltaF/' + recording_name + '_df.csv', 'w', newline='') as fn:
        writer = csv.writer(fn)
        for row in deltaf:
            writer.writerows([row])
            
    with open(output_dir + 'deltaF/' + recording_name + '_Fo.csv', 'w', newline='') as fn:
        writer = csv.writer(fn)
        for row in baselines:
            writer.writerows([row])
        
    if len(area) > 0:
        area_c = deepcopy(area)
        for index, row in enumerate(area_c):
            row.insert(0, f"trace_{index+1}")
        
        peaks_c = deepcopy(peaks)
        for index, row in enumerate(peaks_c):
            row.insert(0, f"trace_{index+1}")
            
        times_c = deepcopy(times) 
        for index, row in enumerate(times_c):
            row.insert(0, f"trace_{index+1}")
            
        thresholds_c = deepcopy(thresholds)
        for index, row in enumerate(thresholds_c):
            row.insert(0, f"trace_{index+1}")
            
        if not os.path.exists(output_dir + 'data-peaks/'):
            os.makedirs(output_dir + 'data-peaks/')
            
        with open(output_dir + 'data-peaks/' + recording_name + '_AREA.csv', 'w', newline='') as fn:
            writer = csv.writer(fn)
            for row in area_c:
                writer.writerows([row])
    
        with open(output_dir + 'data-peaks/' + recording_name + '_PEAKS.csv', 'w', newline='') as fn:
            writer = csv.writer(fn)
            for row in peaks_c:
                writer.writerows([row])
    
        with open(output_dir + 'data-peaks/' + recording_name + '_TIMES.csv', 'w', newline='') as fn:
            writer = csv.writer(fn)
            for row in times_c:
                writer.writerows([row])
                
        with open(output_dir + 'data-peaks/' + recording_name + '_THRESHOLD.csv', 'w', newline='') as fn:
            writer = csv.writer(fn)
            for row in thresholds_c:
                writer.writerows([row])
                
        with open(output_dir + 'data-peaks/' + recording_name + '_RESPONSES.npy', 'wb') as fn:
            np.save(fn, np.array(responses, dtype=object), allow_pickle=True)
                

def mean_stack(stacked):
    
    meaned = np.nanmean(stacked, axis=0, dtype = np.float64)
    
    return  meaned

'''
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
'''

def pad(to_pad):
    longest_resp_len = 0
    most_responses =0
    
    #find the longest response and the cell with the most responses       
    for recording, cells in to_pad.items():
        for cell, responses in enumerate(cells):
            if len(responses)>most_responses:
                    most_responses=len(responses)
            for response, trace in enumerate(responses):
                if len(trace)>longest_resp_len:
                    longest_resp_len=len(trace)
    
    for recording, cells in to_pad.items():
        #print(recording)
        for cell, responses in enumerate(cells):
            #print(cell)
            for response, trace in enumerate(responses):
                trace_pad_width = longest_resp_len-len(trace)
                if trace_pad_width%2 != 0:
                    to_pad[recording][cell][response] = np.pad(np.array(trace),pad_width = (math.floor(trace_pad_width/2), math.ceil(trace_pad_width/2)), mode = 'edge').tolist()
                else:
                    to_pad[recording][cell][response] = np.pad(np.array(trace),pad_width = int(trace_pad_width/2), mode = 'edge').tolist()
                    
            cell_pad_width = most_responses-len(responses)
            if cell_pad_width%2 != 0:
                to_pad[recording][cell] = np.pad(np.array(responses),pad_width = ( (math.floor(cell_pad_width/2), math.ceil(cell_pad_width/2)), (0,0) ), mode = 'constant', constant_values = np.nan)
            else:
                #print(len(responses))
                to_pad[recording][cell] = np.pad(np.array(responses), pad_width = ( ( int(cell_pad_width/2), int(cell_pad_width/2) ), (0,0) ), mode = 'constant', constant_values = np.nan)
                
    for recording, cells in to_pad.items():
        to_pad[recording]=np.array(to_pad[recording])
    
    return to_pad
    

def thresh_filter_responses(resp, thresh): #response arrays should have a shape of cells*stims*time
    #remove cell responses where the max is less than the threshold        
    for recording, cells in resp.items():
        resp[recording]= [cell for cell in cells if cell]
    for rec_thresh, cells_thresh in thresh.items():
        thresh[rec_thresh]=[cell_thresholds for cell_thresholds in cells_thresh if cell_thresholds]
    
    for recording, cells in resp.items():
        for cell, responses in enumerate(cells): 
            resp[recording][cell] = [trace for trace in responses if trace]
    for rec_thresh, cells_thresh in thresh.items():
        for cell_thresh, resp_thresh in enumerate(cells_thresh):
            thresh[rec_thresh][cell_thresh] = [threshold for threshold in resp_thresh if threshold]
                
    for recording, cells in resp.items():
        for cell, responses in enumerate(cells):
            for response, trace in enumerate(responses):
                if max(trace, default = 0) < thresh[recording][cell][response]:
                    resp[recording][cell][response] = np.array(resp[recording][cell][response])
                    resp[recording][cell][response] *= np.nan
                    resp[recording][cell][response].tolist()
                

    fltrd_resp_db = pad(resp)
    
    return fltrd_resp_db


def group_by_treatment(resp_db_to_group, thresh_db_to_group):
    treatments = {"cap_and_train": ["21sep22", "22jun22", "26may22", "13aug22"], "cap_notrain": ["01apr23", "02feb23", "30mar23"],
    "nocap_train": ["07sep22", "15jun22", "20jul22", "27apr22"], "nocap_notrain":["19may23", "23aug23"]}
    
    fltrd_resp_db = thresh_filter_responses(resp_db_to_group, thresh_db_to_group)
    
    grouped_by_treatment ={}
    
    for treatment, experiments in treatments.items():
        grouped_by_treatment[treatment] ={}
        for experiment in experiments:
            kl=[]
            vl=[]
            for exp, traces in fltrd_resp_db.items():
                if experiment in exp:
                    kl.append(exp)
                    vl.append(traces)
            kv = dict(zip(kl,vl))
            grouped_by_treatment[treatment][experiment] = kv
            
    return grouped_by_treatment  


#Function for ploting aligned and averaged responses
def plot_averaged(grouped_by_treatment, time_list):
    
    base_min = {}
    base_max = {}
    ready_to_plot = {}
    
    for i in time_list:
        ready_to_plot[i] = {}
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


def make_plots(grouped_by_treatment):#Function for ploting aligned and averaged responses

    #Prepare variables needed to plot average traces
    time_list=['min15', 'min45', 'min60', 'min80', 'min100']
    
    plot_averaged(grouped_by_treatment, time_list)
    
def process_from_raw_traces(input_dir, output_dir):
    
    files =[f.path[:f.path.rfind('\\')+1] for i in glob.glob(f'{input_dir}/*/**/***/',recursive=True) for f in os.scandir(i) if f.path.endswith('glia_iscell.npy') and "capapplication" not in f.path]
    
    resp_db = {}
    thresh_db = {}
    
    for file in files:
        #find the name of the recording that starts with an A and a number followed by _min and ends with a date in the format of ddmmmyy with dd and yy as ints and mmm as string
        recording_name = re.search("A\d{1}_min.+\d{2}\D{3}\d{2}", file).group()

        #perform deltaF operation
        #try:
        print("Normalizing " + file)
        raw_trace = is_cell_responsive(file, recording_name, output_dir)
        baselines = calculate_baselines(raw_trace)
        deltaf = deltaF(raw_trace, baselines)
        #except Exception:
         #   print(file + " contains no traces, so it was skipped!")
          #  pass
        
        #perform response metric operations 
        try:
            print("Extracting response metrics of " + file)
            blurred = blur(deltaf, sigma = 5)
            
            area, peaks, times, thresholds, responses = find_peaks_in_data(deltaf, blurred, recording_name, threshold_multiplier = 2)
            
            #Agregate the responses and the thresholds then filter them and group them by experimental condition
            resp_db[recording_name] = responses
            thresh_db[recording_name] = thresholds
           
            #Save the data to files
            output_csv(baselines, deltaf, area, peaks, times, thresholds, responses, output_dir, recording_name)
            
        except KeyError:
            print(file + " contains no responses, so it was skipped!")
            pass
        
    grouped_by_treatment = group_by_treatment(resp_db, thresh_db)
        
    
    #Write the grouped data to a txt file for use at another time
    with open(output_dir + 'grouped_glia_responses_by_treatment' + str(date.today()) + '.pkl', 'wb') as of:
        pickle.dump(grouped_by_treatment, of)
        
            
def process_from_responses(input_dir, output_dir):

    response_files =[f.path  for f in os.scandir(input_dir) if f.path.endswith('RESPONSES.npy')]  
    resp_db = {}
    thresh_db = {}
    
    for file in response_files:
        #find the name of the recording that starts with an A and ends with a date in the format of ddmmmyy with dd and yy as ints and mmm as string
        recording_name = re.search("A\d{1}_min.+\d{2}\D{3}\d{2}", file).group()
        responses = np.load(file, allow_pickle=True)
        thresholds = np.array(pd.read_csv(file[:file.rfind('_')] + '_THRESHOLD.csv', header=None).iloc[:, 1:])
        
        #Agregate the responses and the thresholds then filter them and group them by experimental condition
        resp_db[recording_name] = responses
        thresh_db[recording_name] = thresholds
        grouped_by_treatment = group_by_treatment(resp_db, thresh_db)
        
        #Write the grouped data to a txt file for use at another time
        with open(output_dir + 'grouped_glia_responses_by_treatment' + str(date.today()) + '.pkl', 'wb') as of:
            pickle.dump(grouped_by_treatment, of)

#process_from_responses(input_dir = 'C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia_training\\analysis\\new glia activity\\data-peaks\\',
#                       output_dir = "C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia_training\\analysis\\new glia activity\\")
process_from_raw_traces(input_dir = 'E:/glia projects/plasticity/data/', 
                        output_dir = "E:/glia projects/plasticity/analysis/new glia activity/")



