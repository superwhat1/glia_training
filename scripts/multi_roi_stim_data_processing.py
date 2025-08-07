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
    from scipy.signal import sosfiltfilt, butter
    
except ImportError as e:
    package = re.search("\'.+\'", str(e)).group()[1:-1]
    print(package + " not installed")
    subprocess.check_call([sys.executable, '-m', 'conda', 'install', package])


def is_cell_responsive(file, recording_name, output_dir):
    
    #Set folder path and then use it to open the iscell with shape: two columns, first contains 0 or 1 depending on if ROI is a cell and second the probability of that ROI of being a cell according to the suite2p classifier. 
    #And to open F.npy containing an array of fluorescent trace arrays
    is_cell = np.load(file + 'tectal_neuron_iscell.npy', allow_pickle=True)
    F = np.load(file + 'F.npy',allow_pickle=True)
    
    #create array of traces containing only the ROIs that are cells according to iscell
    F_fltrd = np.array([F[i] for i in range(len(is_cell)) if is_cell[i][0] == True])
    fltrd_roi_idxs = [i for i in range(len(is_cell)) if is_cell[i][0] == True]
    
    #write data to csv files
    np.savetxt(output_dir + recording_name +'_neuron_traces.csv', F_fltrd, delimiter=", ", fmt="% 1.4f")
    
    return F_fltrd

def lowpass(data):
    passed = np.zeros(data.shape)
    sos = butter(4,1, "lowpass", output='sos',fs='15',analog=False)
    for idx, trace in enumerate(data):
        passed[idx] = sosfiltfilt(sos,trace)
    return passed

def calculate_baselines(data, percentile=7.5):

    pre_window_size = 50
    post_window_size = 60
    baselines = [[] for i in range(len(data))]

    for cell, frames in enumerate(data): # One cell at a time
    
        for frame, value in enumerate(frames): # One frame at a time

            # Get the index for window/2 before and after frame
            before = frame - (pre_window_size) if (frame > pre_window_size) else 0
            after = frame + (post_window_size) if (frame + post_window_size) < (len(frames) - 1) else len(frames) - 1

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
    
    
def find_peaks_in_data(data, stims, recording_name):

    peaks = [[] for i in range(len(data))]
    times = [[] for i in range(len(data))]
    area = [[] for i in range(len(data))]
    responses = [[] for i in range(len(data))]
    thresholds = [[] for i in range(len(data))]
    first_stim = stims.at[recording_name, "stim 1"]
    
    for row_index, row in enumerate(data):
        
        for i in range(5):
            window = first_stim + i*915
            if window - 15 < 0:
                left = 0
                right =window + 165 - first_stim
            elif window + 150 > 4499:
                left = window - 165 + (4499 - window)
                right = 4499
            else:
                left = window - 15
                right = window + 150
            
            area[row_index].append(auc(list(range(left, right)), list(data[row_index][left:right])))
            
            if left == 0:
                rolling_trgt = 15 - first_stim
                responses[row_index].append(np.roll(np.array(data[row_index][left:right]), rolling_trgt))
            
            elif right == 4499:
                rolling_trgt = 4499 - window - 15
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
            np.save(fn, np.array(responses), allow_pickle=True)
                

def thresh_filter_responses(resp_db, thresh_db): #response arrays should have a shape of cells*stims*time
            
    #remove cell responses where the max is less than the threshold        
    for recording, cells in resp_db.items():
        for cell, responses in enumerate(cells):
            for response, trace in enumerate(responses):
                if trace.max() < thresh_db[recording][cell][response]:
                    resp_db[recording][cell][response] *= np.nan
                    
    fltrd_resp_db = resp_db
    
    return fltrd_resp_db


def mean_stack(stacked):
    
    meaned = np.nanmean(stacked, axis=0, dtype = np.float64)
    
    return  meaned


def group_by_treatment(resp_db, thresh_db):
    treatments = {"cap_and_train": ["21sep22", "22jun22", "26may22", "13aug22"], "cap_notrain": ["01apr23", "02feb23", "30mar23"],
    "nocap_train": ["07sep22", "15jun22", "20jul22", "27apr22"], "nocap_notrain":["19may23", "23aug23"]}
    
    fltrd_resp_db = thresh_filter_responses(resp_db, thresh_db)
    
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
            if treatment =="nocap_notrain":
                to_stack = []
                for animal, times in animals.items():    
                    for time, cells in times.items():
                        for traces in cells:
                            if i in time:
                                to_stack.append(traces)
                                
                stacked = np.stack(to_stack)
                meaned = mean_stack(stacked)
                #meaned_again = mean_stack(meaned)
                
                if i == 'min15':
                    base_min[treatment] = meaned.min()
                    base_max[treatment] = meaned.max()
                normalized = (meaned - base_min[treatment])/ (base_max[treatment] - base_min[treatment])
                ready_to_plot[i][treatment]=normalized
                
                plt.title(i)
                plt.plot(normalized[:],"cyan", label=treatment)
                plt.ylim(-0.1,1.5)
                
            elif treatment == "cap_notrain":
                to_stack = []
                for animal, times in animals.items():    
                    for time, cells in times.items():
                        for traces in cells:
                            if i in time:
                                to_stack.append(traces)
                                
                stacked = np.stack(to_stack)
                meaned = mean_stack(stacked)
                #meaned_again = mean_stack(meaned)
                
                if i == 'min15':
                    base_min[treatment] = meaned.min()
                    base_max[treatment] = meaned.max()
                normalized = (meaned - base_min[treatment])/ (base_max[treatment] - base_min[treatment])
                ready_to_plot[i][treatment]=normalized
                
                plt.title(i)
                plt.plot(normalized[:],"magenta", label=treatment)
                plt.ylim(-0.1,1.5)
        #plt.legend(loc="best")
        
        plt.show()


def make_plots(grouped_by_treatment):#Function for ploting aligned and averaged responses

    #Prepare variables needed to plot average traces
    time_list=['min15', 'min45', 'min60', 'min80', 'min100']
    
    plot_averaged(grouped_by_treatment, time_list)
    
def process_from_raw_traces(input_dir, output_dir, stims_timings):

    files =[f.path[:f.path.rfind('\\')+1] for i in glob.glob(f'{input_dir}/*/**/***/',recursive=True) for f in os.scandir(i) if f.path.endswith('tectal_neuron_iscell.npy') and "capapplication" not in f.path]
    
    stims = pd.read_csv(stims_timings)
    stims = stims.set_index("file")
    
    resp_db = {}
    thresh_db = {}
    for file in files:
        #find the name of the recording that starts with an A and ends with a date in the format of ddmmmyy with dd and yy as ints and mmm as string
        recording_name = re.search("A\d{1}_min.+\d{2}\D{3}\d{2}", file).group()
        print(recording_name)
        #perform deltaF operation
        #try:
        print("Normalizing " + file)
        raw_trace = is_cell_responsive(file, recording_name, output_dir)
        baselines = calculate_baselines(raw_trace)
        deltaf = deltaF(raw_trace, baselines)
        #except Exception:
        #    print(file + " contains no traces, so it was skipped!")
        #    pass
        
        #perform response metric operations 
        try:
            print("Extracting response metrics of " + file)
            area, peaks, times, thresholds, responses = find_peaks_in_data(deltaf, stims, recording_name)
        
            #Save the data to files
            output_csv(baselines, deltaf, area, peaks, times, thresholds, responses, output_dir, recording_name)
            
            #Agregate the responses and the thresholds then filter them and group them by experimental condition
            resp_db[recording_name] = np.array(responses)
            thresh_db[recording_name] = thresholds
            
        except KeyError:
            print(file + " contains no responses, so it was skipped!")
            pass
        
    grouped_by_treatment = group_by_treatment(resp_db, thresh_db)
    
    #Write the grouped data to a txt file for use at another time
    with open(output_dir + 'grouped_neurons' + '_responses_by_treatment' + str(date.today()) + '.pkl', 'wb') as of:
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
        with open(output_dir + 'grouped_neurons' + '_responses_by_treatment' + str(date.today()) + '.pkl', 'wb') as of:
            pickle.dump(grouped_by_treatment, of)

#process_from_responses(input_dir = 'C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia projects\\plasticity\\analysis\\max proj roi activity\\data-peaks\\', output_dir = "C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia projects\\plasticity\\analysis\\")
process_from_raw_traces(input_dir = "E:/glia projects/plasticity/data", output_dir = "C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia projects\\plasticity\\analysis\\new max proj roi activity\\", stims_timings = "C:/Users/BioCraze/Documents/Ruthazer lab/glia projects/plasticity/summaries/stim_timings.csv")

#access roi location to match based on position
#for roi in iscell:
    
#    stat[roi]["med"]

    