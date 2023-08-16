# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 15:15:09 2023

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
    import matplotlib as plt
    import pickle
    from datetime import date
    
except ImportError as e:
    package = re.search("\'.+\'", str(e)).group()[1:-1]
    print(package + " not installed")
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])

def read_in_neuropil_csv(file):
    
    raw = pd.read_csv(file)
    Ftrace = raw["Mean"]
    data = Ftrace.to_list()
    
    return data


def calculate_neuropil_baselines(data, percentile=7.5):

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


def neuropil_deltaF(data, baselines):

    deltaf = []
    count=0

    for frames, value in enumerate(data):
        count+=1
        ratioed = []

        try:
            ratioed.append((float(value) - baselines[frames]) / float(baselines[frames]))

        except IndexError:
            ratioed.append((float(value) - baselines[frames-1]) / float(baselines[frames-1]))
         
        deltaf.append(ratioed)

    Fo_stdev =[ f"trace_{count}",(1.5 * statistics.stdev(baselines)) / statistics.mean(baselines)] #calculate the delta F/Fo of the 1.5*STDEV of Fo 
    return Fo_stdev, deltaf


def find_peaks_in_data(data, stim_times, animal):

    peaks = []
    times = []
    area = []
    threshold = []
    responses = []
    
    first_stim = stim_times.at[animal, "first stim"]

    for i in range(5):
        window = first_stim + i*915
        if window - 250 < 0:
            left = 0
            right =window + 500 - first_stim
        elif window + 250 > 4500:
            left = window - 500 + (4500 - window)
            right = 4500
        else:
            left = window - 250
            right = window + 250
            
        area.append(auc(list(range(left, right)), list(data[left:right])))
            
        responses.append(data[left:right])
            
        peaks.append(max(data[left:right]))

        times.append(data.index(max(data[left:right])) + 1)
        
        try:
            threshold.append(max(data[left-100:left]))
        except:
            threshold.append(max(data[right:right+100]))
                
    return area, peaks, times, threshold, responses


def output_neuropil_csv(Fo, deltaf, area, peaks, times, thresholds, responses, output_dir, recording_name):
      
    print("writing")
    if not os.path.exists(output_dir + 'deltaF/'):
        os.makedirs(output_dir + '/deltaF/')

    with open(output_dir + 'deltaF/' + recording_name + '_df.csv', 'w', newline='') as fn:
        writer = csv.writer(fn)
        for row in deltaf:
            writer.writerows([row])
            
            
    with open(output_dir + 'deltaF/' + recording_name + '_Fo.csv', 'w', newline='') as fn:
        writer = csv.writer(fn)
        writer.writerows(Fo)
        
    if len(area) > 0:
    
        if not os.path.exists(output_dir + 'data-peaks/'):
            os.makedirs(output_dir + 'data-peaks/')
            
        with open(output_dir + 'data-peaks/' + recording_name + '_AREA.csv', 'w', newline='') as fn:
            writer = csv.writer(fn)
            writer.writerow(area)
    
        with open(output_dir + 'data-peaks/' + recording_name + '_PEAKS.csv', 'w', newline='') as fn:
            writer = csv.writer(fn)
            writer.writerow(peaks)
    
    
        with open(output_dir + 'data-peaks/' + recording_name + '_TIMES.csv', 'w', newline='') as fn:
            writer = csv.writer(fn)
            writer.writerow(times)
                
        with open(output_dir + 'data-peaks/' + recording_name + '_THRESHOLD.csv', 'w', newline='') as fn:
            writer = csv.writer(fn)
            writer.writerow(thresholds)
        
        with open(output_dir + 'data-peaks/' + recording_name + '_RESPONSES.npy', 'wb') as fn:
            np.save(fn, np.array(responses), allow_pickle=True)
        
            
def thresh_filter_responses(resp_db, thresh_db): #response arrays should have a shape of cells*stims*time
            
    #remove cell responses where the max is less than the threshold        
    for k, v in resp_db.items():
        for i in range(len(v)):
            if v[i].max() < thresh_db[k][i]:
                resp_db[k][i] *= np.nan
                
    fltrd_resp_db = resp_db
    
    return fltrd_resp_db


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

def group_by_treatment(resp_db, thresh_db):
    treatments = {"cap_and_train": ["21sep22", "22jun22", "26may22", "13aug22"], "cap_notrain": ["01apr23", "02feb23", "30mar23"],
    "nocap_train": ["07sep22", "15jun22", "20jul22", "27apr22"], "nocap_notrain":["19may23"]}
    
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
        
        
def make_plots(grouped_by_treatment):
    
    #Prepare variables needed to plot average traces    
    time_list=['min15', 'min45', 'min60', 'min80', 'min100']
    
    plot_averaged(grouped_by_treatment, time_list)
 
    
def process_from_raw_traces():
    #Run data processing
    input_dir = 'C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia_training\\analysis\\neuropil activity\\raw'
    files =[i.path for i in os.scandir(input_dir) if i.path.endswith('.csv')]
    output_dir = "C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia_training\\analysis\\neuropil activity\\"
    
    stims_timings = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/summaries/stim_timings.csv"
    stims = pd.read_csv(stims_timings)
    stims = stims.set_index("file")
    
    resp_db = {}
    thresh_db = {}
    
    for file in files:
        #perform deltaF operation
        raw_trace = read_in_neuropil_csv(file)
        baselines = calculate_neuropil_baselines(raw_trace)
        Fo, deltaf = neuropil_deltaF(raw_trace, baselines)
        
        #perform response metric operations
        recording_name = re.search("A.+\d{2}\D{3}\d{2}", file).group() #find the name of the recording that starts with an A and ends with a date in the format of ddmmmyy with dd and yy as ints and mmm as string
        
        try:
            area, peaks, times, thresholds, responses = find_peaks_in_data(deltaf, stims, recording_name)
            
        except KeyError:
            print(file + " contains no traces, so it was skipped!")
            pass
        
        #Save the data to files
        output_neuropil_csv(Fo, deltaf, area, peaks, times, thresholds, responses, output_dir, recording_name)
        
        #Agregate the responses and the thresholds then filter them and group them by experimental condition
        resp_db[recording_name] = responses
        thresh_db[recording_name] = thresholds
    
    grouped_by_treatment = group_by_treatment(resp_db, thresh_db)
        
    #Write the grouped data to a txt file for use at another time
    with open(output_dir + 'grouped_neuropil_responses_by_treatment' + str(date.today()) + '.pkl', 'wb') as of:
        pickle.dump(grouped_by_treatment, of)
    

def process_from_responses():
    
    input_dir = 'C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia_training\\analysis\\neuropil activity\\data-peaks'
    response_files =[f.path  for f in os.scandir(input_dir) if f.path.endswith('RESPONSES.npy')]
    output_dir = "C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia_training\\analysis\\neuropil activity\\"
    
    resp_db = {}
    thresh_db = {}
    
    for file in response_files:
        #find the name of the recording that starts with an A and ends with a date in the format of ddmmmyy with dd and yy as ints and mmm as string
        recording_name = re.search("A.+\d{2}\D{3}\d{2}", file).group()
        responses = np.load(file, allow_pickle=True)
        thresholds = pd.read_csv(file[:file.rfind('_')] + '_THRESHOLD.csv', header=None).iloc[0,:].to_list()
        
        #Agregate the responses and the thresholds then filter them and group them by experimental condition
        resp_db[recording_name] = responses
        thresh_db[recording_name] = thresholds
        
    grouped_by_treatment = group_by_treatment(resp_db, thresh_db)
    
        #Write the grouped data to a txt file for use at another time
    with open(output_dir + 'grouped_neuropil_responses_by_treatment' + str(date.today()) + '.pkl', 'wb') as of:
        pickle.dump(grouped_by_treatment, of)

process_from_responses()