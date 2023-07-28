# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 15:15:43 2023

@author: BioCraze
"""
import numpy as np
import pandas as pd
import os, glob, statistics, csv, re
from sklearn.metrics import auc
import matplotlib.pyplot as plt

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
    
    
def find_peaks_in_data(data, stims, recording_name):


    peaks = [[] for i in range(len(data))]

    times = [[] for i in range(len(data))]

    area = [[] for i in range(len(data))]
    
    responses = [[] for i in range(len(data))]
    
    thresholds = [[] for i in range(len(data))]
    
    a=0
    
    first_stim = stims.at[recording_name, "first stim"]


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
    
    
    
def output_csv(Fo, deltaf, area, peaks, times, thresholds, responses, output_dir, recording_name):
      
    print("writing")
    if not os.path.exists(output_dir + '/deltaF/'):

        os.makedirs(output_dir + '/deltaF/')



    with open(output_dir + '/deltaF/' + recording_name + '_df.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in deltaf:
            writer.writerows([row])
            
            
    with open(output_dir + '/deltaF/' + recording_name + '_Fo.csv', 'w', newline='') as fn:

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
            
            
        if not os.path.exists(output_dir + '/data-peaks/'):
    
            os.makedirs(output_dir + '/data-peaks/')
            
            
        with open(output_dir + '/data-peaks/' + recording_name + '_AREA.csv', 'w', newline='') as fn:
    
            writer = csv.writer(fn)
    
            for row in area:
    
                writer.writerows([row])
    
    
        with open(output_dir + '/data-peaks/' + recording_name + '_PEAKS.csv', 'w', newline='') as fn:
    
            writer = csv.writer(fn)
    
            for row in peaks:
    
                writer.writerows([row])
    
    
        with open(output_dir + '/data-peaks/' + recording_name + '_TIMES.csv', 'w', newline='') as fn:
    
            writer = csv.writer(fn)
    
            for row in times:
    
                writer.writerows([row])
                
        with open(output_dir + '/data-peaks/' + recording_name + '_THRESHOLD.csv', 'w', newline='') as fn:
            
            writer = csv.writer(fn)
    
            for row in thresholds:
    
                writer.writerows([row])
                
        with open(output_dir + '/data-peaks/' + recording_name + '_RESPONSES.npy', 'wb') as fn:
            
            np.save(fn, np.array(responses), allow_pickle=True)
                

def thresh_filter_responses(resp_db, thresh_db): #response arrays should have a shape of cells*stims*time
            
    #remove cell responses where the max is less than the threshold        
    for k, v in resp_db.items():
        for w, x in enumerate(v):
            for y, z in enumerate(x):
                if z.max() < thresh_db[k][w][y]:
                    resp_db[k][w][y] *= np.nan
                    
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

def group_by_treatment(treatments, fltrd_resp_db):
    grouped_by_treatment ={}
    for key, values in treatments.items():
        grouped_by_treatment[key] ={}
        for value in values:
            kl=[]
            vl=[]
            for k, v in fltrd_resp_db.items():
                if value in k:
                    kl.append(k)
                    vl.append(v)
            kv = dict(zip(kl,vl))
            grouped_by_treatment[key][value] = kv
    return grouped_by_treatment  


#Set parameters and run code
def plot_averaged(grouped_by_treatment, time_list):
    
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



input_dir = 'C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia_training\\data\\training_19may23'
files =[f.path[:f.path.rfind('\\')+1] for i in glob.glob(f'{input_dir}/*/**/***/',recursive=True) for f in os.scandir(i) if f.path.endswith('iscell.npy')]
output_dir = "C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia_training\\data"

stims_timings = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/summaries/stim_timings.csv"
stims = pd.read_csv(stims_timings)
stims = stims.set_index("file")

resp_db = {}
thresh_db = {}

for file in files:
    #perform deltaF operation
    raw_trace = is_cell_responsive(file)
    percentiles = calculate_percentiles(raw_trace)
    Fo, deltaf = deltaF(raw_trace, percentiles)
        
    #perform response metric operations
    recording_name = re.search("A.+\d{2}\D{3}\d{2}", file).group() #find the name of the recording that starts with an A and ends with a date in the format of ddmmmyy with dd and yy as ints and mmm as string
    
    try:
        area, peaks, times, thresholds, responses = find_peaks_in_data(deltaf, stims, recording_name)
        
    except KeyError:
        print(file)
        pass
    
    output_csv(Fo, deltaf, area, peaks, times, thresholds, responses, output_dir, recording_name)
    
    resp_db[recording_name] = responses
    thresh_db[recording_name] = thresholds


def make_plots(resp_db, thresh_db):
    #Prepare variables needed to plot average traces
    treatments = {"cap_and_train": ["21sep22", "22jun22", "26may22", "13aug22"], "cap_notrain": ["01apr23", "02feb23", "30mar23"],
    "nocap_train": ["07sep22", "15jun22", "20jul22", "27apr22"]}
    time_list=['min15', 'min45', 'min60', 'min80', 'min100']
    
    fltrd_resp_db = thresh_filter_responses(resp_db, thresh_db)
    grouped_by_treatment = group_by_treatment(treatments, fltrd_resp_db)
    
    plot_averaged(grouped_by_treatment, time_list)
    
    
#Run data processing
