import csv, os
from sklearn.metrics import auc
import pandas as pd
import numpy as np


def read_in_csv(file):

    df = pd.read_csv(file, header = None)
    data = df.iloc[:,0].tolist()
        
    return data


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

       
def main():

    stim_times = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/summaries/stim_timings.csv"
    stims = pd.read_csv(stim_times)
    stims = stims.set_index("file")
    
    exp_fold = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/neuropil activity/deltaF/"
    
    files = [i.path for i in os.scandir(exp_fold) if i.path.endswith('df.csv')]

    for file in files:
        animal = file[file.rfind("A"):file.rfind("neu")-1]
        
        data = read_in_csv(file)
        try:
            area, peaks, times, threshold, responses = find_peaks_in_data(data, stims, animal)
        
            output_new_csv(area, peaks, times, threshold, responses, file)
      
        except KeyError:
            pass

#Run script

main()
