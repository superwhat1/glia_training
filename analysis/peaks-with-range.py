import os

import csv

from sklearn.metrics import auc

import numpy as np

import pandas as pd

def read_in_csv(file):

    with open(file) as fn:

        data = []

        reader = csv.reader(fn)

        for row in reader:

            float_row = [float(numeric_string) for numeric_string in row]

            data.append(float_row)

    return data

def roll_and_average(data): 
    
    data = np.array(data)   
            
    meaned = np.mean(data,axis=0).tolist()
    return  meaned

def output_new_csv(area, peaks, times, threshold, responses, meaned, file):
    if len(area) > 0: 
        for index, row in enumerate(area):
    
            row.insert(0, f"trace_{index+1}")
            
            
        for index, row in enumerate(peaks):
    
            row.insert(0, f"trace_{index+1}")
            
    
        for index, row in enumerate(times):
    
            row.insert(0, f"trace_{index+1}")
            
        for index, row in enumerate(threshold):
    
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
    
            for row in threshold:
    
                writer.writerows([row])
                
        with open('data-peaks/' + file[file.find("A"):-4] + '_RESPONSES.npy', 'wb') as fn:
            
            np.save(fn, np.array(responses), allow_pickle=True)
                
        with open('data-peaks/' + file[file.find("A"):-4] + "_ROLLEDandMEANED.csv", 'w', newline = '') as fn:
            
            writer = csv.writer(fn)
            for row in meaned:
                writer.writerow(row)


def find_peaks_in_data(data, stims, animal):


    peaks = [[] for i in range(len(data))]

    times = [[] for i in range(len(data))]

    area = [[] for i in range(len(data))]
    
    threshold = [[] for i in range(len(data))]

    responses = [[] for i in range(len(data))]
    
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
                threshold[row_index].append(max(data[row_index][left-100:left]))
            except:
                threshold[row_index].append(max(data[row_index][right:right+100]))

    return area, peaks, times, threshold, responses


def main():
    stims_timings = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/summaries/stim_timings.csv"
    stims = pd.read_csv(stims_timings)
    stims = stims.set_index("file")
    
    exp_fold = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/max proj roi activity/deltaF/"
    
    files = [i.path for i in os.scandir(exp_fold) if i.path.endswith('df.csv')]

    for file in files:
        animal = file[file.rfind("A"):file.rfind("neu")-1]
        data = read_in_csv(file)
        try:
            area, peaks, times, threshold, responses = find_peaks_in_data(data, stims, animal)
    
            meaned = roll_and_average(responses)
        
            output_new_csv(area, peaks, times, threshold, responses, meaned, file)
        except KeyError:
            print(file)
            pass


if __name__ == "__main__":

    main()
