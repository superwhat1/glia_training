import os

import csv

from sklearn.metrics import auc

import pandas as pd

import numpy as np

def read_in_csv(file):

    df = pd.read_csv(file, header = None)
    data = df.iloc[:,0].tolist()
            
    return data


def output_new_csv(area, peaks, times, threshold, responses, meaned, file):
        

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
            
    with open('data-peaks/' + file[file.find("A"):-4] + "_ROLLEDandMEANED.csv", 'w', newline = '') as fn:
        
        writer = csv.writer(fn)
        
        writer.writerow(meaned)
            
            
def roll_and_average(data): 
    
    data = np.array(data)
    
    target_idx = int(np.floor(data.shape[-1]/2))

    rolled = np.zeros(np.shape(data))

    for session, responses in enumerate(data):
        for response, time in enumerate(responses):
            rolled[session,response] = np.roll(time, target_idx - np.argmax(time))   
            
    meaned = np.mean(rolled,axis=0).tolist()
    return rolled, meaned

def find_peaks_in_data(data):


    peaks = []

    times = []

    area = []
    
    threshold = []

    responses = []
    
    a=0
    
    first_stim = 120
    
    a+=1
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

#    files = [i for i in os.listdir() if i.endswith('.csv')]
#    for file in files:

    file = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/neuropil activity/deltaF/A2_min100_27apr22_df.csv"
    
    data = read_in_csv(file)

    area, peaks, times, threshold, responses = find_peaks_in_data(data)

    rolled, meaned = roll_and_average(responses)
    
    output_new_csv(area, peaks, times, threshold, responses, meaned, file)



if __name__ == "__main__":

    main()
