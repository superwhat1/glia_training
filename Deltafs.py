import os

import csv

import numpy as np

import statistics

import pandas as pd

def read_in_neuropil_csv(file):
    raw = pd.read_csv(file)
    Ftrace = raw["Mean"]
    data = Ftrace.to_list()
    return data


def read_in_csv(file):
    with open(file) as fn:
        data = []
        reader = csv.reader(fn)
        for row in reader:
            float_row = [float(value) for value in row]
            data.append(float_row)
                
    return data



# Modify this function to return a percentile for EVERY cell in a column instead of a single percentile for the entire column

# We will end up with a n * m list where n = the number of cells in a column, and m = the number of columns in our data

# shape of percentiles should be equal to shape of original data (except the header and first column)

def calculate_neuropil_percentiles(data, percentile=7.5):



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


def neuropil_deltaF(data, percentiles):

    deltaf = []

    count=0

    for frames, value in enumerate(data):
        count+=1
        ratioed = []

        try:
            ratioed.append((float(value) - percentiles[frames]) / float(percentiles[frames]))

        except IndexError:
            ratioed.append((float(value) - percentiles[frames-1]) / float(percentiles[frames-1]))
         
        deltaf.append(ratioed)

    Fo_stdev =[ f"trace_{count}",(1.5 * statistics.stdev(percentiles)) / statistics.mean(percentiles)] #calculate the delta F/Fo of the 1.5*STDEV of Fo 
    return Fo_stdev, deltaf


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


def output_neuropil_csv(Fo, data, file):
      
    print("writing")
    if not os.path.exists(file[:file.rfind("/")+1]+'deltaF/'):

        os.makedirs(file[:file.rfind("/")+1]+'deltaF/')



    with open(file[:file.rfind("/")+1] + 'deltaF/' + file[file.rfind("y/")+2:-4] + '_df.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in data:
            writer.writerows([row])
            
            
    with open(file[:file.rfind("/")+1] + 'deltaF/' + file[file.rfind("y/")+2:-4] + '_Fo.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)
        
        writer.writerows(Fo)
    
    
def output_csv(Fo, data, file):
      
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
    


def main():

    files = [i.path for i in os.scandir("C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/neuropil activity/") if i.path.endswith('.csv')]
    
    tectal_area = "neuropil"
    
    if tectal_area == "neuropil":
        for file in files:
            
            data = read_in_neuropil_csv(file)
            
            percentiles = calculate_neuropil_percentiles(data)
            
            Fo_temp, deltaf = neuropil_deltaF(data, percentiles)
            
            output_neuropil_csv(Fo_temp, deltaf, file)
    
    else:
        for file in files:
            data = read_in_csv(file)
            
            percentiles = calculate_percentiles(data)
            
            Fo_temp, deltaf = deltaF(data, percentiles)
            
            output_csv(Fo_temp, deltaf, file)

#Run script
main()

