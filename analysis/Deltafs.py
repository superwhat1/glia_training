import os

import csv

import numpy as np

import statistics

import pandas as pd

def read_in_csv(file):
        
    df = pd.read_csv(file)
    data = df["Mean"]
    data_list = data.to_list()
            
    return data_list





# Modify this function to return a percentile for EVERY cell in a column instead of a single percentile for the entire column

# We will end up with a n * m list where n = the number of cells in a column, and m = the number of columns in our data

# shape of percentiles should be equal to shape of original data (except the header and first column)

def calculate_percentiles(data, percentile=7.5):



    window_size = 150
    
    percentiles = []
    baselines = []

    for frames, value in enumerate(data): # One column at a time
    
        '''
        frames = list(map(float, frames[0:]))
 

        for frame, value in enumerate(frames[0:]): # One cell at a time, in a given column
        '''


            # Get the index for window/2 before and after column_index

        before = frames - (window_size / 2) if (frames > window_size / 2) else 0

        after = frames + (window_size / 2) if (frames + (window_size / 2)) < (len(data) - 1) else len(data) - 1

            #create the sliding window from the values of column with before and after indexes

        window = data[int(before):int(after)]
            
            #Band filter the window to pull the 10th(exclusive) to 5th percentile(exclusive) values

            #btm_prcntl = list(filter(lambda a: a < np.percentile(window,percentile),window))
            
            #baseline_band = list(filter(lambda a: a > np.median(btm_prcntl),btm_prcntl))
            
            
        baseline = np.percentile(window, percentile)
        baselines.append(baseline)
            
    percentiles = baselines

    return percentiles



def deltaF(data, percentiles):

    Fo_stdev = []
    temp = []
    deltaf = []

    count=0

    for frames, value in enumerate(data):
        count+=1
        new_column = []

        try:
            new_column.append((float(value) - percentiles[frames]) / float(percentiles[frames]))

        except IndexError:
            new_column.append((float(value) - percentiles[frames-1]) / float(percentiles[frames-1]))
        
        temp.append([ f"trace_{count}",(1.5 * statistics.stdev(percentiles)) / statistics.mean(percentiles)]) #append the delta F/Fo of the 1.5*STDEV of Fo 
        deltaf.append(new_column)

    Fo_stdev = temp
    return Fo_stdev, deltaf



def output_new_csv(Fo, data, file):
      
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
    


def main():

    files = [i.path for i in os.scandir("C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/neuropil activity/") if i.path.endswith('.csv')]

    for file in files:
        
        data = read_in_csv(file)
        
        percentiles = calculate_percentiles(data)
        
        Fo_temp, deltaf = deltaF(data, percentiles)

        output_new_csv(Fo_temp, deltaf, file)



if __name__ == "__main__":

    main()

