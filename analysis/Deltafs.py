import os

import csv

import numpy as np

import statistics

def read_in_csv(file):

    with open(file) as fn:
        
        data = []

        reader = csv.reader(fn)

        for row in reader:

            data.append(row)
            
    return data





# Modify this function to return a percentile for EVERY cell in a column instead of a single percentile for the entire column

# We will end up with a n * m list where n = the number of cells in a column, and m = the number of columns in our data

# shape of percentiles should be equal to shape of original data (except the header and first column)

def calculate_percentiles(data, percentile=7.5):



    window_size = 150
    


    percentiles = [[] for i in range(len(data[0:]))]

    for cell, frames in enumerate(data): # One column at a time
    

        frames = list(map(float, frames[0:]))
 

        for frame, value in enumerate(frames[0:]): # One cell at a time, in a given column



            # Get the index for window/2 before and after column_index

            before = frame - (window_size / 2) if (frame > window_size / 2) else 0

            after = frame + (window_size / 2) if (frame + (window_size / 2)) < (len(frames) - 1) else len(frames) - 1

            #create the sliding window from the values of column with before and after indexes

            window = frames[int(before):int(after)]
            
            #Band filter the window to pull the 10th(exclusive) to 5th percentile(exclusive) values

            #btm_prcntl = list(filter(lambda a: a < np.percentile(window,percentile),window))
            
            #baseline_band = list(filter(lambda a: a > np.median(btm_prcntl),btm_prcntl))
            
            
            baseline = np.percentile(window, percentile)
            
            
            percentiles[cell].append(baseline)
            
    return percentiles



def do_math(data, percentiles):

    Fo_stdev = []
    deltaf = []

    count=0

    for cell, frames in enumerate(data):
        Fo = []
        count+=1
        new_column = []

        for frame, value in enumerate(frames):

            try:
                Fo.append(percentiles[cell][frame])
                new_column.append((float(value) - percentiles[cell][frame]) / float(percentiles[cell][frame]))

            except IndexError:
                Fo.append(percentiles[cell][frame-1])
                new_column.append((float(value) - percentiles[cell][frame-1]) / float(percentiles[cell][frame-1]))
        
        Fo_stdev.append([ f"trace_{count}",(1.5 * statistics.stdev(percentiles[cell])) / statistics.mean(percentiles[cell])]) #append the delta F/Fo of the 1.5*STDEV of Fo 
        deltaf.append(new_column)


    return Fo_stdev, deltaf



def output_new_csv(Fo, data, file):
      
    print("writing")
    if not os.path.exists(file[:file.rfind("L/")+2]+'deltaF/'):

        os.makedirs(file[:file.rfind("L/")+2]+'deltaF/')



    with open(file[:file.rfind("L/")+2] + 'deltaF/' + file[file.find("L/")+2:-4] + '_df.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in data:
            writer.writerows([row])
            
            
    with open(file[:file.rfind("L/")+2] + 'deltaF/' + file[file.find("L/")+2:-4] + '_Fo.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)
        
        writer.writerows(Fo)
    


def main():

    files = [i.path for i in os.scandir("C:/Users/David/Documents/Ruthazer lab/SOUL/") if i.path.endswith('.csv')]

    for file in files:
        
        data = read_in_csv(file)
        
        percentiles = calculate_percentiles(data)
        
        Fo_temp, deltaf = do_math(data, percentiles)

        output_new_csv(Fo_temp, deltaf, file)



if __name__ == "__main__":

    main()

