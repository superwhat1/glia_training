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

    for data_index, column in enumerate(data): # One column at a time
    

        column = list(map(float, column[0:]))
 

        for column_index, value in enumerate(column[0:]): # One cell at a time, in a given column



            # Get the index for window/2 before and after column_index

            before = column_index - (window_size / 2) if (column_index > window_size / 2) else 0

            after = column_index + (window_size / 2) if (column_index + (window_size / 2)) < (len(column) - 1) else len(column) - 1

            #create the sliding window from the values of column with before and after indexes

            window = column[int(before):int(after)]
            
            #Band filter the window to pull the 10th(exclusive) to 5th percentile(exclusive) values

            #btm_prcntl = list(filter(lambda a: a < np.percentile(window,percentile),window))
            
            #baseline_band = list(filter(lambda a: a > np.median(btm_prcntl),btm_prcntl))
            
            
            baseline = np.percentile(window, percentile)
            
            
            percentiles[data_index].append(baseline)
            
    return percentiles



def do_math(data, percentiles):

    Fo_stdev = []
    new_data_with_math = []

    count=0

    for data_index, column in enumerate(data):
        Fo = []
        count+=1
        new_column = []

        for column_index, value in enumerate(column):

            try:
                Fo.append(percentiles[data_index][column_index-1])
                new_column.append((float(value) - percentiles[data_index][column_index-1]) / float(percentiles[data_index][column_index-1]))

            except IndexError:
                Fo.append(percentiles[data_index][column_index-2])
                new_column.append((float(value) - percentiles[data_index][column_index-2]) / float(percentiles[data_index][column_index-2]))
        
        Fo_stdev.append([ f"trace_{count}",(1.5 * statistics.stdev(Fo)) / statistics.mean(Fo)]) #append the delta F/Fo of the 1.5*STDEV of Fo 
        new_data_with_math.append(new_column)


    return Fo_stdev, new_data_with_math



def output_new_csv(Fo, data, file):
      
    print("writing")
    if not os.path.exists(file+'deltaF/'):

        os.makedirs(file+'deltaF/')



    with open(file + 'deltaF/' + file[file.find("A"):-4] + '_df.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in data:
            writer.writerows([row])
            
            
    with open(file + 'deltaF/' + file[file.find("A"):-4] + '_Fo.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)
        
        writer.writerows(Fo)
    


def main():

    files = [i.path for i in os.scandir("C:/Users/BioCraze/Documents/Ruthazer lab/glial training/analysis/max proj roi activity/") if i.path.endswith('.csv')]

    for file in files:
        data = read_in_csv(file)
        
        percentiles = calculate_percentiles(data)
        
        Fo_temp, new_data = do_math(data, percentiles)

        output_new_csv(Fo_temp, new_data, file)



if __name__ == "__main__":

    main()

