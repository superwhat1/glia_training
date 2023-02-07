import os

import csv

from sklearn.metrics import auc


def read_in_csv(file):

    with open(file) as fn:

        data = []

        reader = csv.reader(fn)

        for row in reader:

            float_row = [float(numeric_string) for numeric_string in row]

            data.append(float_row)

    return data



def output_new_csv(area, peaks, times, file):
    
    for index, row in enumerate(area):

        row.insert(0, f"trace_{index+1}")
        
        
    for index, row in enumerate(peaks):

        row.insert(0, f"trace_{index+1}")
        

    for index, row in enumerate(times):

        row.insert(0, f"trace_{index+1}")
        

    if not os.path.exists('data-peaks/'):

        os.makedirs('data-peaks/')
        
        
    with open('data-peaks/' + file[file.find("A"):-7] + '_AREA.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in area:

            writer.writerows([row])


    with open('data-peaks/' + file[file.find("A"):-7] + '_PEAKS.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in peaks:

            writer.writerows([row])


    with open('data-peaks/' + file[file.find("A"):-7] + '_TIMES.csv', 'w', newline='') as fn:

        writer = csv.writer(fn)

        for row in times:

            writer.writerows([row])
            
            
            
def find_peaks_in_data(data):


    peaks = [[] for i in range(len(data))]

    times = [[] for i in range(len(data))]

    area = [[] for i in range(len(data))]

    windows = [540,
               1450,
               2360,
               3270,
               4180]


    for row_index, row in enumerate(data):

        for window in windows:
            if window - 250 < 0:
                left = 0
                right =window + 250
            elif window + 250 > 4500:
                left = window - 250
                right = 4500
            else:
                left = window - 250
                right = window + 250
            
            area[row_index].append(auc(list(range(left, right)), list(data[row_index][left:right])))

            peaks[row_index].append(max(data[row_index][left:right]))

            times[row_index].append(data[row_index].index(max(data[row_index][left:right])) + 1)
            

    return area, peaks, times     

       

def main(threshold_to_start, threshold_to_end):

#    files = [i for i in os.listdir() if i.endswith('.csv')]
#    for file in files:

    file = "C:/Users/BioCraze/Documents/Ruthazer lab/glial training/analysis/max proj roi activity/deltaF/A2_min80_27apr22_processed_is_cell_traces_df.csv"
    
    data = read_in_csv(file)

    area, peaks, times = find_peaks_in_data(data)

    output_new_csv(area, peaks, times, file)



if __name__ == "__main__":

    main()
