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



def output_new_csv(area, peaks, times, threshold, file):
    
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
            
            
def find_peaks_in_data(data):


    peaks = [[] for i in range(len(data))]

    times = [[] for i in range(len(data))]

    area = [[] for i in range(len(data))]
    
    threshold = [[] for i in range(len(data))]

    first_stim = 380


    for row_index, row in enumerate(data):
        for i in range(5):
            window = first_stim + i*915
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
            
            threshold[row_index].append(max(data[row_index][right:right+150]))

    return area, peaks, times, threshold

       

def main():

#    files = [i for i in os.listdir() if i.endswith('.csv')]
#    for file in files:

    file = "E:/glia training/neuron_notraining/deltaF/A1_cap_notraining_min100_30mar23_neuron_traces_df.csv"
    
    data = read_in_csv(file)

    area, peaks, times, threshold = find_peaks_in_data(data)

    output_new_csv(area, peaks, times, threshold, file)



if __name__ == "__main__":

    main()
