# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 19:35:50 2022

@author: BioCraze
"""
import csv
import glob
import pandas as pd

#files = [i for i in os.listdir("C:/Users/BioCraze/Documents/Ruthazer lab/glial training/scripts/data/data-peaks/") if i.endswith('neuropil.csv')]
rootdir = "C:\\Users\\BioCraze\\Documents\\Ruthazer lab\\glia projects\\plasticity\\analysis\\neuropil activity\\raw\\"

files =[f for f in glob.glob(f'{rootdir}\**', recursive=True) if f.endswith('.csv')]

cumulative = []
for file in files:
    df = pd.read_csv(file)
    data = df["Mean"]
    data_list = data.to_list()
    ready = pd.concat([pd.Series(file), data])
    cumulative.append(list(ready))
    #new_csv = open("E:\\glia training\\neuropil_notraining\\transposed\\" + file[file.find("A"):-4] + ".csv",'w',newline='')
    #write_to = csv.writer(new_csv)
    #write_to.writerow(data_list)
    #new_csv.close()
new_csv = open("C:\\Users\\BioCraze\\Documents\\Ruthazer lab\\glia projects\\plasticity\\analysis\\neuropil activity\\cumululated_neuropil.csv",'w',newline='')
write_to = csv.writer(new_csv)
for row in cumulative:
    write_to.writerows([row])
new_csv.close()