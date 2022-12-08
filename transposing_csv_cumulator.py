# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 19:35:50 2022

@author: BioCraze
"""
import csv
import os

import pandas as pd

#files = [i for i in os.listdir("C:/Users/BioCraze/Documents/Ruthazer lab/glial training/scripts/data/data-peaks/") if i.endswith('neuropil.csv')]
rootdir = 'C:\\Users\\BioCraze\\Documents\\Ruthazer lab\\glial training\\analysis\\neuropil activity'

files =[f.path for f in os.scandir(rootdir) if f.path.endswith('.csv')]

cumulative = []
for file in files:
    df = pd.read_csv(file)
    data = df["Mean"]
    data.transpose
    ready = pd.concat([pd.Series(file),data])
    cumulative.append(list(ready))

        
new_csv = open("C:/Users/BioCraze/Documents/Ruthazer lab/glial training/cumulated_neuropil.csv",'w',newline='')
write_to = csv.writer(new_csv)
for row in cumulative:
    write_to.writerows([row])

new_csv.close()