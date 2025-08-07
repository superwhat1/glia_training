
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 19:35:50 2022

@author: BioCraze
"""
import csv
import os




files = [i.path for i in os.scandir("C:/Users/BioCraze/Documents/Ruthazer lab/glia projects/plasticity/analysis/new max proj roi activity/data-peaks") if i.path.endswith('THRESHOLDS.csv')]



cumulative = []
for file in files:
    f = open(file)
    read = csv.reader(f)
    for row in read:
        row.insert(0, file)
        cumulative.append(row)

        



new_csv = open("C:/Users/BioCraze/Documents/Ruthazer lab/glia projects/plasticity/analysis/new max proj roi activity/cumulated_thresholds.csv",'w',newline='')

write_to = csv.writer(new_csv)
for row in cumulative:
    write_to.writerows([row])

new_csv.close()
