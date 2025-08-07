# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 11:16:53 2024

@author: BioCraze
"""
import numpy as np
import os, glob
import pandas as pd

input_dir = "C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia projects\\plasticity\\data\\"
files =[f.path[:f.path.rfind('\\')+1] for i in glob.glob(f'{input_dir}/**/',recursive=True) for f in os.scandir(i) if f.path.endswith('tectal_neuron_iscell.npy') and "capapplication" not in f.path]

celldb = pd.DataFrame(columns=("file",'iscell','neuron','glia'))
idx=0
for file in files:
    celllist = []
    celllist.append(file)
    celllist.append(np.load(file+"iscell.npy", allow_pickle=True).shape[0])
    celllist.append(np.load(file+"tectal_neuron_iscell.npy", allow_pickle=True).shape[0])
    celllist.append(np.load(file+"glia_iscell.npy", allow_pickle=True).shape[0])
    celldb.loc[idx]=celllist
    idx+=1
    
inconsistant = celldb[(celldb['iscell']!=celldb['neuron']) | (celldb['iscell']!=celldb['glia'])]
