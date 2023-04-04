# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 16:29:30 2022

@author: BioCraze
"""
import os
import numpy as np
import glob
from suite2p import run_s2p

ops = np.load('C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/anatomical_max_roi_detection_ops.npy', allow_pickle=True).item()
rootdir = 'C:\\Users\\BioCraze\\Documents\\Ruthazer lab\\glia_training\\data\\training_02feb23\\'

files =[f.path[:f.path.rfind('\\')+1] for i in glob.glob(f'{rootdir}/*/**/',recursive=True) for f in os.scandir(i) if f.path.endswith('.tif')]

for f in files:
    db = {'data_path':[f]}
    opsEND = run_s2p(ops=ops, db=db)