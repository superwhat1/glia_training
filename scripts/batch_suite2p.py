# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 16:29:30 2022

@author: BioCraze
"""

import os
import numpy as np
import glob
from suite2p import run_s2p

ops = np.load('C:/Users/BioCraze/Documents/Ruthazer lab/glia projects/anatomical_max_roi_detection_ops.npy', allow_pickle=True).item()
rootdir = 'C:\\Users\\BioCraze\\Documents\\Ruthazer lab\\glia projects\\plasticity\\data\\'

batches = [f.path for f in os.scandir(rootdir) if 'training' in  f.path]

for b in batches:
    files = [f.path[:f.path.rfind('\\')+1] for i in glob.glob(f'{b}/**/',recursive=True) for f in os.scandir(i) if f.path.endswith('.tif') and 'min' in  f.path]
    files = files[1:]+files[:1]
    db = {'data_path':[files]}
    opsEND = run_s2p(ops=ops, db=db)