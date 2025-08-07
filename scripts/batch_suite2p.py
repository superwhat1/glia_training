# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 16:29:30 2022

@author: BioCraze
"""

import os
import numpy as np
import glob
from suite2p import run_s2p

ops = np.load('E:/glia projects/anatomical_mean_roi_detection_ops.npy', allow_pickle=True).item()
rootdir = 'E:/glia projects/plasticity/data/training_A2_27apr22/'

files = [f.path[:f.path.rfind('\\')+1] for i in glob.glob(f'{rootdir}/**/',recursive=True) for f in os.scandir(i) if f.path.endswith('.tif') and 'min' in  f.path]
for file in files:
    db = {'data_path':[file]}
    opsEND = run_s2p(ops=ops, db=db)