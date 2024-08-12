# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 14:30:21 2024

@author: David
"""

import numpy as np
from os import scandir
import glob

folders = [i.path for i in scandir('E:\glia projects\plasticity\data')]
output_dir = "C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia projects\\plasticity\\analysis\\new glia activity\\"

transforms = [["time_point","x","y","r"]] #where x and y are the translations to register a to b and r is the rotation.
for folder in folders:
    files =[f.path[:f.path.rfind('\\')+1] for i in glob.glob(f'{folder}/*/**/***/',recursive=True) for f in scandir(i) if f.path.endswith('ops.npy') and "capapplication" not in f.path]
        
    for file in files:
        raw_img = np.load(file,allow_pickle=True).item()["max_proj"]
        thr_img = raw_img[raw_img<raw_img.max()*0.2]=0
        
        #register images to min15 and update a list with 4 columns (time point compared to min15, x,y and rotaional shifts from the registration)
        transforms.append([file,x,y,r])