# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 14:30:21 2024

@author: David
"""

import numpy as np
from os import scandir, path, makedirs
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from tifffile import imwrite
from re import search
folders = [i.path for i in scandir('E:\glia projects\plasticity\data')]
output_dir = "C:\\Users\\Biocraze\\Documents\\Ruthazer lab\\glia projects\\plasticity\\analysis\\new glia activity\\"

transforms = [["time_point","x","y","r"]] #where x and y are the translations to register a to b and r is the rotation.

normalization = Normalize()
for folder in folders:
    files =[f.path[:f.path.rfind('\\')+1] for i in glob.glob(f'{folder}/*/**/***/',recursive=True) for f in scandir(i) if f.path.endswith('ops.npy') and "capapplication" not in f.path]
        
    for file in files:
        raw_img = np.load(file + "ops.npy",allow_pickle=True).item()["meanImg"]
        normed_img = normalization(raw_img)
        thr_img_1 = np.copy(normed_img)
        thr_img_2 = np.copy(normed_img)
        thr_img_1[thr_img_1<0.001]=0
        thr_img_2[thr_img_2<0.005]=0
        
        fig,axes =plt.subplots(nrows=1, ncols=3, figsize = (30,10), layout="constrained")
        axes[0].imshow(normed_img, cmap='gray',vmin=0,vmax=0.15)
        axes[1].imshow(thr_img_1, cmap='gray',vmin=0,vmax=0.15)
        axes[2].imshow(thr_img_2, cmap='gray',vmin=0,vmax=0.15)
        
        if not path.exists(output_dir + 'reg_mean_projs/'):
            makedirs(output_dir + 'reg_mean_proj/')
        imwrite(folder+"reg_mean_projs/"+search("A\d{1}_min.+\d{2}\D{3}\d{2}",file).group())
        #register images to min15 and update a list with 4 columns (time point compared to min15, x,y and rotaional shifts from the registration)
        transforms.append([file,x,y,r]) 