# -*- coding: utf-8 -*-
"""
Created on Fri May 23 15:38:26 2025

@author: David
"""

from skimage.io import imread
from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage.filters import median, threshold_triangle
from skimage.morphology import octagon, binary_erosion,binary_dilation
from os import listdir
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import ceil
import matplotlib as mlp

mlp.rcParams['axes.titlesize']=20
folder='E:/transfection/staining/cropped_cldu_staining//'
files = [i for i in listdir(folder) if "tif" in i]
files.sort()
data = {"animal":[],"percent green":[], "slices":[]}
octagon = octagon(5,2,decomposition='sequence')

patch_kw = dict(
patch_size=4,       # 3x3 patches
patch_distance=8,  # 13x13 search area
)                   #Bigger is better but slower


for idx in range(0,len(files),2):
    print(files[idx])
    red_chan = imread(folder+files[idx])
    green_chan = imread(folder+files[idx+1])
    slices = red_chan.shape[0]
    
    animal_min = np.min([red_chan,green_chan],axis=(0,2,3)).reshape((slices,1,1))
    animal_max = np.max([red_chan,green_chan],axis=(0,2,3)).reshape((slices,1,1))
    n_red = ((red_chan-animal_min)/(animal_max-animal_min))*255
    n_green = ((green_chan-animal_min)/(animal_max-animal_min))*255

    mf_red = median(n_red)
    mf_green = median(n_green)
    
    sigma_red = np.array(estimate_sigma(mf_red,channel_axis=0))
    sigma_green = np.array(estimate_sigma(n_green,channel_axis=0))
    dn_red=np.copy(mf_red)
    dn_green=np.copy(mf_green)
    
    for i in range(len(sigma_red)):    
        dn_red[i]= denoise_nl_means(dn_red[i,:,:], h=0.8 * sigma_red[i], sigma=sigma_red[i], fast_mode=False, **patch_kw)
        dn_green[i]=denoise_nl_means(dn_green[i], h=0.8 * sigma_green[i], sigma=sigma_green[i], fast_mode=False, **patch_kw)
    
    gminusr = dn_green-dn_red
    gminusr[gminusr<0]=0
    
    mask_red=np.zeros_like(dn_red)
    mask_green = np.zeros_like(gminusr)
    
    for i in range(slices):
        mask_red[i][dn_red[i]>threshold_triangle(dn_red[i],nbins=256)]=1
        mask_green[i][gminusr[i]>threshold_triangle(gminusr[i],nbins=256)]=1
        
        mask_red = binary_erosion(mask_red, mode='min')
        mask_red = binary_dilation(mask_red,footprint=octagon,mode='min')
        mask_red = binary_erosion(mask_red,footprint=octagon,mode='min')

#calculate the green over red overlaping area percentage
    percent_green = np.mean(np.sum(np.matmul(mask_green,mask_red),axis=(1,2))/np.sum(mask_red,axis=(1,2)))

    data["animal"].append(files[idx])
    data["percent green"].append(percent_green)
    data["slices"].append(len(red_chan))

#make figures
    fifths=ceil(slices/5)
    exs = [fifths,fifths*2,fifths*3]
    cldu_fig, cldu_axes = plt.subplots(4,3,figsize=(30,40),layout="tight")
    cldu_fig.suptitle("Cldu processing")
    for i in range(3):
        cldu_axes[0,i].imshow(n_green[exs[i]],cmap="gray",vmin=0,vmax=10,interpolation=None)
        cldu_axes[0,i].set_title("Normalized slice " + str(exs[i]))
        
        cldu_axes[1,i].imshow(dn_green[exs[i]],cmap="gray",vmin=0,vmax=10,interpolation=None)
        cldu_axes[1,i].set_title("Median filtered & Normed slice " + str(exs[i]))
        
        cldu_axes[2,i].imshow(gminusr[exs[i]],cmap="gray",vmin=0,vmax=10,interpolation=None) 
        cldu_axes[2,i].set_title("Green with red subtracted slice " + str(exs[i]))
        
        cldu_axes[3,i].imshow(mask_green[exs[i]],cmap="gray",vmin=0,vmax=1,interpolation=None)
        cldu_axes[3,i].set_title("Li thresholded slice " + str(exs[i]))
            
    tub_fig, tub_axes = plt.subplots(3,3,figsize=(30,30),layout="tight")
    tub_fig.suptitle("Tubulin processing")
    for i in range(3):
        tub_axes[0,i].imshow(n_red[exs[i]],cmap="gray",vmin=0,vmax=10,interpolation=None)
        tub_axes[0,i].set_title("Normalized slice " + str(exs[i]))
        
        tub_axes[1,i].imshow(dn_red[exs[i]],cmap="gray",vmin=0,vmax=10,interpolation=None)
        tub_axes[1,i].set_title("Median filtered & Normed slice " + str(exs[i]))
        
        tub_axes[2,i].imshow(mask_red[exs[i]],cmap="gray",vmin=0,vmax=1,interpolation=None)
        tub_axes[2,i].set_title("Li thresholded slice " + str(exs[i]))
    
    plt.savefig("E:/transfection/staining/figures/cldu_processing_"+files[idx])
    plt.savefig("E:/transfection/staining/figures/tubulin_processing_"+files[idx])
    
    
    
out_path="E:/transfection/staining//"
df = pd.DataFrame(data)
df.to_csv(out_path+'cldu_percentage.csv', index=False)




