# -*- coding: utf-8 -*-
"""
Created on Fri May 23 15:38:26 2025

@author: David
"""

from skimage.io import imread
from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage import filters
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


for idx in range(0,len(files),2):
    print(files[idx])
    red_chan = imread(folder+files[idx])
    green_chan = imread(folder+files[idx+1])
    slices = red_chan.shape[0]
    
    animal_min = np.min([red_chan,green_chan],axis=(0,2,3)).reshape((slices,1,1))
    animal_max = np.max([red_chan,green_chan],axis=(0,2,3)).reshape((slices,1,1))
    n_red = ((red_chan-animal_min)/(animal_max-animal_min))*255
    n_green = ((green_chan-animal_min)/(animal_max-animal_min))*255
    
    #red_min = np.min(red_chan,axis=(1,2)).reshape((slices,1,1))
    #green_min = np.min(green_chan,axis=(1,2)).reshape((slices,1,1))
    #red_max = np.max(red_chan,axis=(1,2)).reshape((slices,1,1))
    #green_max = np.max(green_chan,axis=(1,2)).reshape((slices,1,1))
    
    #c_red = (red_chan-red_min)/(red_max-red_min)
    #c_green = (green_chan-green_min)/(green_max-green_min)
    
    #fn_red = filters.median(n_red)
    #fn_green = filters.median(n_green)
    
    patch_kw = dict(
    patch_size=5,  # 5x5 patches
    patch_distance=6,  # 13x13 search area
    )
    
    sigma_red = np.array(estimate_sigma(n_red,channel_axis=0))
    sigma_green = np.array(estimate_sigma(n_green,channel_axis=0))
    fn_red=np.copy(n_red)
    fn_green=np.copy(n_green)
    for i in range(len(sigma_red)):    
        fn_red[i]= denoise_nl_means(fn_red[i,:,:], h=1.2 * sigma_red[i], sigma=sigma_red[i], fast_mode=False, **patch_kw)
        fn_green[i]=denoise_nl_means(fn_green[i], h=1.2 * sigma_green[i], sigma=sigma_green[i], fast_mode=False, **patch_kw)
    
    gminusr = fn_green-fn_red
    #gminusr[gminusr<0]=0
    
    li_red=np.copy(fn_red)
    li_gminusr=np.copy(gminusr)

    for i in range(slices):
        '''
        r_thresholds = np.linspace(np.min(fn_red[i]) + 0.001, np.mean(fn_red[i]), num=50)
        g_thresholds = np.linspace(np.min(gminusr[i]) + 0.001, np.mean(gminusr[i]), num=50)
        try:
            r_entropies = [_cross_entropy(fn_red[i], t) for t in r_thresholds]
            optimal_red_threshold = r_thresholds[np.nanargmin(r_entropies)]
            red_threshs.append(optimal_red_threshold)
        except:
            try:
                optimal_red_threshold = np.nanmean(red_threshs)
                print(f"failed to compute red entropies for slice {i} using mean of priors instead")
            except:
                optimal_red_threshold=np.max(fn_red)*0.02
                print(f"failed to compute red thresh from prior threshs for slice {i} using 2% of max")
        try: 
            g_entropies = [_cross_entropy(gminusr[i], t) for t in g_thresholds]
            optimal_gminusr_threshold = g_thresholds[np.nanargmin(g_entropies)]
            green_threshs.append(optimal_gminusr_threshold)
        except:
            try:
                optimal_gminusr_threshold = np.nanmean(green_threshs)
                print(f"failed to compute gminusr entropies for slice {i} using mean of priors instead")
            except:
                optimal_gminusr_threshold=np.max(gminusr)*0.02
                print(f"failed to compute gminusr thresh from prior threshs for slice {i} using 2% of max")

        brute_red[i] = fn_red[i]>optimal_red_threshold
        brute_gminusr[i] = gminusr[i]>optimal_gminusr_threshold
        '''
        gminusr_li_thresh = filters.threshold_li(gminusr)
        red_li_thresh = filters.threshold_li(fn_red)
        
        li_red[i]=li_red[i]>red_li_thresh
        li_gminusr[i] = gminusr[i]>gminusr_li_thresh
    
    fifths=ceil(slices/5)
    exs = [fifths,fifths*2,fifths*3]
    cldu_fig, cldu_axes = plt.subplots(4,3,figsize=(30,40),layout="tight")
    cldu_fig.suptitle("Cldu processing")
    for i in range(3):
        cldu_axes[0,i].imshow(n_green[exs[i]],cmap="gray",vmin=0,vmax=10,interpolation=None)
        cldu_axes[0,i].set_title("Normalized slice " + str(exs[i]))
        
        cldu_axes[1,i].imshow(fn_green[exs[i]],cmap="gray",vmin=0,vmax=10,interpolation=None)
        cldu_axes[1,i].set_title("Median filtered & Normed slice " + str(exs[i]))
        
        cldu_axes[2,i].imshow(gminusr[exs[i]],cmap="gray",vmin=0,vmax=10,interpolation=None) 
        cldu_axes[2,i].set_title("Green with red subtracted slice " + str(exs[i]))
        
        cldu_axes[3,i].imshow(li_gminusr[exs[i]],cmap="gray",vmin=0,vmax=1,interpolation=None)
        cldu_axes[3,i].set_title("Li thresholded slice " + str(exs[i]))
            
    tub_fig, tub_axes = plt.subplots(3,3,figsize=(30,30),layout="tight")
    tub_fig.suptitle("Tubulin processing")
    for i in range(3):
        tub_axes[0,i].imshow(n_red[exs[i]],cmap="gray",vmin=0,vmax=10,interpolation=None)
        tub_axes[0,i].set_title("Normalized slice " + str(exs[i]))
        
        tub_axes[1,i].imshow(fn_red[exs[i]],cmap="gray",vmin=0,vmax=10,interpolation=None)
        tub_axes[1,i].set_title("Median filtered & Normed slice " + str(exs[i]))
        
        tub_axes[2,i].imshow(li_red[exs[i]],cmap="gray",vmin=0,vmax=1,interpolation=None)
        tub_axes[2,i].set_title("Li thresholded slice " + str(exs[i]))
    
    plt.savefig("E:/transfection/staining/figures/cldu_processing_"+files[idx])
    plt.savefig("E:/transfection/staining/figures/tubulin_processing_"+files[idx])
    
    percent_green = np.mean(np.sum(li_gminusr,axis=(1,2))/np.sum(li_red,axis=(1,2)))

    data["animal"].append(files[idx])
    data["percent green"].append(percent_green)
    data["slices"].append(len(red_chan))

out_path="E:/transfection/staining//"
df = pd.DataFrame(data)
df.to_csv(out_path+'cldu_percentage.csv', index=False)




