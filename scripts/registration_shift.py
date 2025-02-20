# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 14:30:21 2024

@author: David
"""

import numpy as np
from os import scandir, path, makedirs
import glob, re
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from tifffile import imwrite
from re import search
from pystackreg import StackReg
from skimage import transform, io, exposure
import pandas as pd
from math import exp

def composite_images(imgs, equalize=False, aggregator=np.mean):

    if equalize:
        imgs = [exposure.equalize_hist(img) for img in imgs]

    imgs = [img / img.max() for img in imgs]

    if len(imgs) < 3:
        imgs += [np.zeros(shape=imgs[0].shape)] * (3-len(imgs))

    imgs = np.dstack(imgs)

    return imgs

def overlay_images(imgs, equalize=False, aggregator=np.mean):

    if equalize:
        imgs = [exposure.equalize_hist(img) for img in imgs]

    imgs = np.stack(imgs, axis=0)

    return aggregator(imgs, axis=0)


folders = [i.path for i in scandir('E:\glia projects\plasticity\data')]
output_dir = "E:/glia projects/plasticity/analysis/registration/"
roi_db={}
img_db={}
grid_search={"transform":[],"percent":[],"vlim":[],"cell_count":[]}
transformations = {
    'TRANSLATION': StackReg.TRANSLATION,
    'RIGID_BODY': StackReg.RIGID_BODY,
    'SCALED_ROTATION': StackReg.SCALED_ROTATION,
    'AFFINE': StackReg.AFFINE,
    'BILINEAR': StackReg.BILINEAR
}
for i, (name, tf) in enumerate(transformations.items()):
    for percent in range(0,50,10):
        for vlim in range(0,5):
            total_count=0
            for folder in folders[:5]:
                animal = re.search("A\d{1}_\d{2}\D{3}\d{2}", folder).group()
                files =[f.path[:f.path.rfind('\\')+1] for i in glob.glob(f'{folder}/*/**/***/',recursive=True) for f in scandir(i) if f.path.endswith('ops.npy') and "capapplication" not in f.path]
                files=files[[idx for idx, file in enumerate(files) if 'min15' in file][0]:]+files[:[idx for idx, file in enumerate(files) if 'min15' in file][0]]
                
                to_reg=[]
                roi_coords=pd.DataFrame()
                for file in files:
                    raw_img = np.load(file + "ops.npy",allow_pickle=True).item()["meanImg"]
                    #normed_img = (raw_img-np.min(raw_img))/(np.max(raw_img)-np.min(raw_img))
                    
                    normed_img = np.copy(raw_img)
                    norm=Normalize(vmin=np.min(raw_img),vmax=np.min(raw_img)+(vlim*0.1)*np.max(raw_img),clip=True)
                    normed_img[raw_img<np.percentile(raw_img,percent)]=0
                    
                    to_reg.append(norm(normed_img))
                    stat = np.load(file+"stat.npy",allow_pickle=True)
                    roi_coords["med_"+re.search("min\d{3}|min\d{2}", file).group()]=pd.Series([stat[i]["med"]for i in range(len(stat))])
                    roi_coords["radius_"+re.search("min\d{3}|min\d{2}", file).group()]=pd.Series([stat[i]["radius"] for i in range(len(stat))])
                    roi_coords["id"+re.search("min\d{3}|min\d{2}",file).group()]=pd.Series([])
                
                to_reg = np.stack(to_reg,axis=0)
                sr = StackReg(tf)
        
                reference = 'first' if name == 'BILINEAR' else 'previous'
        
                tmat = sr.register_stack(to_reg, axis=0, reference=reference)
                
                tr_coords=roi_coords.copy(deep=True)
                #sr = StackReg(StackReg.AFFINE)
                #tmat = sr.register_stack(to_reg,reference='previous')
                img_db[animal] = sr.transform_stack(to_reg)
                
                for img in range(int(len(roi_coords.columns)/3)):
                    c_name = roi_coords.columns[img*3]
                    tr_coords.insert(tr_coords.columns.get_loc(c_name)+1,"trsfrmd_"+c_name, pd.Series(transform.matrix_transform(np.array([i for i in roi_coords.loc[:,c_name] if ~np.any(np.isnan(i))]),tmat[img]).tolist())) #using np.any because i is a coord pair [x,y]
                    
                for trcs in range(int(len(tr_coords.columns)/4)):
                    idx=trcs*4+1
                    if idx < int(len(tr_coords.columns)-3):     #find the last time point and skip comparing to nothing as it is the last one
                        src=np.array([i for i in tr_coords.iloc[:,idx] if ~np.any(np.isnan(i))])
                        dst=np.array([i for i in tr_coords.iloc[:,idx+4] if ~np.any(np.isnan(i))])
                        dst_matrix = np.linalg.norm(src[:,None] - dst,axis=-1)
                        #src is reshaped to be [N,1,2] in order to be braodcastable with dst [M,2] and then linalg.norm computes the Frobenius norm to finish calculating the eucl distance (i.e. sqrt(sum([xdiff,ydiff]**2)))
                        dst_matrix[dst_matrix>=np.array(tr_coords.iloc[:,idx+1][~np.isnan(tr_coords.iloc[:,idx+1])])[:,None]]=np.inf #all distances above the radius of the current time point for each cell is set to infinity. 
                        match_matrix=np.full_like(dst_matrix,np.inf)
                        
                        for src in range(dst_matrix.shape[0]):
                            for dst in range(dst_matrix.shape[1]):
                                if [src,np.argmin(dst_matrix[src,:])]==[np.argmin(dst_matrix[:,dst]),dst]:
                                    match_matrix[src,dst]=1 #check each value in the array to determine if it the lowest value in its row and column, if so that the coresponding index in a masking array get set to one, otherwise it's inf.
                        fltrd_dst_matrix = dst_matrix*match_matrix #do the element wise product of the distance matrix and the mask array to get the unique closest matches.
            
                        for src in range(fltrd_dst_matrix.shape[0]):
                            id_col = trcs*4+3
                            if "min15" in tr_coords.columns[id_col]:
                                if np.any(fltrd_dst_matrix[src,:]!=np.inf):
                                    dst = np.argmin(fltrd_dst_matrix[src,:])
                                    tr_coords.iloc[src,id_col]=src
                                    tr_coords.iloc[dst,id_col+4]=src
                            elif id_col<int(len(tr_coords.columns))-4:
                                if np.any(fltrd_dst_matrix[src,:]!=np.inf) and tr_coords.iloc[src,id_col]!=np.nan:
                                    tr_coords.iloc[dst,id_col+4]=tr_coords.iloc[src,id_col]
                
                roi_db[animal]=tr_coords.loc[:,[i for i in tr_coords.columns if "idmin" in i]]
                total_count+=tr_coords.iloc[:,-1].count()
            grid_search["transform"].append(name)
            grid_search["percent"].append(percent)
            grid_search["vlim"].append(vlim)
            grid_search["cell_count"].append(total_count)
            print("tf: "+str(name)+"\n", "percent: "+str(percent)+"\n","vlim: "+str(vlim)+"\n", "count"+str(total_count))

'''
    fig,axes =plt.subplots(nrows=1, ncols=5, figsize = (30,10), layout="constrained")
    axes=axes.ravel()
    for i in range(out.shape[0]):
        if i==0:
            axes[i].imshow(out[0], cmap='gray', vmin=0,vmax=0.15)
        else:
            axes[i].imshow(out[i], cmap='gray',vmin=0,vmax=0.15)
    
'''

#############################################################################################################################################
#TESTING BEST REGISTRATION METHOD
transformations = {
    'TRANSLATION': StackReg.TRANSLATION,
    'RIGID_BODY': StackReg.RIGID_BODY,
    'SCALED_ROTATION': StackReg.SCALED_ROTATION,
    'AFFINE': StackReg.AFFINE,
    'BILINEAR': StackReg.BILINEAR
}

f, ax = plt.subplots(2, int(np.ceil((len(transformations)+1)/2)), figsize=(20, 12))
ax = ax.ravel()

ax[0].imshow(overlay_images(to_reg, aggregator=np.mean), cmap='gray', vmin=0, vmax=0.04)
ax[0].set_title('Original (overlay)')
ax[0].axis('off')

# store transformation matrices for later use in this variable
tmats = []

for i, (name, tf) in enumerate(transformations.items()):
    sr = StackReg(tf)

    reference = 'first' if name == 'BILINEAR' else 'previous'

    tmat = sr.register_stack(to_reg, axis=0, reference=reference, verbose=True)
    reg = sr.transform_stack(to_reg)

    tmats.append(tmat)

    ax[i+1].imshow(overlay_images(reg, aggregator=np.mean), cmap='gray', vmin=0, vmax=0.04)
    ax[i+1].set_title(name + ' (overlay)')
    ax[i+1].axis('off')

#AFFINE LOOKS BEST

f, ax = plt.subplots(1,3, figsize=(20, 12))
ax = ax.ravel()

ax[0].imshow(overlay_images(to_reg, aggregator=np.mean), cmap='gray', vmin=0, vmax=0.04)
ax[0].set_title('Original (overlay)')
ax[0].axis('off')

ar=StackReg(StackReg.AFFINE)
tmat = ar.register_stack(to_reg, axis=0, reference='previous', verbose=True)
reg_ap = ar.transform_stack(to_reg)
ax[1].imshow(overlay_images(reg_ap, aggregator=np.mean), cmap='gray', vmin=0, vmax=0.04)
ax[1].set_title('AFFINE previous')
ax[1].axis('off')


tmat = ar.register_stack(to_reg, axis=0, reference='first', verbose=True)
reg_af = ar.transform_stack(to_reg)
ax[2].imshow(overlay_images(reg_af, aggregator=np.mean), cmap='gray', vmin=0, vmax=0.04)
ax[2].set_title('AFFINE first')
ax[2].axis('off')

plt.show()

f, ax =plt.subplots(2,5, figsize=(20,12))
ax=ax.ravel()
for i in range(len(ax)):
    if i <5:
        ax[i].imshow(reg_ap[i],cmap="gray", vmin=0,vmax=0.04)
    else:
        ax[i].imshow(reg_af[i-5],cmap="gray",vmin=0,vmax=0.04)
plt.show()

#AFFINE PREVIOUS LOOKS BETTER THAN AFFINE FIRST


test = np.copy(raw_img)
norm=Normalize(vmin=np.min(test),vmax=np.min(test)+np.max(test)*0.1,clip=True)
test[raw_img<np.percentile(raw_img,45)]=0
plt.imshow(norm(normed_img),cmap='gray')
    
plt.imshow(normed_img,cmap='gray')
