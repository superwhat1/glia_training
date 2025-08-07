# -*- coding: utf-8 -*-
"""
Created on Fri May 23 12:57:48 2025

@author: David
"""
from os import listdir, rename
from numpy import random
folder = 'E:\cellpose_model_training_data\models\immuno_model\\'
rng=random.default_rng()
nums = rng.choice( 9999,len(listdir(folder)), replace=False)
for idx, file in enumerate(listdir(folder)):
    old = folder+file
    new = folder+str(nums[idx])+file
    rename(src=old,dst=new)
