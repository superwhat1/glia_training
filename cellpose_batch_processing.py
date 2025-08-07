# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:08:24 2023

@author: David
"""

from cellpose import models
from cellpose.io import imread
import glob

# Put the PATH to your trained model
model = models.Cellpose(pretrained_model = 'MODEL PATH')

# PUT PATH TO YOUR FILES HERE!
rootdir = "FOLDER WITH FILES"
# list of files
files = [glob.glob(f'{rootdir}/*/**/*.tif', recursive=True)]

# list of images
imgs = [imread(f) for f in files]
nimg = len(imgs)

# define CHANNELS to run segementation on
# grayscale=0, R=1, G=2, B=3
# channels = [cytoplasm, nucleus]
# if NUCLEUS channel does not exist, set the second channel to 0
channels = [[0,0]]
# IF ALL YOUR IMAGES ARE THE SAME TYPE, you can give a list with 2 elements
# channels = [0,0] # IF YOU HAVE GRAYSCALE
# channels = [2,3] # IF YOU HAVE G=cytoplasm and B=nucleus
# channels = [2,1] # IF YOU HAVE G=cytoplasm and R=nucleus

# if diameter is set to None, the size of the cells is estimated on a per image basis
# you can set the average cell `diameter` in pixels yourself (recommended)
# diameter can be a list or a single number for all images

masks, flows, styles, diams = model.eval(imgs, diameter=None, channels=channels)
