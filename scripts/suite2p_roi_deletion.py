# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 11:05:53 2023

@author: BioCraze
"""

import numpy as np

stat = np.load("E:/glia projects/plasticity/data/training_A1_30mar23/A1_min80_30mar23/suite2p/plane0/stat.npy", allow_pickle=True)
F = np.load("E:/glia projects/plasticity/data/training_A1_30mar23/A1_min80_30mar23/suite2p/plane0/F.npy", allow_pickle=True)
iscell = np.load("E:/glia projects/plasticity/data/training_A1_30mar23/A1_min80_30mar23/suite2p/plane0/iscell.npy", allow_pickle=True)
spks = np.load("E:/glia projects/plasticity/data/training_A1_30mar23/A1_min80_30mar23/suite2p/plane0/spks.npy", allow_pickle=True)
fneu = np.load("E:/glia projects/plasticity/data/training_A1_30mar23/A1_min80_30mar23/suite2p/plane0/Fneu.npy", allow_pickle=True)

fneu0 = np.delete(fneu, [120,31],0)
F0 = np.delete(F, [120,31],0)
stat0 = np.delete(stat, [120,31],0)
iscell0 = np.delete(iscell, [120,31],0)
spks0 = np.delete(spks, [120,31],0)

np.save("E:/glia projects/plasticity/data/training_A1_30mar23/A1_min80_30mar23/suite2p/plane0/fneu.npy", fneu0,allow_pickle=True)
np.save("E:/glia projects/plasticity/data/training_A1_30mar23/A1_min80_30mar23/suite2p/plane0/F.npy", F0,allow_pickle=True)
np.save("E:/glia projects/plasticity/data/training_A1_30mar23/A1_min80_30mar23/suite2p/plane0/iscell.npy", iscell0,allow_pickle=True)
np.save("E:/glia projects/plasticity/data/training_A1_30mar23/A1_min80_30mar23/suite2p/plane0/stat.npy", stat0,allow_pickle=True)
np.save("E:/glia projects/plasticity/data/training_A1_30mar23/A1_min80_30mar23/suite2p/plane0/spks.npy", spks0,allow_pickle=True)
