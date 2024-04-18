# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 11:05:53 2023

@author: BioCraze
"""

import numpy as np

stat = np.load("C:/Users/BioCraze/Documents/Ruthazer lab/glia projects/plasticity/data/training_23aug23/A1_min15_23aug23/suite2p/plane0/stat.npy", allow_pickle=True)
F = np.load("C:/Users/BioCraze/Documents/Ruthazer lab/glia projects/plasticity/data/training_23aug23/A1_min15_23aug23/suite2p/plane0/F.npy", allow_pickle=True)
iscell = np.load("C:/Users/BioCraze/Documents/Ruthazer lab/glia projects/plasticity/data/training_23aug23/A1_min15_23aug23/suite2p/plane0/iscell.npy", allow_pickle=True)
spks = np.load("C:/Users/BioCraze/Documents/Ruthazer lab/glia projects/plasticity/data/training_23aug23/A1_min15_23aug23/suite2p/plane0/spks.npy", allow_pickle=True)
fneu = np.load("C:/Users/BioCraze/Documents/Ruthazer lab/glia projects/plasticity/data/training_23aug23/A1_min15_23aug23/suite2p/plane0/Fneu.npy", allow_pickle=True)

fneu0 = np.delete(fneu, [75],0)
F0 = np.delete(F, [75],0)
stat0 = np.delete(stat, [75],0)
iscell0 = np.delete(iscell, [75],0)
spks0 = np.delete(spks, [75],0)

np.save("C:/Users/BioCraze/Documents/Ruthazer lab/glia projects/plasticity/data/training_23aug23/A1_min15_23aug23/suite2p/plane0/fneu.npy", fneu,allow_pickle=True)
np.save("C:/Users/BioCraze/Documents/Ruthazer lab/glia projects/plasticity/data/training_23aug23/A1_min15_23aug23/suite2p/plane0/F.npy", F,allow_pickle=True)
np.save("C:/Users/BioCraze/Documents/Ruthazer lab/glia projects/plasticity/data/training_23aug23/A1_min15_23aug23/suite2p/plane0/iscell.npy", iscell,allow_pickle=True)
np.save("C:/Users/BioCraze/Documents/Ruthazer lab/glia projects/plasticity/data/training_23aug23/A1_min15_23aug23/suite2p/plane0/stat.npy", stat,allow_pickle=True)
np.save("C:/Users/BioCraze/Documents/Ruthazer lab/glia projects/plasticity/data/training_23aug23/A1_min15_23aug23/suite2p/plane0/spks.npy", spks,allow_pickle=True)

