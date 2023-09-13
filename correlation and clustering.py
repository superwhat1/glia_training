# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 16:48:23 2023

@author: David
"""


import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl

#K-SHAPE METHOD for clustering time series
'''
from tslearn.clustering import KShape
from tslearn.datasets import CachedDatasets
from tslearn.preprocessing import TimeSeriesScalerMeanVariance

seed = 0
np.random.seed(seed)
X_train, y_train, X_test, y_test = CachedDatasets().load_dataset("Trace")
# Keep first 3 classes and 50 first time series

X_train = X_train[:50]
np.random.shuffle(X_train)
# For this method to operate properly, prior scaling is required
X_train = TimeSeriesScalerMeanVariance().fit_transform(X_train)
sz = X_train.shape[1]

# Instantiate k-Shape model
ks = KShape(n_clusters=3, verbose=True, random_state=seed)

# Train
ks.fit(X_train)

# Save model
ks.to_hdf5('./ks_trained.hdf5')

# Load model
trained_ks = KShape.from_hdf5('./ks_trained.hdf5')

y_pred = trained_ks.predict(X_test)

plt.figure()
for yi in range(3):
    plt.subplot(3, 1, 1 + yi)
    for xx in  X_test[y_pred == yi]:
        plt.plot(xx.ravel(), "k-", alpha=.2)
    plt.plot(ks.cluster_centers_[yi].ravel(), "r-")
    plt.xlim(0, sz)
    plt.ylim(-4, 4)
    plt.title("Cluster %d" % (yi + 1))

plt.tight_layout()
plt.show()
'''

#import the Responses dictionary
grouped_by_treatment = pkl.load(open("C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/max proj roi activity/grouped_neurons_responses_by_treatment2023-09-06.pkl", 'rb',))

time_list=['min15', 'min45', 'min60', 'min80', 'min100']
ready_to_analyze = {}

for i in time_list:
    ready_to_analyze[i] = {}
    for treatment, animals in grouped_by_treatment.items():
        to_add = []
        for j in range(5):
            to_stack = []
            for animal, times in animals.items():    
                for time, cells in times.items():
                    for cell in cells:
                        if i in time:
                            to_stack.append(cell[j,:]) #with data of shape cells, n=5 responses, traces of length 350
                try:
                    to_add.append(np.stack(to_stack))
                except ValueError:
                    pass
        ready_to_analyze[i][treatment] = to_add

#EUCLIDIAN DISTANCE method for determining time series similarity
from sklearn.metrics.pairwise import nan_euclidean_distances
import matplotlib as mpl

for times, treatments in ready_to_analyze.items():
    for treatment, responses in treatments.items():
        for response in np.arange(stop = len(responses), step = 5):
            concatenated = np.concatenate(responses[response:response+5])
            nans_removed = concatenated[~np.isnan(concatenated).all(axis=1)]
            
            euc_dist = nan_euclidean_distances(nans_removed)
            #df = pd.DataFrame(euc_dist, index=recording_name, columns=recording_name)
            fig,ax = plt.subplots()
            edm = ax.imshow(euc_dist, cmap = 'RdBu_r')
            fig.colorbar(edm)
            ax.set_title(times + treatment + " responses " + str(response))


#PEARSON CORRELATION method for determining time series similarity
'''
corr = np.corrceff()
fig,ax = plt.subplot()
ax.imshow(corr, vmin = -1, vmax = 1)
'''