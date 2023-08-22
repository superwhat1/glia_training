# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 16:48:23 2023

@author: David
"""

#K-SHAPE METHOD for clustering time series
import numpy
import matplotlib.pyplot as plt

from tslearn.clustering import KShape
from tslearn.datasets import CachedDatasets
from tslearn.preprocessing import TimeSeriesScalerMeanVariance

seed = 0
numpy.random.seed(seed)
X_train, y_train, X_test, y_test = CachedDatasets().load_dataset("Trace")
# Keep first 3 classes and 50 first time series

X_train = X_train[:50]
numpy.random.shuffle(X_train)
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


#EUCLIDIAN DISTANCE method for determining time series similarity
from sklearn.metrics.pairwise import euclidean_distances

euc_dist = euclidean_distances(data)
pd.DataFrame(euc_dist, index=recording_name, columns=recording_name)

#PEARSON CORRELATION method for determining time series similarity
np.corrceff(x,y)