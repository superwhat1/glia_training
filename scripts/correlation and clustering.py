# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 16:48:23 2023

@author: David
"""


import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
from sklearn.metrics.pairwise import nan_euclidean_distances
import scipy.ndimage as sn
import pandas as pd
    
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
grouped_neurons = pkl.load(open("C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/max proj roi activity/grouped_neurons_responses_by_treatment2023-09-06.pkl", 'rb',))
grouped_glia = pkl.load(open("C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/glia activity/grouped_glia_responses_by_treatment2023-09-15.pkl",'rb'))
grouped_neuropil = pkl.load(open("C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/neuropil activity/grouped_neuropil_responses_by_treatment2023-09-18.pkl",'rb'))
time_list=['min15', 'min45', 'min60', 'min80', 'min100']

def plot_neuron_distance(grouped_by_treatment):
    ready_to_analyze = {}
    output_dir = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/max proj roi activity/average responses/"
    
    for i in time_list:
        ready_to_analyze[i] = {}
        for treatment, animals in grouped_by_treatment.items():
            ready_to_analyze[i][treatment] = {}
            for animal, times in animals.items():
                to_add = []
                for j in range(5):
                    to_stack = []
                    for time, cells in times.items():
                        for cell in cells:
                            if i in time:
                                if treatment !="nocap_notrain": #Temporary code to align older longer traces with newer shorter traces until rois are fixed for certain neuron iscells
                                    to_stack.append(sn.gaussian_filter1d(cell[j, 140:490], sigma = 2) ) #with data of shape cells, n=5 responses, traces of length 350
                                else:
                                    to_stack.append(sn.gaussian_filter1d(cell[j,:], sigma = 2) )
                    try:
                        to_add.append(np.stack(to_stack))
                    except ValueError:
                        pass
                
                ready_to_analyze[i][treatment][animal] = to_add
    
    #EUCLIDIAN DISTANCE method for determining time series similarity
    euc_dist_list = []
    for time, treatments in ready_to_analyze.items():
        for treatment, animals in treatments.items():
            for animal, responses in animals.items():
                for response_ind, response in enumerate(responses):
                    nans_removed = response[~np.isnan(response).all(axis=1)]
                    euc_dist = nan_euclidean_distances(nans_removed)
                    euc_dist_list.append({"time":time, "treatment":treatment, "animal":animal, "response":response_ind, "average euc dist":np.mean(euc_dist[np.tril_indices_from(euc_dist, k=-1)]) })
                    fig, ax = plt.subplots()
                    edm = ax.imshow(euc_dist, cmap = 'RdBu_r')
                    fig.colorbar(edm)
                    ax.set_title(time + treatment + animal + " responses " + str(response_ind))
                    plt.savefig(output_dir + time +treatment + animal + "euc_dist_responses_"+str(response_ind))
    euc_dist_df = pd.DataFrame(euc_dist_list, columns=["time", "treatment", "animal", "response", "average euc dist"])
    
    return ready_to_analyze, euc_dist_df

def plot_glia_distance(grouped_by_treatment):
    ready_to_analyze = {}
    output_dir = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/glia activity/average responses/"
    
    for i in time_list:
        ready_to_analyze[i] = {}
        for treatment, animals in grouped_by_treatment.items():
            ready_to_analyze[i][treatment]={}
            for animal, times in animals.items():
                to_stack = []
                for time, cells in times.items():
                    for cell in cells:
                        for response in cell:
                            if i in time:
                                to_stack.append(sn.gaussian_filter1d(response, sigma = 2) )
                try:
                    stacked = np.stack(to_stack)
                    ready_to_analyze[i][treatment][animal] = stacked[~np.isnan(stacked).all(axis=1)]
                except ValueError:
                    pass

    #EUCLIDIAN DISTANCE method for determining time series similarity
    euc_dist_list = []
    for time, treatments in ready_to_analyze.items():
        for treatment, animals in treatments.items():
            for animal, responses in animals.items():
                euc_dist = nan_euclidean_distances(responses)
                euc_dist_list.append({"time":time, "treatment":treatment, "animal":animal, "average euc dist":np.mean(euc_dist[np.tril_indices_from(euc_dist, k=-1)]) })
                fig, ax = plt.subplots()
                edm = ax.imshow(euc_dist, cmap = 'RdBu_r')
                fig.colorbar(edm)
                ax.set_title(time + treatment + animal+ " glia responses")
                plt.savefig(output_dir + time +treatment + animal + "euc dist glia responses")
        euc_dist_df = pd.DataFrame(euc_dist_list, columns=["time", "treatment", "animal", "average euc dist"])
        
    return ready_to_analyze, euc_dist_df


def plot_neuropil_distance(grouped_by_treatment):
    ready_to_analyze = {}
    output_dir = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/neuropil activity/average responses/"
    
    for i in time_list:
        ready_to_analyze[i] = {}
        for treatment, animals in grouped_by_treatment.items():
            ready_to_analyze[i][treatment]={}
            for animal, times in animals.items():
                to_stack = []
                for time, responses in times.items():
                    for response in responses:
                        if i in time:
                            to_stack.append(sn.gaussian_filter1d(response, sigma = 2) )
                try:
                    stacked = np.stack(to_stack)
                    ready_to_analyze[i][treatment][animal] = stacked[~np.isnan(stacked).all(axis=1)]
                except ValueError:
                    pass

    #EUCLIDIAN DISTANCE method for determining time series similarity
    euc_dist_list = []
    for time, treatments in ready_to_analyze.items():
        for treatment, animals in treatments.items():
            for animal, responses in animals.items():
                euc_dist = nan_euclidean_distances(responses)
                euc_dist_list.append({"time":time, "treatment":treatment, "animal":animal, "average euc dist":np.mean(euc_dist[np.tril_indices_from(euc_dist, k=-1)]) })
                fig, ax = plt.subplots()
                edm = ax.imshow(euc_dist, cmap = 'RdBu_r')
                fig.colorbar(edm)
                ax.set_title(time + treatment + animal+ " glia responses")
                plt.savefig(output_dir + time +treatment + animal + "euc dist neuropil responses")
        euc_dist_df = pd.DataFrame(euc_dist_list, columns=["time", "treatment", "animal", "average euc dist"])
        
    return ready_to_analyze, euc_dist_df


def mean_stack(stacked):
    
    meaned = np.nanmean(stacked, axis=0, dtype = np.float64)
    
    return  meaned


#Function for ploting aligned and averaged responses
def plot_neuron_averaged(grouped_by_treatment):
    output_dir = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/max proj roi activity/average responses/"
    time_list=['min15', 'min45', 'min60', 'min80', 'min100']
    base_min = {}
    base_max = {}
    rolling_amount = {}
    ready_to_plot = {}
    
    for i in time_list:
        ready_to_plot[i] = {}
        for treatment, animals in grouped_by_treatment.items():
            to_stack = []
            for animal, times in animals.items():    
                for time, cells in times.items():
                    for traces in cells:
                        if i in time:
                            to_stack.append(sn.gaussian_filter1d(traces, sigma = 2) )
                            
            stacked = np.stack(to_stack)
            meaned = mean_stack(stacked)
            meaned_again = mean_stack(meaned)
            
            if i == "min15":
                base_min[treatment] = meaned_again.min()
                base_max[treatment] = meaned_again.max()
                rolling_amount[treatment] = 100-np.argmax(meaned_again)    
            normalized = np.roll((meaned_again - base_min[treatment])/ (base_max[treatment] - base_min[treatment]), rolling_amount[treatment])
            ready_to_plot[i][treatment]=normalized
            
        fig, ax = plt.subplots()
        ax.set_title("Average responses for "+i)
        ax.set_ylim(0,1.2)
        ax.set_xlim(0,350)
        ax.set_ylabel('normalized fold change')
        ax.set_xlabel('frames')
        for treatment_name, treatment_data in ready_to_plot[i].items():
            legend_label = treatment_name + "; n=" + str(len(grouped_by_treatment[treatment_name]))
            ax.plot(treatment_data, label =  legend_label)
        ax.legend()
        plt.savefig(output_dir+"average_response_"+i)
        

def plot_neurpil_averaged(grouped_by_treatment):
    output_dir = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/neuropil activity/average responses/"
    time_list=['min15', 'min45', 'min60', 'min80', 'min100']
    base_min = {}
    base_max = {}
    rolling_amount ={}
    ready_to_plot = {}
    
    for i in time_list:
        ready_to_plot[i] = {}
        for treatment, animals in grouped_by_treatment.items():
            to_stack = []
            for animal, times in animals.items():    
                for time, trace in times.items():
                    if i in time:
                        to_stack.append(sn.gaussian_filter1d(trace, sigma = 2) )
                            
            stacked = np.stack(to_stack)
            meaned = mean_stack(stacked)
            meaned_again = mean_stack(meaned)

            
            if i == "min15":
                base_min[treatment] = meaned_again.min()
                base_max[treatment] = meaned_again.max()
                rolling_amount[treatment] = 100-np.argmax(meaned_again)    
            normalized = np.roll((meaned_again - base_min[treatment])/ (base_max[treatment] - base_min[treatment]), rolling_amount[treatment])
            
            ready_to_plot[i][treatment]=normalized
            
        fig, ax = plt.subplots()
        ax.set_title("Average responses for "+i)
        ax.set_ylim(0,1.2)
        ax.set_xlim(0,350)
        ax.set_ylabel('normalized fold change')
        ax.set_xlabel('frames')
        for treatment_name, treatment_data in ready_to_plot[i].items():
            legend_label = treatment_name + "; n=" + str(len(grouped_by_treatment[treatment_name]))
            ax.plot(treatment_data, label =  legend_label)
        ax.legend()
        plt.savefig(output_dir+"average_response_"+i)
        
#PEARSON CORRELATION method for determining time series similarity
def plot_neuron_corr(grouped_by_treatment):
    
    ready_to_analyze = {}
    output_dir = "C:/Users/BioCraze/Documents/Ruthazer lab/glia_training/analysis/max proj roi activity/average responses/"
    
    for i in time_list:
        ready_to_analyze[i] = {}
        for treatment, animals in grouped_by_treatment.items():
            ready_to_analyze[i][treatment] = {}
            for animal, times in animals.items():
                to_add = []
                for j in range(5):
                    to_stack = []
                    for time, cells in times.items():
                        for cell in cells:
                            if i in time:
                                if treatment !="nocap_notrain": #Temporary code to align older longer traces with newer shorter traces until rois are fixed for certain neuron iscells
                                    to_stack.append(sn.gaussian_filter1d(cell[j, 140:490], sigma = 2) ) #with data of shape cells, n=5 responses, traces of length 350
                                else:
                                    to_stack.append(sn.gaussian_filter1d(cell[j,:], sigma = 2) )
                    try:
                        to_add.append(np.stack(to_stack))
                    except ValueError:
                        pass
                
                ready_to_analyze[i][treatment][animal] = to_add
    
    p_corr_list = []
    for time, treatments in ready_to_analyze.items():
        for treatment, animals in treatments.items():
            for animal, responses in animals.items():
                for response_ind, response in enumerate(responses):
                    nans_removed = response[~np.isnan(response).all(axis=1)]
                    if nans_removed.shape[0]>1:
                        p_corr = np.corrcoef(nans_removed)
                        p_corr_list.append({"time":time, "treatment":treatment, "animal":animal, "response":response_ind, "average pcorr":np.mean(p_corr[np.tril_indices_from(p_corr, k=-1)]) })
                    '''
                    fig, ax = plt.subplots()
                    edm = ax.imshow(p_corr, cmap = 'RdBu_r')
                    fig.colorbar(edm)
                    ax.set_title(time + treatment + animal + " responses " + str(response_ind))
                    plt.savefig(output_dir + time +treatment + animal + "euc_dist_responses_"+str(response_ind))
                    '''
    p_corr_df = pd.DataFrame(p_corr_list, columns=["time", "treatment", "animal", "response", "average pcorr"])
    
    return ready_to_analyze, p_corr_df
