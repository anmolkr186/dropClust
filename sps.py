import numpy as np
import sys
import pandas as pd
import scanpy as sc
import math
from GenSA import *
# import scipy
from scipy.optimize import dual_annealing
from scipy.sparse import csr_matrix
from annoy import AnnoyIndex
import random
from igraph import *
from sklearn.preprocessing import normalize

def annPartition(data):
    f = data.shape[0]
    # print("datashpae",data.shape)
    print("Building graphs with ", data.shape[1], " nodes...")
    t = AnnoyIndex(f, 'angular')
    for i in range(data.shape[1]):
        v = data[:,i]
        t.add_item(i, v)
    
    t.build(30)
    
    get_nn = lambda x: t.get_nns_by_item(x, 6)

    indices = list(map(get_nn, np.arange(data.shape[1])))
    
    indices = np.array(indices)
    
    fin = []
    for i in range(indices.shape[0]):
        for j in indices[i]:
            fin.append((i,j))
    fin = np.array(fin)
    
    g = Graph(fin)
    G = g.simplify(multiple=False, loops = False)
    
    print("Louvain Partition...")
    partition = G.community_leiden(objective_function = "modularity")    
    dataMatrix = np.c_[np.arange(data.shape[1]),np.array(partition.membership)]
    print("Partitioning Done.")

    return dataMatrix

def optimized_param(partition, nsamples = 500):
    if(nsamples> partition.shape[0]):
        nsamples = partition.shape[0]
    
    global_min = 0
    tol = 1e-3
    max_time = 20
    lower = np.array([0.05, 0.9, 500])
    upper = np.array([0.1, 0.95, 4000])
    
    params = None
    
    def pin_find(params, max_c = nsamples):
        output = np.array(partition)
        pinit = params[0]
        pfin = params[1]
        K = params[2]
        
        unique_elements, counts_elements = np.unique(output[:,1], return_counts=True)
        cluster_freq = np.asarray((counts_elements), dtype = int)
        prop = np.round((pinit - np.exp(-cluster_freq/K) * (pinit - pfin) )* cluster_freq)
        cluster_freq = np.vstack((cluster_freq,prop)).T
                
        subsamples_lovain = np.empty((0))
        
        for i in range(len(prop)):
            subsamples_lovain = np.concatenate((subsamples_lovain, np.random.choice( output[output[:,1]==i,0], size = int(prop[i]), replace = False)), axis = None)

        return np.abs(max_c - subsamples_lovain.shape[0]) 

    #out = gensa(func = pin_find, x0 = params, bounds = list(zip(lower, upper)))
    out = dual_annealing(func = pin_find, x0 = params, bounds = list(zip(lower, upper)))
    
    # print(out)    
    # returning best set of parameters
    return out.x

        
def sampling(adata, axis = 0, nsamples=500, method = "sps", optm_parameters=True, pinit=0.195, pfin = 0.9, K=500):
    
    ob = adata.X
    ob = csr_matrix.toarray(ob)
    #sampling rows
    if(axis == 0):
        ob = ob.T
    # print(ob.shape)
    if(nsamples>=ob.shape[1]):
        print("Number of samples are greater than number of columns. Sampling cant be done")
        exit(0)
    
    no_samples = ob.shape[1]
    init = no_samples if no_samples < 20000 else min(20000,round(no_samples/3))

    # random sample of ids from sample = 0 to no_samples - 1 of size init
    sample_ids = np.random.choice(list(range(0, no_samples,1)), init) 

    data = normalize(ob)
    data = np.take(ob, sample_ids, axis = 1)
    
    partition = annPartition(data)

    if(optm_parameters==True):
        param = optimized_param(partition, nsamples)
        pinit = param[0]
        pfin = param[1]
        K = param[2]
        print("Optimized parameters: ", param,"\n")

    unique_elements, counts_elements = np.unique(partition[:,1], return_counts=True)
    cluster_freq = np.asarray((counts_elements), dtype = int)
    # print(cluster_freq.shape)
    prop = np.round((pinit - np.exp(-cluster_freq/K) * (pinit - pfin) )* cluster_freq)
    cluster_freq = np.vstack((cluster_freq,prop)).T
    subsamples = np.empty((0))
    
    for i in range(len(prop)):
        subsamples = np.concatenate((subsamples, np.random.choice(partition[partition[:,1]==i,0], size = int(prop[i]), replace = False)), axis = None)

    subsamples = np.asarray(subsamples, dtype = int)

    print(len(subsamples), "Samples extracted. Returning indices of samples")

    return subsamples
    