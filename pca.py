## Enter the Gene Expression Matrix After SPS
## Rows should be cells
## Columns should be Gene Expression Values

import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def pca(data,n_components):

    sc = StandardScaler()
    sc_data = sc.fit_transform(data)
    pc = PCA(n_components=n_components)
    reduced_data = pc.fit_transform(sc_data)
   
    evr = pc.explained_variance_ratio_
    sv = pc.singular_values_
    return reduced_data,evr,sv
    
#output values are the reduced data, explained variance ratio of the PCs, and Singular values
#use reduced data for further clustering
#evr and sv are just supplementary data, to plot later if required
