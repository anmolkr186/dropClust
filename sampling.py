import numpy as np
import pandas as pd
import scanpy as sc
import math

# Sampling Primary Clusters
# Desciption: Performs sampling from the primary clusters in an inverse exponential order of cluster size.

# Details: 
# Sampling in inverse proportion of cluster size following a exponential decay equation.
# To ensure selection of sufficient representative transcriptomes from small clusters,
# an exponential decay function  is used to determine the proportion of transcriptomes to be sampled from each
# cluster. For $i^{th}$ cluster, the proportion of expression profiles $p_i$ was {formula from paper}

# Parameters: 
# object: A SingleCellExperiment object containing normalized expression values in \code{"normcounts"}.
# nsamples: integer, total number of samples to return post sampling; ignored when \code{optm_parameters = FALSE}.
# method: character, one of c("sps","random"). Structure Preserving Sampling (sps) selects proportional number of members from each cluster obtained from partitioning an approximate nearest neighbour graph.
# optm_parameters: logical, when TRUE the parameters (\code{pinit, pfin, K}) are optimized such that exactly \code{nsamples} are returned. Optimization is performed using simulated annealing
# pinit: numeric [0,0.5], minimum probability of that sampling occurs from a cluster, ignored when \code{optm_parameters = TRUE}.
# pfin: numeric [0.5,1], maximum probability of that sampling occurs from a cluster, ignored when \code{optm_parameters = TRUE}.
# K: numeric, scaling factor analogous to Boltzmann constant, ignored when \code{optm_parameters = TRUE}.

# Return: A SingleCellExperiment object with an additional column named \code{Sampling} in \code{colData} column. The column stores a a logical value against each cell  to indicate if it has been sampled.

adata = sc.read_10x_mtx('hg19/', var_names='gene_symbols', cache=True)   
print(adata.X)


# objects needs to be an AnnData object 
def sampling(object, nsamples=500, method = "sps", optm_parameters=FALSE, pinit=0.195, pfin = 0.9, K=500)):
    global nsamples
    global method 
    global optm_parameters
    global pinit
    global pfin

    """
    if(nsamples>=ncol(object)) {
    SummarizedExperiment::colData(object)$Sampling=  rep(TRUE, ncol(object))
    return(object)
    """

    if(method not in ["random", "sps"]):
        print("Method not found")
        exit()
    
    no_samples = object.X.shape[1]
    init = no_samples if no_samples < 20000 else min(20000,round(no_samples/3)
    # random sample of ids from sample = 0 to no_samples - 1 of size init 
    sample_ids = np.random.choice(np.array([range(0, no_samples-1,1)]), init) 

    if(method=="sps"):
        """
        if(!any(reducedDimNames(object)=="CComponents"))
        data = Log2Normalize(normcounts(object)[SingleCellExperiment::rowData(object)$HVG, sample_ids],return.sparse = FALSE)
        else
        data = as.matrix(normcounts(object)[, sample_ids])
        # data = as.matrix(assay(object,i = "logcounts")[SingleCellExperiment::rowData(object)$HVG,idx])
        """
        partition = annPartition(data)

        # ---- How to convert partition-----
        # data = pd.Series([1, 1, 1, 2, 3, 3, 3, 3, 4, 4, 5])
        # data.value_counts()


        if(optm_parameters==True):
            param = optimized_param(partition, nsamples)
            pinit = param[1]
            pfin = param[2]
            K = param[3]
            print("Optimized parameters:\n", param,"\n")

        #old seed
        cluster_freq = table(partition[,2])
        prop = round((pinit - exp(-cluster_freq/K) * (pinit - pfin) )* cluster_freq)
        cluster_freq = np.vstack(cluster_freq, prop).T
        
        # prop = reshape2::melt(prop)$value
        
        prt = np.asarray(partition)
        subsample = []

        for i in range(len(prop)):
            subsample = subsample + np.random.choice(prt[np.nonzero(prt[:][2]==i)], size = prop[k], replace = False)

        subsample = np.asarray(subsample)

        """
        
        .Random.seed = oldseed
        SummarizedExperiment::colData(object)$Sampling = rep(FALSE, ncol(object))
        SummarizedExperiment::colData(object)$Sampling[sample_ids[subsamples]] =  TRUE
        """

        print(len(subsamples), "Samples extracted.\n")

    elif(method=="random"):
        """
        oldseed = .Random.seed
        subsamples = sample(sample_ids, nsamples)
        .Random.seed = oldseed
        SummarizedExperiment::colData(object)$Sampling = rep(FALSE, ncol(object))
        SummarizedExperiment::colData(object)$Sampling[subsamples] =  TRUE

        """
    else:
        print("Invalid Sampling. Fallback to all samples")


    # object@metadata[["dropClust"]] = c(unlist(object@metadata[["dropClust"]]),"Sampling")

def annData():

def optimized_param(partition, nsamples=500):
    
