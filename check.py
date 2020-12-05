from sps import * 
from pca import *
import scanpy as sc

adata = sc.read_10x_mtx('hg19/', var_names='gene_symbols', cache=True)   
print(adata)
# print(adata.shape)
# print(adata.var.shape)

# if axis = 0, then rows will be samples, if axis = 1 then columns will be sampled
out = sampling(adata, axis = 0)

print(out)

ob = adata.X
ob = csr_matrix.toarray(ob)
ob = np.take(ob, out, axis = 0)
print("Sampled matrix", ob.shape)

reduced_data,evr,sv = pca(ob, 200)
print(reduced_data.shape)

