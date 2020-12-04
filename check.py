from sps import * 
import scanpy as sc

adata = sc.read_10x_mtx('hg19/', var_names='gene_symbols', cache=True)   
# if axis = 0, then rows will be samples, if axis = 1 then columns will be sampled
out = sampling(adata, axis = 1)

print(out)