#!/usr/bin/env python
# coding: utf-8
import sklearn.decomposition
import scanpy as sc
import anndata
import schema
import schema.utils
import pandas as pd
import numpy as np
import getopt
import sys
import re

dataset="TEA"

try:
    #options,args = getopt.getopt(sys.argv[1:],"dn")
    options,args = getopt.getopt(sys.argv[1:],"d:")
except getopt.GetoptError:
    print("Erorr Parametes")
    sys.exit()
for name,value in options:
    if name in "-d":
        dataset = value
        print("dataest is: ", value)




adata = sc.read_h5ad(f"{dataset}/save/rna.h5ad")
bdata = sc.read_h5ad(f"{dataset}/save/atac.h5ad")
cdata = sc.read_h5ad(f"{dataset}/save/adt.h5ad")

#### add 1st
chrs = list(range(1, 22))
if dataset == "skin":
    chrs = list(range(1, 19))

chrs.extend(["X", "Y"])
chrs = [f"chr{i}" for i in chrs]
if dataset == "pbmc":
    index = [i for i in bdata.var.index if re.split(':|-', i)[0] in chrs and len(re.split(':|-', i)) == 3]
else:
    index = [i for i in bdata.var.index if i.split('-')[0] in chrs and len(i.split('-')) == 3]
bdata._inplace_subset_var(index)

adata.uns['atac.X'] = bdata.X
adata.uns['atac.type'] = "peak"
adata.uns['atac.var'] = np.array([[idx] + re.split(":|-", i) for idx, i in enumerate(bdata.var['features'])])
adata.uns['atac.var.columns'] = np.array(['id', 'chr', 'start', 'end'])
adata.uns['atac.var.index'] = np.array(bdata.var['features'])
adata.uns['names'] = np.array(['rna', 'atac'])
adata.uns['rna.type'] = np.array(['gene'])


#### add 2nd
adata.uns['adt.X'] = cdata.X
adata.uns['adt.type'] = "ADT"
adata.uns['adt.var'] = np.array([[idx] + [i] for idx, i in enumerate(cdata.var['features'])])
adata.uns['adt.var.columns'] = np.array(['id', 'protein'])
adata.uns['adt.var.index'] = np.array(cdata.var['features'])
adata.uns['names'] = np.array(['rna', 'adt'])
adata.uns['rna.type'] = np.array(['gene'])


n_components_atac = 50
svd2 = sklearn.decomposition.TruncatedSVD(n_components= n_components_atac, random_state = 17)
H2 = svd2.fit_transform(adata.uns["atac.X"])


n_components_adt = 30
svd3 = sklearn.decomposition.TruncatedSVD(n_components= n_components_adt, random_state = 17)
H3 = svd3.fit_transform(adata.uns["adt.X"])



sqp99 = schema.SchemaQP(0.99, mode='affine', params= {"decomposition_model":"nmf",
                                                      "num_top_components":min(n_components_adt, n_components_atac),
                                                      "do_whiten": 0,
                                                      "dist_npairs": 5000000})
dz99 = sqp99.fit_transform(adata.X, [H2, H3], ['feature_vector', 'feature_vector'], [1, 1])


#embed = sqp99._decomp_mdl.transform(adata.X.todense())


df = pd.DataFrame(dz99)
df.index = adata.obs.index
df.columns = [f"schema_{i+1}" for i in df.columns]
df.to_csv(f"{dataset}/save/schema.csv")
