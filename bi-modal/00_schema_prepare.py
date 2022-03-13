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

dataset="pbmc"

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



secondary="ATAC"
if dataset.startswith("cite"):
    secondary = "ADT"

isecondary = secondary.lower()


adata = sc.read_h5ad(f"{dataset}/save/rna.h5ad")
bdata = sc.read_h5ad(f"{dataset}/save/{isecondary}.h5ad")

if isecondary == "adt":
    adata.uns['adt.X'] = bdata.X
    adata.uns['adt.type'] = "ADT"
    adata.uns['adt.var'] = np.array([[idx] + [i] for idx, i in enumerate(bdata.var['features'])])
    adata.uns['adt.var.columns'] = np.array(['id', 'protein'])
    adata.uns['adt.var.index'] = np.array(bdata.var['features'])
    adata.uns['names'] = np.array(['rna', 'adt'])
    adata.uns['rna.type'] = np.array(['gene'])
    n_components = 30 ## cite_elife
    if dataset == "cite_30k":
        n_components = 24 ## cite_elife
else:
    chrs = list(range(1, 22)) ## human or mm need to change
    if dataset == "skin":
        chrs = list(range(1, 19)) ## human or mm need to change
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
    n_components = 50


svd2 = sklearn.decomposition.TruncatedSVD(n_components= n_components, random_state = 17)
H2 = svd2.fit_transform(adata.uns[f"{isecondary}.X"])

sqp99 = schema.SchemaQP(0.99, mode='affine', params= {"decomposition_model":"nmf",
                                                      "num_top_components":n_components,
                                                      "do_whiten": 0,
                                                      "dist_npairs": 5000000})
dz99 = sqp99.fit_transform(adata.X, [H2], ['feature_vector'], [1])


df = pd.DataFrame(dz99)
df.index = adata.obs.index
df.columns = [f"schema_{i+1}" for i in df.columns]
df.to_csv(f"{dataset}/save/schema.csv")
