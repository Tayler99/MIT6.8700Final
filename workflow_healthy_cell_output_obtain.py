# 64GB, 24 cores

import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

DATA_FOLDER="/n/holyscratch01/stephenson_lab/Users/taylerli/SCENIC/data"
RESOURCES_FOLDER="/n/holyscratch01/stephenson_lab/Users/taylerli/SCENIC/data/resources"
DATABASE_FOLDER="/n/holyscratch01/stephenson_lab/Users/taylerli/SCENIC/data/cisTarget/"
OUTPUT_FOLDER="/n/holyscratch01/stephenson_lab/Users/taylerli/SCENIC/data/BM_normal_output"
SCHEDULER="123.122.8.24:8786"
DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "hg19-*.mc9nr.genes_vs_motifs.rankings.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'allTFs_hg38.txt')
SC_EXP_FNAME = os.path.join(DATA_FOLDER, "BM_normal_ex_matrix.pkl.gz")
ADJ_FNAME = os.path.join(OUTPUT_FOLDER, "BM_normal_adjacencies.pkl.gz")
REGULONS_FNAME = os.path.join(OUTPUT_FOLDER, "BM_normal_regulons.p")
MOTIFS_FNAME = os.path.join(OUTPUT_FOLDER, "BM_normal_motifs.csv")
AUCELL_FNAME = os.path.join(OUTPUT_FOLDER, "BM_normal_aucell_values.pkl.gz")

ex_matrix = pd.read_pickle(SC_EXP_FNAME)
tf_names = load_tf_names(MM_TFS_FNAME)

adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)
adjacencies.to_pickle(ADJ_FNAME)

db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
dbs = [dbs[5]]

modules = list(modules_from_adjacencies(adjacencies, ex_matrix))

with ProgressBar():
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

regulons = df2regulons(df)

df.to_csv(MOTIFS_FNAME)
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)
    
auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
auc_mtx.to_pickle(AUCELL_FNAME)