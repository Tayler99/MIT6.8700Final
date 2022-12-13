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
OUTPUT_FOLDER="/n/holyscratch01/stephenson_lab/Users/taylerli/SCENIC/data/BM_AML_normal_output"
SCHEDULER="123.122.8.24:8786"
DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "hg19-*.mc9nr.genes_vs_motifs.rankings.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'allTFs_hg38.txt')
SC_EXP_FNAME = os.path.join(DATA_FOLDER, "BM_*_ex_matrix.pkl.gz")
ADJ_FNAME = os.path.join(OUTPUT_FOLDER, "BM_AML_normal_adjacencies.pkl.gz")
REGULONS_FNAME = os.path.join(OUTPUT_FOLDER, "BM_AML_normal_regulons.p")
MOTIFS_FNAME = os.path.join(OUTPUT_FOLDER, "BM_AML_normal_motifs.csv")
AUCELL_FNAME = os.path.join(OUTPUT_FOLDER, "BM_AML_normal_aucell_values.pkl.gz")

names = glob.glob(SC_EXP_FNAME)
names = [names[i] for i in [0, 2, 1, 3]]

ex_matrix = pd.read_pickle(names[0])
for i in range(1, len(names)):
    ex_matrix_new = pd.read_pickle(names[i])
    ex_matrix = pd.concat([ex_matrix, ex_matrix_new], axis = 0)
tf_names = load_tf_names(MM_TFS_FNAME)

adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)
adjacencies.to_pickle(ADJ_FNAME)

db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
dbs = [dbs[5]]

modules = list(modules_from_adjacencies(adjacencies, ex_matrix))

# Calculate a list of enriched motifs and the corresponding target genes for all modules.
with ProgressBar():
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

# Create regulons from this table of enriched motifs.
regulons = df2regulons(df)

# Save the enriched motifs and the discovered regulons to disk.
df.to_csv(MOTIFS_FNAME)
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)
    
auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
auc_mtx.to_pickle(AUCELL_FNAME)