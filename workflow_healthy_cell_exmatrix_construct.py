import os
import glob
import pandas as pd
import numpy as np

DATA_FOLDER="/n/holyscratch01/stephenson_lab/Users/taylerli/SCENIC/data"
SC_EXP_FNAME = os.path.join(DATA_FOLDER, "GSM35*_BM*.dem.txt.gz")

names = glob.glob(SC_EXP_FNAME)

bm_names = [names[i] for i in [0, 5, 4, 2, 1]]
bm_ex_matrix = pd.read_csv(bm_names[0], sep='\t', header=0, index_col=0).T
for i in range(1, len(bm_names)):
    bm_ex_matrix_new = pd.read_csv(bm_names[i], sep='\t', header=0, index_col=0).T
    bm_ex_matrix = pd.concat([bm_ex_matrix, bm_ex_matrix_new], axis = 0)
bm5_pm_ex_matrix = pd.read_csv(names[3], sep = '\t', header = 0, index_col=0).T
bm5_pm_anno = pd.read_csv(os.path.join(DATA_FOLDER, "GSM3588003_BM5-34p38n.anno.txt.gz"), sep = '\t', header = 0, index_col = 0)

kept_cells = bm5_pm_anno[bm5_pm_anno['CellType'].isnull() == False]
kept_cells_names = list(kept_cells.index)
bm5_pm_kept = bm5_pm_ex_matrix.loc[kept_cells_names]
bm_ex_matrix = pd.concat([bm_ex_matrix, bm5_pm_kept], axis = 0)

bm_ex_matrix.to_pickle("BM_normal_ex_matrix.pkl.gz")