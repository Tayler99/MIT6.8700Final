import os
import glob
import pandas as pd
import numpy as np

DATA_FOLDER="/n/holyscratch01/stephenson_lab/Users/taylerli/SCENIC/data"

FNAME_556 = os.path.join(DATA_FOLDER, '*_AML556-D*.dem.txt.gz')
names_556 = glob.glob(FNAME_556)
names_556 = [names_556[i] for i in [1, 0, 2]]
aml556_ex_matrix = pd.read_csv(names_556[0], sep='\t', header=0, index_col=0).T
for i in range(1, len(names_556)):
    aml556_ex_matrix_new = pd.read_csv(names_556[i], sep='\t', header=0, index_col=0).T
    aml556_ex_matrix = pd.concat([aml556_ex_matrix, aml556_ex_matrix_new], axis = 0)
    
FNAME_707b = os.path.join(DATA_FOLDER, '*_AML707B-D*.dem.txt.gz')
names_707b = glob.glob(FNAME_707b)
names_707b = [names_707b[i] for i in [1, 0, 3, 2, 4]]
aml707b_ex_matrix = pd.read_csv(names_707b[0], sep='\t', header=0, index_col=0).T
for i in range(1, len(names_707b)):
    aml707b_ex_matrix_new = pd.read_csv(names_707b[i], sep='\t', header=0, index_col=0).T
    aml707b_ex_matrix = pd.concat([aml707b_ex_matrix, aml707b_ex_matrix_new], axis = 0)
    
aml921a_ex_matrix = pd.read_csv(os.path.join(DATA_FOLDER, "GSM3587990_AML921A-D0.dem.txt.gz"), sep='\t', header=0, index_col=0).T

aml556_ex_matrix.to_pickle("BM_AML556_ex_matrix.pkl.gz")
aml707b_ex_matrix.to_pickle("BM_AML707B_ex_matrix.pkl.gz")
aml921a_ex_matrix.to_pickle("BM_AML921A_ex_matrix.pkl.gz")