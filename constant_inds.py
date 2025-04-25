import pandas as pd
import numpy as np

def constant_inds():
    file1 = 'CBB_constants.xlsx'
    tablek = pd.read_excel(file1)
    knames = tablek['name'].values
    k = tablek['value'].values

    file3 = 'CBB_model.xlsx'
    Rknames = pd.read_excel(file3, header=None).iloc[:, 1].values

    kconst = np.zeros(len(Rknames))

    rate_constants_reactions = Rknames

    for i, rate_const in enumerate(knames):
        rate_const_idcs = np.where(rate_constants_reactions == rate_const)[0]
        if rate_const_idcs.size > 0:
            kconst[rate_const_idcs] = i  # Adjusting index to be 1-based

    return k, kconst

