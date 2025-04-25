import numpy as np

def ODEs(t, y, krxn, k, rate_inds, S, Rknames, species, kidcs):
    cidcs = [0, 2, 4, 5]
    gapidx = 7

    nrxn = len(rate_inds)
    r = np.zeros(nrxn)

    for irxn in range(nrxn):
        r[irxn] = krxn[irxn] * np.prod(y[np.array(rate_inds[irxn])])

    r[6] = r[1] / 3
    dydt = S @ r

    return dydt

