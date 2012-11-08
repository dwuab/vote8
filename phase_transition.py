import pandas as pd
from pandas import DataFrame, Series
from numpy import linspace

N = 900
K = N
#p_list=[0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99]
#p_list=[0.96,0.97,0.98,0.983,0.985,0.987,0.99,0.995]
p_list = [0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99]
eta = 20

p_order = {p:0.0 for p in p_list}

for p in p_list:
    fname = 'data8/N%d_K_%d_p_%.3f_eta_%d_m.dat' % (N,K,p,eta)
    data = pd.read_table(fname, header=None)
    data.columns = ['m']
    tmp = (data['m']<1).value_counts()
    if True not in tmp:
        p_order[p] = 1.0
    elif False not in tmp:
        p_order[p] = 0.0
    else:
        p_order[p] = (1.0*(tmp.sum()-tmp[True])/tmp.sum())

Series(p_order).plot()
