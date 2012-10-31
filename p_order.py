import pandas as pd
from pandas import DataFrame, Series
from sys import argv,exit

if len(argv)==1:
    print("please supply a file name")
    exit()

data = pd.read_table(argv[1], header=None)
data.columns = ['m']
tmp = (data['m']<1).value_counts()
if True not in tmp:
    print(1.0)
elif False not in tmp:
    print(0.0)
else:
    print((1.0*(tmp.sum()-tmp[True])/tmp.sum()))
