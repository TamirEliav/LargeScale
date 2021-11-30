#%% imports
import pandas as pd
import numpy as np
import xarray as xr
import os
import sys
import exp_decode
import time
from glob import glob
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import normalize

#%%
def read_netcdfs(files, dim):
    # glob expands paths with * to a list of files, like the unix shell
    paths = sorted(glob(files))
    datasets = [xr.open_dataset(p) for p in paths]
    combined = xr.concat(datasets, dim)
    return combined

#%% 
exp_ID = 'b9861_d180527'
param_opt = 1
session_type = 'flight'
exp = exp_decode.exp_decode(exp_ID)
exp.load_data()

#%%
main_dir_IN = "D:\sequences\seq_uri_eden\decoded"
file_template = Path(main_dir_IN) / session_type / exp_ID/ "opt_{}".format(param_opt)
file_template = file_template.glob('b*_d*_res.nc')

#%%
res = xr.open_mfdataset(file_template, parallel=True,concat_dim='time')
# file_template = 'X:\\sequences\\data_decoded\\718383\\b*_d*_flight_*_res.nc'
# res = xr.open_mfdataset(file_template, parallel=True,concat_dim='time')

#%%
res_pos = res['acausal_posterior'].sum('state')
res_pos = res_pos.coarsen(time=10,boundary='trim').mean()
pos_MAP = res_pos.position[res_pos.argmax('position')]

#%%
res_pos.plot(x='time', robust=True, cmap="bone_r")
# plt.plot(pos_MAP.time, pos_MAP,'.r')
# exp.data.position.plot(linewidth=0.5,color='r')
# plt.plot(exp.data.position.data,'.r')

#%% compare
pos_real, pos_pred = xr.align(exp.data.position, pos_MAP, join='left')
plt.scatter(pos_real, pos_pred)

#%%

#%%
n = 1000
x = np.random.randn(n)
y = 0.5*x + np.random.randn(n)*0.1
plt.scatter(x,y)
plt.plot([-3,3],[-3,3],color='k')

#%%
df = pd.DataFrame(x,y)

#%%
bins = np.arange(-3,3,0.1)
cm, bin_edges, bin_edges  = np.histogram2d(x, y, bins=bins)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
cm = normalize(cm, axis=0, norm='l1')
# plt.imshow(cm)
# confusion_matrix
# H = H.T
# plt.imshow()
sns.heatmap(cm,xticklabels=bin_centers,yticklabels=bin_centers)









#%%