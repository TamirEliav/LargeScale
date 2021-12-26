# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 10:14:06 2021

@author: tamire
"""

#%% imports
import pandas as pd
import numpy as np
import xarray as xr
import os
import sys
import time
from glob import glob
from pathlib import Path

#%%
def join_jobs_results(path):
    file_template = path.glob('b*_d*_res.nc')
    res = xr.open_mfdataset(file_template, parallel=True,concat_dim='time')
    return res

#%% general definitions
dir_IN = "X:\\sequences\\decoding"
# dir_IN = "X:\\sequences\\decoding_without_likelihood_saved"
dir_OUT = "F:\sequences\decoded"
exp_IDs = [
    # 'b9861_d180519',
    # 'b9861_d180521',
    # 'b9861_d180522',
    # 'b9861_d180523',
    # 'b9861_d180524',
    # 'b9861_d180525',
    # 'b9861_d180526',
    # 'b9861_d180527',
    # 'b9861_d180529',
    # 'b9861_d180530',
    
    # 'b0184_d191128',
    # 'b0184_d191129',
    # 'b0184_d191130',
    # 'b0184_d191201',
    # 'b0184_d191202',
    # 'b0184_d191203',
    # 'b0184_d191204',
    # 'b0184_d191205',
    'b0184_d191208',
    # 'b0184_d191209',
    # 'b0184_d191210',
    # 'b0184_d191211',
    # 'b0184_d191212',
    # 'b0184_d191215',
    # 'b0184_d191216',
    # 'b0184_d191217',
    # 'b0184_d191220',
    # # 'b0184_d191222',
    # 'b0184_d191224',
    # 'b0184_d191225',
    # 'b0184_d191226',
    # 'b0184_d200101',
    # 'b0184_d200102',
    
    # 'b0034_d180305',
    # 'b0034_d180306',
    # 'b0034_d180308',
    # 'b0034_d180310',
    # 'b0034_d180311',
    # 'b0034_d180312',
    # 'b0034_d180313',
    # 'b0034_d180314',
    # 'b0034_d180315',
    
    # 'b2289_d180518',
    # 'b2289_d180520',
    
    # 'b0148_d170606', 
    # 'b0148_d170607', 
    # 'b0148_d170608',  
    # 'b0148_d170625',
    # 'b0148_d170626', 
    # 'b0148_d170627', 
    # 'b0148_d170703', 
    # 'b0148_d170718',  
    # 'b0148_d170720', 
    # 'b0148_d170723', 
    # 'b0148_d170801', 
    # 'b0148_d170802', 
    # 'b0148_d170806', 
    # 'b0148_d170807',
    ]
session_type = 'flight'
# session_type = 'sleep'
# session_type = 'rest'
# opt_params_list = [4]
opt_params_list = [10]
# opt_params_list = [1,2,3,4,5,6]
# opt_params_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
# opt_params_list = [8,9,10,11,12,13,14]
# opt_params_list = [14]
# opt_params_list = [15,16,17,18,19,20,21]
# opt_params_list = [8,9,10,11,12,13,14,15,16,17,18,19,20,21]


#%% 
def main():
    for exp_ID in exp_IDs:
        for opt_params in opt_params_list:
            print('exp: {}, session_type: {}, opt_params: {}'.format(exp_ID,session_type,opt_params))
            path = Path(dir_IN) / session_type / exp_ID/ f'opt_{opt_params}'
            res = join_jobs_results(path)
            file_OUT = Path(dir_OUT) / session_type / exp_ID/ f'{exp_ID}_{session_type}_opt_{opt_params}.nc'
            file_OUT.parent.mkdir(parents=True,exist_ok=True)
            res.attrs
            res.to_netcdf(file_OUT)

#%%
if __name__ == "__main__":
    main()
