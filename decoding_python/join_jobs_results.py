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
    # # 'b9861_d180526',
    # # 'b9861_d180527',
    # 'b9861_d180529',
    # 'b9861_d180530',
    # 'b9861_d180531',
    # 'b9861_d180601',
    # 'b9861_d180603',
    # 'b9861_d180604',
    # 'b9861_d180605',
    # 'b9861_d180606',
    # 'b9861_d180607',
    # 'b9861_d180608',
    # 'b9861_d180609',
    # 'b9861_d180610',
    # 'b9861_d180611',
    # 'b9861_d180612',
    # 'b9861_d180613',
    # 'b9861_d180614',
    # 'b9861_d180615',
    # 'b9861_d180616',
    # 'b9861_d180617',
    # 'b9861_d180618',
    # 'b9861_d180619',
    # 'b9861_d180620',
    # 'b9861_d180621',
    # 'b9861_d180622',
    # 'b9861_d180623',
    # 'b9861_d180624',
    # 'b9861_d180625',
    # 'b9861_d180626',
    # 'b9861_d180627',
    # 'b9861_d180628',
    # 'b9861_d180629',
    # 'b9861_d180630',
    # 'b9861_d180701',
    # 'b9861_d180702',
    # 'b9861_d180703',
    # 'b9861_d180704',
    # 'b9861_d180705',
    # 'b9861_d180706',
    # 'b9861_d180708',
    # 'b9861_d180709',
    # 'b9861_d180710',
    # 'b9861_d180711',
    
    # 'b0184_d191128',
    # 'b0184_d191129',
    # 'b0184_d191130',
    # 'b0184_d191201',
    # 'b0184_d191202',
    # 'b0184_d191203',
    # 'b0184_d191204',
    # 'b0184_d191205',
    # 'b0184_d191208',
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
    # 'b2289_d180523',
    # 'b2289_d180524',
    # 'b2289_d180525',
    # 'b2289_d180528',
    # 'b2289_d180529',
    # 'b2289_d180531',
    # 'b2289_d180619',
    # 'b2289_d180620',
    
    # 'b0079_d160908',
    # 'b0079_d160909',
    # 'b0079_d160911',
    # 'b0079_d160912',
    # 'b0079_d160913',
    # 'b0079_d160914',
    # 'b0079_d160915',
    # 'b0079_d160916',
    # 'b0079_d160918',
    # 'b0079_d160919',
    # 'b0079_d160920',
    # 'b0079_d160921',
    # 'b0079_d160923',
    # 'b0079_d160924',
    # 'b0079_d160925',
    # 'b0079_d160926',
    # 'b0079_d160927',
    # 'b0079_d160928',
    # 'b0079_d160929',
    # 'b0079_d160930',
    # 'b0079_d161003',
    # 'b0079_d161004',
    # 'b0079_d161005',
        
    # 'b0148_d170606',
    # 'b0148_d170607',
    # 'b0148_d170608',
    # 'b0148_d170611',
    # 'b0148_d170612',
    # 'b0148_d170613',
    # 'b0148_d170614',
    # 'b0148_d170615',
    # 'b0148_d170618',
    # 'b0148_d170619',
    # 'b0148_d170620',
    # 'b0148_d170621',
    # 'b0148_d170622',
    # 'b0148_d170625',
    # 'b0148_d170626',
    # 'b0148_d170627',
    # 'b0148_d170628',
    # 'b0148_d170703',
    # 'b0148_d170704',
    # 'b0148_d170705',
    # 'b0148_d170710',
    # 'b0148_d170711',
    # 'b0148_d170712',
    # 'b0148_d170713',
    # 'b0148_d170716',
    # 'b0148_d170717',
    # 'b0148_d170718',
    # 'b0148_d170720',
    # 'b0148_d170723',
    # 'b0148_d170801',
    # 'b0148_d170802',
    # 'b0148_d170803',
    # 'b0148_d170806',
    # 'b0148_d170807',
    
    ]
session_type = 'flight'
# session_type = 'sleep'
# session_type = 'rest'
opt_params_list = [4]
# opt_params_list = [10]
# opt_params_list = [1,2,3,4,5,6]
# opt_params_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
# opt_params_list = [8,9,10,11,12,13,14]
# opt_params_list = [11]
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
