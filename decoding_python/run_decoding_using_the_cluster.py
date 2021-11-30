#%% imports
import pandas as pd
import numpy as np
import xarray as xr
import os
import sys
import exp_decode
import time

#%% main function
if __name__ == "__main__":
    exp_ID = 'b9861_d180527'
    job_index = int(os.environ['LSB_JOBINDEX'])
    print(type(job_index))
    # t = pd.read_csv('exp_table.csv')
    # t = t.set_index('job_index')
    # print(t)
    # exp_ID = t.loc[job_index].exp_ID
    print('Running decoding on the cluster!!!')
    print('exp ID:', exp_ID)
    print('job number:', job_index)

    exp = exp_decode.exp_decode(exp_ID)
    exp.load_data(folder='..//data_pre_proc')
    # exp.replay_fit(place_bin_size=0.5,pos_BW=1,replay_speed=80)
    # exp.replay_fit(place_bin_size=0.5,max_mark_value=1000,ds_mark_values=10)
    # exp.replay_fit(place_bin_size=1,use_marks_as_integers=True,max_mark_value=500,ds_mark_values=10)
    exp.replay_fit(place_bin_size=1,use_marks_as_integers=False)
    
    # ti = exp.PE_ti[440,:]
    # ti = ti + np.array([-1,1])*40e6
    
    # divide the data into 100 overlapping chunks (1 per job)
    # ti = 
    
    
    start = time.time()
    exp.sleep_predict(ti=ti)
    # exp.sleep_predict()
    end = time.time()
    print('analysis time:',end-start,'s')
    exp.save_sleep_res(folder='..//data_decoded')
    