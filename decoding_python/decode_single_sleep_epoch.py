#%% imports
import argparse
import pandas as pd
import numpy as np
import xarray as xr
import os
import sys
import exp_decode
import time
from pathlib import Path
import numba


#%% main function
def main():
    
    print('hello!')
    
    parser = argparse.ArgumentParser(description='decode during a specific sleep epoch (chunk)')
    parser.add_argument("--exp_ID", dest='exp_ID', type=str, default='b9861_d180527', help="exp_ID")
    # parser.add_argument("--job_index", dest='job_index', type=int, default=1, help="flight index (starting from 1)")
    parser.add_argument('--pos_bin_size', dest='pos_bin_size', type=float, default=1.0, help='position bin size (in meters)')
    parser.add_argument('--pos_std', dest='pos_std', type=float, default=2.0, help='position std (in meters)')
    parser.add_argument('--mark_std', dest='mark_std', type=float, default=20.0, help='marks std (in uVolts)')
    parser.add_argument('--state_decay_timescale', dest='state_decay_timescale', type=float, default=0.1, help='state_decay_timescale (in seconds)')
    parser.add_argument('--replay_speed', dest='replay_speed', type=int, default=20, help='replay speed')
    parser.add_argument('--res_dir', dest='res_dir', type=str, default='.', help='results folder')
    parser.add_argument("--num_epochs", dest='num_epochs', type=int, default=100, help="to how many epochs should we divide the entire sleep time")
    args = parser.parse_args()
    args.job_id = os.getenv('LSB_JOBID', 'None')
    args.job_index = int(os.getenv('LSB_JOBINDEX', 1))
    args.LSB_MAX_NUM_PROCESSORS = int(os.getenv('LSB_MAX_NUM_PROCESSORS', 1))
    args.LSB_DJOB_NUMPROC = int(os.getenv('LSB_DJOB_NUMPROC', 1))
    args.OMP_NUM_THREADS = int(os.getenv('OMP_NUM_THREADS', 1))
    args.epoch_idx = args.job_index-1
    for arg in vars(args):
        print(arg, ':', getattr(args, arg))

    print('numba num threads: ', numba.get_num_threads())
    
    exp = exp_decode.exp_decode(args.exp_ID)
    exp.load_data(folder='..//data_pre_proc')
    start = time.time()
    exp.replay_fit(state_decay_timescale = args.state_decay_timescale,
                   pos_bin_size = args.pos_bin_size,
                   replay_speed = args.replay_speed,
                   pos_std = args.pos_std,
                   mark_std = args.mark_std,
                   use_marks_as_integers = False)
    exp.sleep_predict_epoch(args.num_epochs, args.epoch_idx)
    end = time.time()
    print('analysis time:',end-start,'s')
    print('saving results')
    exp.sleep_results.attrs.update({'exp_ID' : args.exp_ID,
                                   'epoch_idx' : args.epoch_idx,
                                   'num_epochs' : args.num_epochs,
                                   'pos_bin_size' : args.pos_bin_size,
                                   'pos_std' : args.pos_std,
                                   'mark_std' : args.mark_std,
                                   'state_decay_timescale' : args.state_decay_timescale,
                                   'replay_speed' : args.replay_speed,
                                   'job_id' : args.job_id})
    
    folder = Path(args.res_dir)
    folder.mkdir(parents=True, exist_ok=True)
    exp.save_sleep_epoch_res(folder=folder)
    print('results saved')
   
#%%
if __name__ == "__main__":
    main()

    



#%%











