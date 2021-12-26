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
    
    parser = argparse.ArgumentParser(description='decode during a specific flight')
    parser.add_argument("--exp_ID", dest='exp_ID', type=str, default='b9861_d180527', help="exp_ID")
    # parser.add_argument("--job_index", dest='job_index', type=int, default=1, help="flight index (starting from 1)")
    parser.add_argument('--pos_bin_size', dest='pos_bin_size', type=float, default=1.0, help='position bin size (in meters)')
    parser.add_argument('--pos_std', dest='pos_std', type=float, default=2.0, help='position std (in meters)')
    parser.add_argument('--mark_std', dest='mark_std', type=float, default=20.0, help='marks std (in uVolts)')
    parser.add_argument('--state_decay_timescale', dest='state_decay_timescale', type=float, default=1.0, help='state_decay_timescale (in seconds)')
    parser.add_argument('--replay_speed', dest='replay_speed', type=int, default=1, help='replay speed')
    parser.add_argument('--res_dir', dest='res_dir', type=str, default='.', help='results folder')
    args = parser.parse_args()
    args.job_id = os.getenv('LSB_JOBID', 'None')
    args.job_index = int(os.getenv('LSB_JOBINDEX', 'None'))
    args.LSB_MAX_NUM_PROCESSORS = int(os.getenv('LSB_MAX_NUM_PROCESSORS', 1))
    args.LSB_DJOB_NUMPROC = int(os.getenv('LSB_DJOB_NUMPROC', 1))
    args.flight_idx = args.job_index-1
    for arg in vars(args):
        print(arg, ':', getattr(args, arg))

    print('numba num threads: ', numba.get_num_threads())
    
    exp = exp_decode.exp_decode(args.exp_ID)
    exp.load_data(folder='..//data_pre_proc')
    start = time.time()
    exp.decode_behave(replay_speed=args.replay_speed,
                      state_decay_timescale = args.state_decay_timescale,
                      pos_bin_size = args.pos_bin_size,
                      pos_std = args.pos_std,
                      mark_std = args.mark_std,
                      flight_idx = args.flight_idx)
    end = time.time()
    print('analysis time:',end-start,'s')
    print('saving results')
    exp.behave_results.attrs.update({'exp_ID' : args.exp_ID,
                                     'flight_idx' : args.flight_idx,
                                     'pos_bin_size' : args.pos_bin_size,
                                     'pos_std' : args.pos_std,
                                     'mark_std' : args.mark_std,
                                     'state_decay_timescale' : args.state_decay_timescale,
                                     'replay_speed' : args.replay_speed,
                                     'job_id' : args.job_id})
    # folder = Path('.').resolve().parent.joinpath('data_decoded',args.job_id)
    folder = Path(args.res_dir)
    folder.mkdir(parents=True, exist_ok=True)
    exp.save_flight_res(folder=folder)
    print('results saved')
   
#%%
if __name__ == "__main__":
    main()

    



#%%











