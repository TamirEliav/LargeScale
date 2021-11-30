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
import re
from subprocess import run
from argparse import Namespace

#%%
paramsets = [
    Namespace(pos_bin_size = 0.5, pos_std = 1.0, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=1),     # 1
    Namespace(pos_bin_size = 0.5, pos_std = 1.0, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=10),
    Namespace(pos_bin_size = 0.5, pos_std = 1.0, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=20),
    Namespace(pos_bin_size = 0.5, pos_std = 1.0, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=40),
    Namespace(pos_bin_size = 0.5, pos_std = 1.0, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=80),
    Namespace(pos_bin_size = 0.5, pos_std = 1.0, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=160),
    Namespace(pos_bin_size = 0.5, pos_std = 1.0, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=200),
    
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=1),     # 8
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=10),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=20),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=40),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=80),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=160),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=200),

    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.05, replay_speed=1),    # 15
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.05, replay_speed=10),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.05, replay_speed=20),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.05, replay_speed=40),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.05, replay_speed=80),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.05, replay_speed=160),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.05, replay_speed=200),
    ]

#%%
def get_epochs_to_decode(exp_ID,opt_params,num_epochs):
    # get number of sleep  epochs
    exp = exp_decode.exp_decode(exp_ID)
    exp.load_data(folder='..//data_pre_proc')
    # look for existing decoded sleep files
    folder = Path('..') / 'decoding' / 'sleep' / exp_ID / 'opt_{}'.format(opt_params)
    folder.mkdir(parents=True, exist_ok=True)
    file_template = exp_ID+"_sleep_*_res.nc"
    existing_files = sorted(folder.glob(file_template))
    # print('existing decoded files:')
    # [print(file) for file in existing_files]
    decoded_epochs= [int(re.findall(r'\d+', file.stem)[-1]) for file in existing_files]
    print('existing decoded sleep epochs:', decoded_epochs)
    epochs_to_decode = np.arange(num_epochs)+1
    epochs_to_decode = [x for x in epochs_to_decode if x not in decoded_epochs]
    return epochs_to_decode, folder

#%%
def create_job_submission_str(exp_ID, opt_params, num_epochs, epochs_to_decode, folder):
    # build job submission cmd
    python_function = 'decode_single_sleep_epoch.py'
    params = paramsets[opt_params-1]
    [print(x,getattr(params,x)) for x in vars(params)]
    res_dir = folder.as_posix()
    python_cmd = (
        f'{python_function}'
        f' --exp_ID {exp_ID}'
        # f' --job_index ${{LSB_JOBINDEX}}'
        f' --pos_bin_size {params.pos_bin_size}'
        f' --pos_std {params.pos_std}'
        f' --mark_std {params.mark_std}'
        f' --state_decay_timescale {params.state_decay_timescale}'
        f' --replay_speed {params.replay_speed}'
        f' --res_dir {res_dir}'
        f' --num_epochs {num_epochs}'
        )
    # bsub options
    queue_name = 'new-short'
    memory_usage = 4096
    # num_cores = 1
    NUMBA_NUM_THREADS = 8
    wall_time = '01:00'
    # job_slot_limit = 30
    job_list_str = str(epochs_to_decode).replace(' ','')
    bsub_opts = (
        f' -q {queue_name}'
        f' -J jobs{job_list_str}'
        # f' -J jobs{job_list_str}%{job_slot_limit}'
        f' -R rusage[mem={memory_usage}]'
        # f' -n {num_cores} -R "span[hosts=1]"'
        # f' -R "affinity[core({num_cores})]"'
        f' -R "affinity[thread({NUMBA_NUM_THREADS})]"'
        # f' -R "affinity[thread*{NUMBA_NUM_THREADS}]"'
        f' -We {wall_time}'
        f' -outdir "../jobs_output/%J/"'
        f' -o "../jobs_output/%J/%J_%I.txt"'
        f' -env "all, NUMBA_NUM_THREADS={NUMBA_NUM_THREADS}"'
        # f' -env "NUMBA_NUM_THREADS={num_cores}"'
        # f' -env "OMP_NUM_THREADS={num_cores}"'
        # f' -env "MKL_NUM_THREADS={num_cores}"'
        )
    cmd_line_script = f'bsub {bsub_opts} python {python_cmd}'
    return cmd_line_script 

#%%
def submit_job(exp_ID,opt_params,num_epochs):

    epochs_to_decode, folder = get_epochs_to_decode(exp_ID,opt_params,num_epochs)
    if len(epochs_to_decode)==0:
        print('no sleep epochs to decode! all epochs were already decoded... no job submission!')
        return
    print('sleep epochs numbers to decode:', epochs_to_decode)
    
    cmd_line_script = create_job_submission_str(exp_ID, opt_params, num_epochs, epochs_to_decode, folder)
    print('job submission ccommand string:')
    print(cmd_line_script)
    run(cmd_line_script, shell=True)
    print('job submitted')
    return 

#%% main function
def main():
    parser = argparse.ArgumentParser(description='decode during a specific sleep epoch (chunk)')
    parser.add_argument("--exp_ID", dest='exp_ID', type=str, nargs="+", default=['b9861_d180527'], help="exp_ID")
    parser.add_argument("--opt_params", dest='opt_params', type=int, nargs="+", default=[1], help="Decoding paramset opt(s)")
    parser.add_argument("--num_epochs", dest='num_epochs', type=int, default=100, help="to how many epochs should we divide the entire sleep time")
    args = parser.parse_args()
    for arg in vars(args):
        print(arg, ':', getattr(args, arg))
    
    for exp_ID in args.exp_ID:
        for opt_params in args.opt_params:
            print('exp_ID:',exp_ID, 'opt_params:',opt_params)
            submit_job(exp_ID, opt_params, args.num_epochs)
    
#%%
if __name__ == "__main__":
    main()

    



#%%











