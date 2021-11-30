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
    Namespace(pos_bin_size = 0.5, pos_std = 1.0, mark_std = 15.0, state_decay_timescale = 2.0),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 15.0, state_decay_timescale = 2.0),
    Namespace(pos_bin_size = 0.5, pos_std = 1.0, mark_std = 20.0, state_decay_timescale = 2.0),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 2.0),
    Namespace(pos_bin_size = 0.5, pos_std = 1.0, mark_std = 25.0, state_decay_timescale = 2.0),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 25.0, state_decay_timescale = 2.0),
    Namespace(pos_bin_size = 0.2, pos_std = 0.4, mark_std = 15.0, state_decay_timescale = 2.0),
    Namespace(pos_bin_size = 0.2, pos_std = 0.4, mark_std = 20.0, state_decay_timescale = 2.0),
    Namespace(pos_bin_size = 0.2, pos_std = 0.4, mark_std = 25.0, state_decay_timescale = 2.0),
    ]

#%%
def get_flights_to_decode(exp_ID,opt_params):
    # get number of flights
    exp = exp_decode.exp_decode(exp_ID)
    exp.load_data(folder='..//data_pre_proc')
    n_FE = exp.FE_ti.shape[0]
    # look for existing decoded flight files
    folder = Path('..') / 'decoding' /'flight' / exp_ID / 'opt_{}'.format(opt_params)
    folder.mkdir(parents=True, exist_ok=True)
    file_template = exp_ID+"_flight_*_res.nc"
    existing_files = sorted(folder.glob(file_template))
    # print('existing decoded files:')
    # [print(file) for file in existing_files]
    decoded_flights = [int(re.findall(r'\d+', file.stem)[-1]) for file in existing_files]
    print('existing decoded flights:', decoded_flights)
    flights_to_decode = np.arange(n_FE)+1
    flights_to_decode = [x for x in flights_to_decode if x not in decoded_flights]
    return flights_to_decode, folder

#%%
def create_job_submission_str(exp_ID, opt_params, flights_to_decode, folder):
    # build job submission cmd
    python_function = 'decode_single_flight_epoch.py'
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
        f' --res_dir {res_dir}'
        )
    # bsub options
    queue_name = 'new-short'
    memory_usage = 1024
    num_cores = 24
    wall_time = '06:00'
    # job_slot_limit = 30
    job_list_str = str(flights_to_decode).replace(' ','')
    bsub_opts = (
        f' -q {queue_name}'
        f' -J jobs{job_list_str}'
        # f' -J jobs{job_list_str}%{job_slot_limit}'
        f' -R rusage[mem={memory_usage}]'
        f' -n {num_cores} -R "span[hosts=1]"'
        f' -W {wall_time}'
        f' -outdir "../jobs_output/%J"'
        f' -o "../jobs_output/%J/%J_%I.txt"'
        f' -env "all, NUMBA_NUM_THREADS={num_cores}"'
        # f' -env "NUMBA_NUM_THREADS={num_cores}"'
        # f' -env "OMP_NUM_THREADS={num_cores}"'
        # f' -env "MKL_NUM_THREADS={num_cores}"'
        )
    cmd_line_script = f'bsub {bsub_opts} python {python_cmd}'
    return cmd_line_script 

#%%
def submit_job(exp_ID,opt_params):

    flights_to_decode, folder = get_flights_to_decode(exp_ID,opt_params)
    if len(flights_to_decode)==0:
        print('no flights to decode! all flights were already decoded... no job submission!')
        return
    print('flights numbers to decode:', flights_to_decode)
    
    cmd_line_script = create_job_submission_str(exp_ID, opt_params, flights_to_decode, folder)
    print('job submission ccommand string:')
    print(cmd_line_script)
    run(cmd_line_script, shell=True)
    print('job submitted')
    return 

#%% main function
def main():
    parser = argparse.ArgumentParser(description='decode during a specific flight')
    parser.add_argument("--exp_ID", dest='exp_ID', type=str, nargs="+", default=['b9861_d180527'], help="exp_ID")
    parser.add_argument("--opt_params", dest='opt_params', type=int, nargs="+", default=[1], help="Decoding paramset opt(s)")
    # parser.add_argument("--queue_name", dest='queue_name', type=str, default="new-short", help="which queue to use")
    args = parser.parse_args()
    for arg in vars(args):
        print(arg, ':', getattr(args, arg))
    
    for exp_ID in args.exp_ID:
        for opt_params in args.opt_params:
            print('exp_ID:',exp_ID, 'opt_params:',opt_params)
            submit_job(exp_ID,opt_params)
    
#%%
if __name__ == "__main__":
    main()

    



#%%











