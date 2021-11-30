#%% imports
import argparse
from argparse import ArgumentParser
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
    parser.add_argument("-i", "--job_index", dest='job_index', type=int, default=1, help="flight index (starting from 1)")
    
    args = parser.parse_args()
    job_index = args.job_index
    job_id = os.getenv('LSB_JOBID', 'None')
    
    print('job id:', job_id)
    print('job number:', job_index)
    
    # return 
    print('num cores:', numba.config.NUMBA_NUM_THREADS)
    
    
#%%
if __name__ == "__main__":
    main()

    
#%%
python_function = 'decode_single_flight.py'
exp_ID = 'b9861_d180527'
epoch = 'flight'
python_cmd = (
    f'{python_function} {animal} {day} {epoch}'
    f' --data_type {args.data_type}'
    f' --dim {args.dim}'
    f' --n_workers {args.n_workers}'
    f' --threads_per_worker {args.threads_per_worker}')
directives = ' '.join(
        [f'-l h_rt={args.wall_time}', f'-pe omp {args.n_cores}',
         '-P braincom', '-notify',
         '-v OPENBLAS_NUM_THREADS', '-v NUMBA_NUM_THREADS',
         '-v OMP_NUM_THREADS'])

queue_cmd = f'qsub {directives} -j y -o {log_file} -N {job_name}'
cmd_line_script = f'echo python {python_cmd}  | {queue_cmd}'
# run(cmd_line_script, shell=True)


#%%











