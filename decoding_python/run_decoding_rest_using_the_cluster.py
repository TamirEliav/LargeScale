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


    Namespace(pos_bin_size = 0.2, pos_std = 0.4, mark_std = 20.0, state_decay_timescale = 0.1, replay_speed=40),   # 15
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 0.5, replay_speed=40),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 1.0, replay_speed=40),
    Namespace(pos_bin_size = 0.3, pos_std = 0.3, mark_std = 20.0, state_decay_timescale = 1.0, replay_speed=40),
    Namespace(pos_bin_size = 0.4, pos_std = 0.4, mark_std = 20.0, state_decay_timescale = 1.0, replay_speed=40),
    Namespace(pos_bin_size = 0.5, pos_std = 0.5, mark_std = 20.0, state_decay_timescale = 1.0, replay_speed=40), # 20
    ]


#%% 
exp_ID_list = [
    # 'b0034_d180312',
    # 'b0034_d180313',
    # 'b0148_d170606',
    # 'b0148_d170613',
    # 'b0148_d170614',
    # 'b0148_d170619',
    # 'b0148_d170620',
    # 'b0148_d170621',
    # 'b0148_d170622',
    # 'b0148_d170625',
    # 'b0148_d170626',
    # 'b0148_d170627',
    # 'b0148_d170628',
    # 'b0148_d170703',
    # 'b0148_d170718',
    # 'b0148_d170720',
    # 'b0148_d170802',
    # 'b2289_d180518',
    # 'b9861_d180524',
    # 'b9861_d180525',
    # 'b9861_d180526',
    # 'b9861_d180527',
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
    # 'b0184_d191220',
    # 'b0184_d191225',
    # 'b0184_d200102',
    
    # 'b2382_d190620',
    # 'b2382_d190623',
    # 'b2382_d190624',
    # 'b2382_d190627',
    # 'b2382_d190703',
    # 'b2382_d190709',
    # 'b2382_d190712',
    # 'b2382_d190714',
    # 'b2382_d190715',
    # 'b2382_d190716',
    # 'b2382_d190718',
    # 'b2382_d190721',
    # 'b2382_d190722',
    # 'b2382_d190725',
    # 'b2382_d190728',
    # 'b2382_d190729',
    # 'b2382_d190730',
    # 'b2382_d190731',
    # 'b2382_d190801',
    # 'b2382_d190804',
    # 'b2382_d190805',
    # 'b2382_d190807',
    # 'b2382_d190808',
    # 'b2382_d190811',
    # 'b2382_d190812',
    # 'b2382_d190813',
    # 'b2382_d190814', # this has a very long rest epoch (number 4), need to run seperately with more memory/cores/CPU-time
    # 'b0194_d180503',
    # 'b0194_d180505',
    # 'b0194_d180507',
    # 'b0194_d180508',
    # 'b0194_d180509',
    # 'b0194_d180510',
    # 'b0194_d180513',
    # 'b0194_d180514',
    # 'b0194_d180515',
    # 'b0194_d180516',
    # 'b0194_d180520',
    # 'b0194_d180521',
    # 'b0194_d180604',
    # 'b0194_d180605',
    # 'b0194_d180606',
    # 'b0194_d180612',
    # 'b0194_d180614',    
    
    # 'b2299_d191124',
    # 'b2299_d191125',
    # 'b2299_d191126',
    # 'b2299_d191127',
    # 'b2299_d191128',
    # 'b2299_d191201',
    # 'b2299_d191202',
    # 'b2299_d191203',
    # 'b2299_d191204',
    # 'b2299_d191205',
    # 'b2299_d191208',
    # 'b2299_d191209',
    # 'b2299_d191210',
    # 'b2299_d191211',
    # 'b2299_d191212',
    # 'b2299_d191213',
    
    # 6m
    'b0184_d191124',
    'b0184_d191125',
    'b0184_d191126',
    'b0184_d191127',
    
    ]

# opt_params_list = [8,9,10,11,12,13,14]
# opt_params_list = [11]
# opt_params_list = [15]
opt_params_list = [19,20]
# opt_params_list = [16,17]
# opt_params_list = [25]
# opt_params_list = [8,9,10,12,13,14]

#%%
def get_epochs_to_decode(exp_ID,opt_params):
    # get number of rest epochs
    exp = exp_decode.exp_decode(exp_ID)
    exp.load_data(folder='..//data_pre_proc')
    n_epochs = exp.rest_ti.shape[0]
    # look for existing decoded rest files
    folder = Path('..') / 'decoding' / 'rest' / exp_ID / 'opt_{}'.format(opt_params)
    folder.mkdir(parents=True, exist_ok=True)
    file_template = exp_ID+"_rest_*_res.nc"
    existing_files = sorted(folder.glob(file_template))
    # print('existing decoded files:')
    # [print(file) for file in existing_files]
    decoded_epochs= [int(re.findall(r'\d+', file.stem)[-1]) for file in existing_files]
    print('existing decoded rest epochs:', decoded_epochs)
    epochs_to_decode = np.arange(n_epochs)+1
    epochs_to_decode = [x for x in epochs_to_decode if x not in decoded_epochs]
    return epochs_to_decode, folder

#%%
def create_job_submission_str(exp_ID, opt_params, epochs_to_decode, folder):
    # build job submission cmd
    python_function = 'decode_single_rest_epoch.py'
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
        )
    # bsub options
    queue_name = 'new-short'
    # queue_name = 'new-medium'
    # memory_usage = 4096
    memory_usage = 16384
    # memory_usage = 65536
    # memory_usage = 8192
    # num_cores = 1
    # NUMBA_NUM_THREADS = 8
    NUMBA_NUM_THREADS = 16
    # NUMBA_NUM_THREADS = 32
    wall_time = '01:00'
    # job_slot_limit = 30
    job_list_str = str(epochs_to_decode).replace(' ','')
    bsub_opts = (
        f' -q {queue_name}'
        f' -J rest_{exp_ID}{job_list_str}'
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
def submit_job(exp_ID,opt_params):

    epochs_to_decode, folder = get_epochs_to_decode(exp_ID,opt_params)
    if len(epochs_to_decode)==0:
        print('no rest epochs to decode! all epochs were already decoded... no job submission!')
        return
    print('rest epochs numbers to decode:', epochs_to_decode)
    
    cmd_line_script = create_job_submission_str(exp_ID, opt_params, epochs_to_decode, folder)
    print('job submission ccommand string:')
    print(cmd_line_script)
    run(cmd_line_script, shell=True)
    print('job submitted')
    return 

#%% main function
def main():
    parser = argparse.ArgumentParser(description='decode during a specific rest epoch')
    parser.add_argument("-e", "--exp_ID", dest='exp_ID', type=str, nargs="+", default=['b9861_d180527'], help="exp_ID")
    parser.add_argument("-o", "--opt_params", dest='opt_params', type=int, nargs="+", default=[1], help="Decoding paramset opt(s)")
    parser.add_argument("--use_lists_in_code", dest='use_lists_in_code', type=int, default=0, help="Use exp/parrams list from code (0/1)")
    args = parser.parse_args()
    if args.use_lists_in_code:
        args.exp_ID = exp_ID_list
        args.opt_params = opt_params_list
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











