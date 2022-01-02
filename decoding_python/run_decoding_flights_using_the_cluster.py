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
    Namespace(pos_bin_size = 0.5, pos_std = 1.0, mark_std = 15.0, state_decay_timescale = 2.0, replay_speed=1),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 15.0, state_decay_timescale = 2.0, replay_speed=1),
    Namespace(pos_bin_size = 0.5, pos_std = 1.0, mark_std = 20.0, state_decay_timescale = 2.0, replay_speed=1),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 20.0, state_decay_timescale = 2.0, replay_speed=1),
    Namespace(pos_bin_size = 0.5, pos_std = 1.0, mark_std = 25.0, state_decay_timescale = 2.0, replay_speed=1),
    Namespace(pos_bin_size = 1.0, pos_std = 2.0, mark_std = 25.0, state_decay_timescale = 2.0, replay_speed=1),
    Namespace(pos_bin_size = 0.2, pos_std = 0.4, mark_std = 15.0, state_decay_timescale = 2.0, replay_speed=1),
    Namespace(pos_bin_size = 0.2, pos_std = 0.4, mark_std = 20.0, state_decay_timescale = 2.0, replay_speed=1),
    Namespace(pos_bin_size = 0.2, pos_std = 0.4, mark_std = 25.0, state_decay_timescale = 2.0, replay_speed=1),
    
    Namespace(pos_bin_size = 0.2, pos_std = 0.4, mark_std = 20.0, state_decay_timescale = 4.0, replay_speed=40),
    ]

#%% 
exp_ID_list = [
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
    # 'b0184_d191222',
    # 'b0184_d191224',
    # 'b0184_d191225',
    # 'b0184_d191226',
    # 'b0184_d200101',
    # 'b0184_d200102',

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

    'b0148_d170606',
    'b0148_d170607',
    'b0148_d170608',
    'b0148_d170611',
    'b0148_d170612',
    'b0148_d170613',
    'b0148_d170614',
    'b0148_d170615',
    'b0148_d170618',
    'b0148_d170619',
    'b0148_d170620',
    'b0148_d170621',
    'b0148_d170622',
    'b0148_d170625',
    'b0148_d170626',
    'b0148_d170627',
    'b0148_d170628',
    'b0148_d170703',
    'b0148_d170704',
    'b0148_d170705',
    'b0148_d170710',
    'b0148_d170711',
    'b0148_d170712',
    'b0148_d170713',
    'b0148_d170716',
    'b0148_d170717',
    'b0148_d170718',
    'b0148_d170720',
    'b0148_d170723',
    'b0148_d170801',
    'b0148_d170802',
    'b0148_d170803',
    'b0148_d170806',
    'b0148_d170807',

    ]

opt_params_list = [4]

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
        f' --replay_speed {params.replay_speed}'
        f' --res_dir {res_dir}'
        )
    # bsub options
    queue_name = 'new-short'
    memory_usage = 8192
    NUMBA_NUM_THREADS = 16
    wall_time = '02:00'
    job_list_str = str(flights_to_decode).replace(' ','')
    bsub_opts = (
        f' -q {queue_name}'
        f' -J flight_{exp_ID}{job_list_str}'
        # f' -J jobs{job_list_str}%{job_slot_limit}'
        f' -R rusage[mem={memory_usage}]'
        # f' -n {num_cores} -R "span[hosts=1]"'
        # f' -R "affinity[core({num_cores})]"'
        f' -R "affinity[thread({NUMBA_NUM_THREADS})]"'
        # f' -R "affinity[thread*{NUMBA_NUM_THREADS}]"'
        f' -We {wall_time}'
        f' -outdir "../jobs_output/%J"'
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
    parser.add_argument("--use_lists_in_code", dest='use_lists_in_code', type=int, default=[0], help="Use exp/parrams list from code (0/1)")
    # parser.add_argument("--queue_name", dest='queue_name', type=str, default="new-short", help="which queue to use")
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











