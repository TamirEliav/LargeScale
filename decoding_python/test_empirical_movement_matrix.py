# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 14:16:54 2021

@author: tamire
"""

#%%
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sb
import exp_decode
from pathlib import Path
from scipy.io import savemat

#%%
exp_ID='b0184_d191201'
# exp_ID='b9861_d180525'
exp=exp_decode.exp_decode(exp_ID)
exp.load_data(folder='..//data_pre_proc')
dir_out = 'F:\\sequences\\results\\test_empirical_movement_matrix'
path = Path(dir_out)

#%%
#for replay_speed in [1,10,20,40,80,160,200]:
for replay_speed in [40]:
    exp.replay_fit(pos_bin_size=1.0,replay_speed=replay_speed )
    M = exp.replay_classifier.continuous_state_transition_
    filename = f'{exp_ID}_continuous_state_transition_replay_speed_{replay_speed}.mat'
    fileout = path / filename
    savemat(fileout,{'M':M})
    # np.save(fileout, M)
    
    # fig, ax = plt.subplots(figsize=(11, 9))
    # x = exp.replay_classifier.place_bin_centers_
    # y = exp.replay_classifier.place_bin_centers_
    # M = exp.replay_classifier.continuous_state_transition_[0,0,:,:]
    # sb.heatmap(M, annot=True, fmt=".2f", cmap='bone_r')
    
    # ax.xlim = [0,150]
    # ax.ylim = [0,150]
    # ax.set_title(f'replay speed: {replay_speed}')

    # filename = f'{exp_ID}_replay_speed_{replay_speed}.jpg'
    # fileout = path / filename
    # fig.savefig(fileout)
