# -*- coding: utf-8 -*-
"""
Created on Sun Jun 27 08:38:02 2021

@author: tamire
"""

#%% import
import logging
import os
# import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd
import xarray as xr
from matplotlib.colors import LogNorm
from sklearn.metrics import confusion_matrix
import pickle
import joblib

from replay_trajectory_classification import ClusterlessClassifier
from replay_trajectory_classification.state_transition import estimate_movement_var
from replay_trajectory_classification.misc import NumbaKDE

# %%
# exp = exp_decode.exp_decode('b9861_d180526')
# exp.load_data()
# exp.decode_sleep(win_s=0.05)

#%%
def is_in_ti(t,ti):
    t = t[:,np.newaxis]
    ti = np.atleast_2d(ti)
    return ((t>ti[:,0]) & (t < ti[:,1])).any(axis=1) 

#%%
class exp_decode():
    #%%
    def __init__(self,exp_ID="b0034_d180313"):
        self.exp_ID = exp_ID
       
    #%%
    def load_data(self,folder = "F:\\sequences\\proc\\"):
    
        #%% load data
        data_file_name = os.path.join(folder,self.exp_ID+".mat")
        mat = scipy.io.loadmat(data_file_name)
        
        # get struct fields
        fs = mat['fs'][0][0]
        time = mat['t'][0,:]
        time_2d = time[:,np.newaxis]
        multiunits = mat['multiunits']
        # workaround to solve matlab not allowing singeltons in last dim
        if multiunits.ndim == 2:
            multiunits = multiunits[:,:,None]
        multiunits_spikes = mat['multiunits_spikes'].astype(bool)
        position = mat['position'][0,:]
        direction = mat['direction'][0,:]
        FE_ti = mat['FE_ti']
        sleep_ti = mat['sleep_ti']
        # PE_ts = mat['PE_ts'][0,:]
        # PE_ti = mat['PE_ti']
        rest_ti = mat['rest_ti']
        
        #%% define data labels
        is_FE = ((time_2d>FE_ti[:,0]) & (time_2d<FE_ti[:,1])).any(axis=1)
        is_sleep = ((time_2d>sleep_ti[:,0]) & (time_2d<sleep_ti[:,1])).any(axis=1)
        encoding_labels = np.where(is_sleep, 'Sleep', 'None')
        encoding_labels = encoding_labels.astype('<U8')
        encoding_labels[is_FE & (direction==1)] = 'Outbound'
        encoding_labels[is_FE & (direction==-1)] = 'Inbound'
        
        # is_inbound = direction==-1
        # encoding_labels = np.where(is_inbound, 'Inbound', 'Outbound')
        
        #%% invalidate behavior outside flights
        position[~is_FE] = np.nan
        direction[~is_FE] = np.nan

        #%% create xarray dataset
        data = xr.Dataset({
            "multiunits":(["time","ch","TT"],multiunits),
            "multiunits_spikes":(["time","TT"],multiunits_spikes),
            "position":(["time"], position),
            "direction":(["time"], direction),
            "is_FE":(["time"], is_FE),
            "is_sleep":(["time"], is_sleep),
            "encoding_labels":(["time"], encoding_labels),
                        }, 
            coords={"time":time,
                    "TT":range(multiunits.shape[2]),
                    "ch":range(multiunits.shape[1]),
                       },
            attrs={"fs":fs,
                   })
        
        #%% assign everything to self
        self.data = data
        # self.PE_ts = PE_ts
        # self.PE_ti = PE_ti
        self.rest_ti = rest_ti
        self.FE_ti = FE_ti
        self.sleep_ti = sleep_ti
        
    #%%
    def replay_fit(self,  state_decay_timescale = 0.1,
                         pos_bin_size=0.5,
                         replay_speed=20,
                         pos_std=2.0,
                         mark_std=20.0,
                         use_marks_as_integers=False,
                         max_mark_value=500,
                         ds_mark_values=5):
        # Sleep decoding - train on data from all flights
        training_ti = self.FE_ti
        is_training = is_in_ti(self.data.time.values,training_ti)
        
        if use_marks_as_integers:
            # If your marks are integers, use this algorithm because it is much faster
            # the following line should be done inside the classifier code
            self.data['multiunits'] = xr.ufuncs.minimum(self.data.multiunits,max_mark_value)
            self.data['multiunits'] = np.round(self.data['multiunits'] / ds_mark_values)
            max_mark_value = int(self.data['multiunits'].max())
            clusterless_algorithm = 'multiunit_likelihood_integer'
            clusterless_algorithm_params = {
                'mark_std' : mark_std / ds_mark_values,
                'position_std' : pos_std,
                'dtype' : np.int32,
                'max_mark_value' : max_mark_value
            }
        else:
            ## If your marks are floats, use this algorithm
            clusterless_algorithm = 'multiunit_likelihood'
            clusterless_algorithm_params = {
                'model': NumbaKDE,
                'model_kwargs': {
                      'bandwidth': np.array([mark_std, mark_std, mark_std, mark_std, 
                                             pos_std]) # amplitude 1, amplitude 2, amplitude 3, amplitude 4, position       
                }
            }
        
        continuous_transition_types = [
            ['empirical_movement', 'uniform',  'uniform',   'uniform', 'uniform', 'uniform'],
            ['uniform',            'identity', 'uniform',   'uniform', 'uniform', 'uniform'],
            ['uniform',            'uniform',  'uniform',   'uniform', 'uniform', 'uniform'],
            ['uniform',            'uniform',  'uniform',   'empirical_movement', 'uniform', 'uniform'],
            ['uniform',            'uniform',  'uniform',   'uniform',            'identity', 'uniform'],
            ['uniform',            'uniform',  'uniform',   'uniform',            'uniform', 'uniform'],
            ]
        encoding_group_to_state = ['Inbound', 'Inbound', 'Inbound', 'Outbound', 'Outbound', 'Outbound']
        
        dt = 1/self.data.attrs['fs']
        state_diag_val = 1-1/(state_decay_timescale/dt)
        classifier = ClusterlessClassifier(replay_speed=replay_speed,
                                           place_bin_size=pos_bin_size,
                                           clusterless_algorithm=clusterless_algorithm,
                                           clusterless_algorithm_params=clusterless_algorithm_params,
                                           continuous_transition_types=continuous_transition_types,
                                           discrete_transition_diag=state_diag_val)
        classifier.fit(self.data.position,
                       self.data.multiunits,
                       encoding_group_labels=self.data.encoding_labels,
                       encoding_group_to_state=encoding_group_to_state,
                       is_training=is_training)
        self.replay_classifier = classifier
    
    #%%
    def replay_predict(self, ti):
        is_testing = is_in_ti(self.data.time.values,ti)
        results = self.replay_classifier.predict(
            self.data.multiunits.isel({"time":is_testing}),
            self.data.time.isel({"time":is_testing}))
        self.replay_results = results
    
    #%%
    def sleep_predict(self, win_s=None, ti=None, is_testing=None):
        # Sleep decoding - train on data from all flights
        data = self.data
        if win_s==None:
            # opt 1) all sleep times
            testing_ti = self.sleep_ti
        else:
            # opt 2) use win only around PE
            win = np.array([-1,1]) * win_s * 1e6
            # testing_ti = self.PE_ts[:,np.newaxis] + win[np.newaxis,:]
            testing_ti = self.PE_ti + win[np.newaxis,:]
        if ti is not None:
            testing_ti = ti
        if is_testing is not None:
            is_testing = is_in_ti(data.time.values,testing_ti)
        
        results = self.replay_classifier.predict(
            data.multiunits.isel({"time":is_testing}),
            data.time.isel({"time":is_testing}))
        self.sleep_results = results
    
    #%%
    def sleep_predict_epoch(self, num_epochs, epoch_idx):
        self.epoch_idx = epoch_idx
        
        data = self.data.isel({"time":self.data.is_sleep})
        L = data.time.shape[0]
        epoch = np.arange(L)//(np.ceil(L/num_epochs)).astype(int)
        data = data.assign({"epoch":(["time"],epoch)})
        
        IX1 = data.epoch == epoch_idx
        IX2 = (data.epoch == epoch_idx-1) | (data.epoch == epoch_idx+1)
        IX3 = IX1 | IX2
        IX4 = IX1.where(IX3,drop=True)
        
        results = self.replay_classifier.predict(
            data.multiunits.isel({"time":IX3}),
            data.time.isel({"time":IX3}))
        results = results.where(IX4, drop=True)
        self.sleep_results = results
    
    #%%
    def rest_predict(self, rest_idx=0):
        self.rest_idx = rest_idx
        is_testing = is_in_ti(self.data.time.values,self.rest_ti[rest_idx,:])
        
        results = self.replay_classifier.predict(
            self.data.multiunits.isel({"time":is_testing}),
            self.data.time.isel({"time":is_testing}))
        self.rest_results = results
    
    
    #%% 
    def save_sleep_epoch_res(self,folder = "D:\\sequences\\seq_uri_eden\\proc\\"):
        filename = os.path.join(folder,self.exp_ID+"_sleep_"+str(self.epoch_idx+1)+"_res.nc")
        results = self.sleep_results
        # results = results.drop_vars("likelihood")
        results = results.drop_vars("causal_posterior")
        results.to_netcdf(filename)
        
    #%% 
    def save_sleep_res(self,folder = "D:\\sequences\\seq_uri_eden\\proc\\"):
        filename = os.path.join(folder,self.exp_ID+"_sleep_res.nc")
        results = self.sleep_results
        # results = results.drop_vars("likelihood")
        results = results.drop_vars("causal_posterior")
        results.to_netcdf(filename)
        
    #%% 
    def save_rest_res(self,folder = "D:\\sequences\\seq_uri_eden\\proc\\"):
        filename = os.path.join(folder,self.exp_ID+"_rest_"+str(self.rest_idx+1)+"_res.nc")
        results = self.rest_results
        # results = results.drop_vars("likelihood")
        results = results.drop_vars("causal_posterior")
        results.to_netcdf(filename)

    #%% 
    def save_flight_res(self,folder = "D:\\sequences\\seq_uri_eden\\proc\\"):
        filename = os.path.join(folder,self.exp_ID+"_flight_"+str(self.flight_idx+1)+"_res.nc")
        results = self.behave_results
        # results = results.drop_vars("likelihood")
        results = results.drop_vars("causal_posterior")
        results.to_netcdf(filename)
        
    #%%
    def save_data(self,folder = "D:\\sequences\\seq_uri_eden\\proc\\"):
        filename = os.path.join(folder,self.exp_ID+"_data.pkl")
        joblib.dump(self.data, filename)    
    
    #%%
    def decode_behave(self, 
                      replay_speed = 1,
                      state_decay_timescale = 2,
                      pos_bin_size=0.5,
                      pos_std=2.0,
                      mark_std=20.0,
                      flight_idx=0,
                      win_s = None):
        self.flight_idx = flight_idx
        # Sleep decoding - train on data from all flights except one (leave one out)
        data = self.data
        FE_ti = self.FE_ti
        training_ti = FE_ti[flight_idx != np.arange(FE_ti.shape[0]),:]
        testing_ti = FE_ti[flight_idx,:]
        if win_s is not None:
            testing_ti = testing_ti.mean() + win_s*1e6*np.array([-1,1])
        is_training = is_in_ti(data.time.values,training_ti)
        is_testing = is_in_ti(data.time.values,testing_ti)
                
        clusterless_algorithm = 'multiunit_likelihood'
        clusterless_algorithm_params = {
            'model': NumbaKDE,
            'model_kwargs': {
                  'bandwidth': np.array([mark_std, mark_std, mark_std, mark_std, 
                                         pos_std]) # amplitude 1, amplitude 2, amplitude 3, amplitude 4, position       
            }
        }
        
        continuous_transition_types = [
            ['empirical_movement', 'uniform',  'uniform',   'uniform', 'uniform', 'uniform'],
            ['uniform',            'identity', 'uniform',   'uniform', 'uniform', 'uniform'],
            ['uniform',            'uniform',  'uniform',   'uniform', 'uniform', 'uniform'],
            ['uniform',            'uniform',  'uniform',   'empirical_movement', 'uniform', 'uniform'],
            ['uniform',            'uniform',  'uniform',   'uniform',            'identity', 'uniform'],
            ['uniform',            'uniform',  'uniform',   'uniform',            'uniform', 'uniform'],
            ]
        encoding_group_to_state = ['Inbound', 'Inbound', 'Inbound', 'Outbound', 'Outbound', 'Outbound']
        
        dt = 1/data.attrs['fs']
        state_diag_val = 1-1/(state_decay_timescale/dt)
        classifier = ClusterlessClassifier(replay_speed=replay_speed,
                                           place_bin_size=pos_bin_size,
                                           clusterless_algorithm=clusterless_algorithm,
                                           clusterless_algorithm_params=clusterless_algorithm_params,
                                           continuous_transition_types=continuous_transition_types,
                                           discrete_transition_diag=state_diag_val)
        classifier.fit(data.position,
                       data.multiunits,
                       encoding_group_labels=data.encoding_labels,
                       encoding_group_to_state=encoding_group_to_state,
                       is_training=is_training)
        
        #  ------------ predict (decode) ------------
        results = classifier.predict(data.multiunits.isel({"time":is_testing}),
                                     data.time.isel({"time":is_testing}))
        results.attrs.update({'exp_ID':self.exp_ID,
                              'flight_idx':flight_idx,
                              'pos_bin_size':pos_bin_size,
                              'pos_std':pos_std,
                              'mark_std':mark_std,
                              'state_decay_timescale':state_decay_timescale,
                              'replay_speed':replay_speed})
        self.behave_results = results
        self.behave_classifier = classifier
        
#%%
def main():
    global exp
    exp = exp_decode('b9861_d180527')
    exp.load_data()
    # exp.decode_sleep(win_s=0.05)
if __name__ == "__main__":
    main()






#%%









#%%