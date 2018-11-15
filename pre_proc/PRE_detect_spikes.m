function PRE_detect_spikes(exp_ID)

%% get exp info
[exp_path,exp_info] = DS_get_path(exp_ID);
% active_TT_channels = reshape(eval(exp_info.activeChannels)',[],1)';
active_TT_channels = eval(exp_info.activeChannels);

%% set params
dir_IN = exp_path.spikes_raw;
dir_OUT = exp_path.spikes_detection;
clear detect_params
params.thr_uV = 50;
params.ref_ch = eval(exp_info.refCh);
params.t_start_end = [];
% params.t_start_end = [74190421000	74737334000]; % b0148_d170720 pre-sleep
% params.t_start_end = [33197649000	33847970000]; % b0079_d160913 pre-sleep
params.use_neg_thr = 0;
params.TT_to_use = eval(exp_info.TT_to_use);
% params.merge_thr_crs_width = 4;
% params.lib_spike_shapes = 'library_of_acceptable_spike_shapes.mat';
% params.lib_corr_thr = 0.8;
params.active_TT_channels = active_TT_channels;

%% run!
Nlx_detect_spikes_CSC(dir_IN,dir_OUT,params)

end