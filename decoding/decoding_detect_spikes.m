function decoding_detect_spikes(exp_ID,forcecalc)

%% defaults
if nargin==1; forcecalc = 0; end

%% get exp info
exp = exp_load_data(exp_ID,'details','path');

%% set params
dir_IN = exp.path.spikes_raw;
dir_OUT = exp.path.decoding_spikes_detection;
clear params
params.ref_ch = exp.details.refCh;
params.active_TT_channels = exp.details.activeChannels;
params.TT_to_use = exp.details.TT_to_use;
params.t_start_end = [];
% params.t_start_end = exp_get_sessions_ti(exp_ID,'Sleep1');
% params.thr_type = 'uVolt';
params.thr_type = 'median'; % using median(abs(signal)) * factor
% params.thr = 50;
params.thr = 6;
params.use_neg_thr = 0;
% params.merge_thr_crs_width = 4;
params.lib_spike_shapes = 'library_of_acceptable_spike_shapes_new.mat';
params.lib_corr_thr = 0.8;
params.min_sep_events = 24;
params.CD_thr = 6;
% params.CD_detect_win_len = 32;
params.CD_detect_win_len = 4;
params.CD_invalid_win_len = 32*2;
if length(exp.details.TT_to_use) >= 8
    params.CD_n_TT_thr = round(0.5*length(exp.details.TT_to_use));
else
    params.CD_n_TT_thr = min(4, length(exp.details.TT_to_use));
end
params.is_save_artifacts = 0;

%% detect!
Nlx_detect_spikes_CSC3(dir_IN,dir_OUT,params,forcecalc);

end