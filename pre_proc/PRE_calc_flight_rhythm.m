function PRE_calc_flight_rhythm(exp_ID)

%% get exp info
[exp_path exp_info] = DS_get_path(exp_ID);
p = PARAMS_GetAll();

%%
activeChannels = eval(exp_info.activeChannels)';
ch = find(~activeChannels,1);
if isempty(ch)
	ch = find(activeChannels,1);
end
% convert to TT channel mapping
TT = floor((ch-1)/4)+1;
ch = mod((ch-1),4)+1;
file_IN = fullfile(exp_path.LFP, sprintf('LFP_%s_TT%d_ch%d.ncs', exp_ID,TT,ch));
file_OUT = fullfile(exp_path.LFP, ['flight_rhythm_' exp_ID '.ncs']);
t_start_end = [];
filter_params.passband = p.LFP.flight_rhythm.band;
filter_params.resample_fs = p.LFP.flight_rhythm.resample_fs;
Nlx_filter_CSC(file_IN, file_OUT, t_start_end, filter_params);

end