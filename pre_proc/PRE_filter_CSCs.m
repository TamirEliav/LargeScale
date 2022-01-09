function PRE_filter_CSCs(exp_ID, forcecalc, opts)
arguments 
    exp_ID
    forcecalc = 0
    opts.passband_LFP = [0.5 400]
    opts.passband_spikes = [600 6000]
    opts.LFP_resamlpe_fs = 2000
    opts.run_LFP_filtering = 1
    opts.run_SPIKES_filtering = 1
end

%% get exp info
exp = exp_load_data(exp_ID,'details','path');
active_channels = exp.details.activeChannels;
active_channels = active_channels';
active_channels = active_channels(:);
is_CSC_name_by_TT = ~isempty(dir(fullfile(exp.path.nlx,'CSC_TT*')));

%% first, check what we need to filter
run_LFP_filtering = opts.run_LFP_filtering;
if run_LFP_filtering && exist(exp.path.LFP,'dir')
    if forcecalc
        % delete existing output dir
        warning('LFP output dir already existing and you chose to override it, deleting old LFP dir!')
        rmdir(exp.path.LFP,'s')
    else
        warning('LFP output folder already exist, use forcecalc to override it!');
        run_LFP_filtering = 0;
    end
end
run_SPIKES_filtering = opts.run_SPIKES_filtering;
if run_SPIKES_filtering && exist(exp.path.spikes_raw,'dir')
    if forcecalc
        % delete existing output dir
        warning('spikes_raw output dir already existing and you chose to override it, deleting old spikes_raw dir!')
        rmdir(exp.path.spikes_raw,'s')
    else
        warning('spikes_raw output folder already exist, use forcecalc to override it!');
        run_SPIKES_filtering = 0;
    end
end

if ~run_LFP_filtering && ~run_SPIKES_filtering 
    return
end

%% create out folders
mkdir(exp.path.LFP);
mkdir(exp.path.spikes_raw);

%% create filters
t_start_end = [];
clear LFP_filter_params spikes_filter_params
LFP_filter_params.passband  = opts.passband_LFP;
LFP_filter_params.resample_fs = opts.LFP_resamlpe_fs;
switch 3
    case 1 % original from Nachum post-doc
        spikes_filter_params.passband  = opts.passband_spikes;
    case 2 % from Aronov
        spikes_filter_params.type = 'custom';
        stopfreq = 750;
        passfreq = 1000;
        spikes_filter_params.designfilt_input = {'highpassfir',...
            'StopbandFrequency', stopfreq,...
            'PassbandFrequency', passfreq,...
            'StopbandAttenuation', 60,...
            'PassbandRipple', 1};
    case 3 % simply the original filter without the lowpass (no cutoff @ 6kHz)
        spikes_filter_params.type = 'highpassfir1';
        spikes_filter_params.passband  = opts.passband_spikes(1);
end

%% filter!
fprintf('Start filtering raw data for %s\n', exp_ID);
parfor ii_ch = 1:length(active_channels)
    TT = ceil(ii_ch/4);
    ch_num = mod(ii_ch-1,4)+1;
    if is_CSC_name_by_TT
        file_IN = fullfile(exp.path.nlx,['CSC_TT' num2str(TT) '_' num2str(ch_num) '.ncs']);
    else
        file_IN = fullfile(exp.path.nlx,['CSC' num2str(ii_ch-1) '.ncs']);
    end
    LFP_file_OUT = fullfile(exp.path.LFP,           ['LFP_'    exp_ID '_TT' num2str(TT) '_ch' num2str(ch_num) '.ncs'])
    spikes_file_OUT = fullfile(exp.path.spikes_raw, ['spikes_' exp_ID '_TT' num2str(TT) '_ch' num2str(ch_num) '.ncs'])

    if run_LFP_filtering 
        Nlx_filter_CSC2(file_IN, LFP_file_OUT, t_start_end, LFP_filter_params);
    end
    if run_SPIKES_filtering
        Nlx_filter_CSC2(file_IN, spikes_file_OUT, t_start_end, spikes_filter_params);
    end
end




end

