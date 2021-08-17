function PRE_filter_CSCs(exp_ID, forcecalc)

%% defaults
if nargin==1; forcecalc = 0; end

%% params
passband_spikes   = [600 6000];        % Filter for spikes
passband_LFP      = [0.5 400];         % Filter for LFPs
LFP_resamlpe_fs     = 2000;

%% get exp info
exp = exp_load_data(exp_ID,'details','path');
prm = PARAMS_GetAll();
active_channels = exp.details.activeChannels;
active_channels = active_channels';
active_channels = active_channels(:);
is_CSC_name_by_TT = ~isempty(dir(fullfile(exp.path.nlx,'CSC_TT*')));

%% extract LFPs and save them
run_LFP_filtering = 1;
if exist(exp.path.LFP,'dir')
    if forcecalc
        % delete existing output dir
        warning('LFP output dir already existing and you chose to override it, deleting old LFP dir!')
        rmdir(exp.path.LFP,'s')
    else
        warning('LFP output folder already exist, use forcecalc to override it!');
        run_LFP_filtering = 0;
    end
end
if run_LFP_filtering
    % at this point we should not have the LFP output dir, so let's create it!
    mkdir(exp.path.LFP)
    fprintf('Filtering LFPs for %s', exp_ID)
    
    t_start_end = [];
    clear filter_params
    filter_params.passband  = passband_LFP;
    filter_params.resample_fs = LFP_resamlpe_fs;
    filter_params
    
    parfor ii_ch = 1:length(active_channels)
        TT = ceil(ii_ch/4);
        ch_num = mod(ii_ch-1,4)+1;
        if is_CSC_name_by_TT
            file_IN = fullfile(exp.path.nlx,['CSC_TT' num2str(TT) '_' num2str(ch_num) '.ncs']);
        else
            file_IN = fullfile(exp.path.nlx,['CSC' num2str(ii_ch-1) '.ncs']);
        end
        file_OUT = fullfile(exp.path.LFP, ['LFP_' exp_ID '_TT' num2str(TT) '_ch' num2str(ch_num) '.ncs'])
        
        Nlx_filter_CSC2(file_IN, file_OUT, t_start_end, filter_params)
    end
end

%% exctract high-pass for spike detection
run_SPIKES_filtering = 1;
if exist(exp.path.spikes_raw,'dir')
    if forcecalc
        % delete existing output dir
        warning('spikes_raw output dir already existing and you chose to override it, deleting old spikes_raw dir!')
        rmdir(exp.path.spikes_raw,'s')
    else
        warning('spikes_raw output folder already exist, use forcecalc to override it!');
        run_SPIKES_filtering = 0;
    end
end
if run_SPIKES_filtering
    % at this point we should not have the LFP output dir, so let's create it!
    mkdir(exp.path.spikes_raw)
    fprintf('Filtering spikes_raw for %s', exp_ID)
    
    t_start_end = [];
    clear filter_params
    switch 3
        case 1 % original from Nachum post-doc
            filter_params.passband  = passband_spikes;
        case 2 % from Aronov
            filter_params.type = 'custom';
            stopfreq = 750;
            passfreq = 1000;
            filter_params.designfilt_input = {'highpassfir',...
                'StopbandFrequency', stopfreq,...
                'PassbandFrequency', passfreq,...
                'StopbandAttenuation', 60,...
                'PassbandRipple', 1};
        case 3 % simply the original filter without the lowpass (no cutoff @ 6kHz)
            filter_params.type = 'highpassfir1';
            filter_params.passband  = passband_spikes(1);
    end
    
    parfor ii_ch = 1:length(active_channels)
        if ~active_channels(ii_ch)
            continue;
        end
        TT = ceil(ii_ch/4);
        ch_num = mod(ii_ch-1,4)+1;
        if is_CSC_name_by_TT
            file_IN = fullfile(exp.path.nlx,['CSC_TT' num2str(TT) '_' num2str(ch_num) '.ncs']);
        else
            file_IN = fullfile(exp.path.nlx,['CSC' num2str(ii_ch-1) '.ncs']);
        end
        file_OUT = fullfile(exp.path.spikes_raw, ['spikes_' exp_ID '_TT' num2str(TT) '_ch' num2str(ch_num) '.ncs'])
        
        Nlx_filter_CSC2(file_IN, file_OUT, t_start_end, filter_params)
    end
end


end

