function PRE_filter_LFP_bands(exp_ID, forcecalc)

%% defaults
if nargin==1; forcecalc = 0; end

%% params
bands_names = ["delta","theta","low_gamma","high_gamma","ripple"];
bands_freqs = [ 
    0.5 4;
    4 10;
    30 60;
    60 100;
    100 200];

%% get exp info
exp = exp_load_data(exp_ID,'details','path');
prm = PARAMS_GetAll();
active_channels = exp.details.activeChannels;
active_channels = active_channels';
active_channels = active_channels(:);

%% check if input files exist
if ~exist(exp.path.LFP,'dir')
    warning('LFP files do not exist!');
    return;
end

%% extract LFPs and save them
run_bands_filtering = 1;
if exist(exp.path.LFP_bands,'dir')
    if forcecalc
        % delete existing output dir
        warning('LFP bands output dir already existing and you chose to override it, deleting old LFP bands dir!')
        rmdir(exp.path.LFP_bands,'s');
    else
        warning('LFP bands output folder already exist, use forcecalc to override it!');
        run_bands_filtering = 0;
    end
end
if run_bands_filtering
    % at this point we should not have the LFP output dir, so let's create it!
    mkdir(exp.path.LFP_bands)
    fprintf('Filtering LFPs for %s', exp_ID)
    
    for ii_band = 1:length(bands_freqs)
        t_start_end = [];
        clear filter_params
        filter_params.passband  = bands_freqs(ii_band,:);

        parfor ii_ch = 1:length(active_channels)
            TT = ceil(ii_ch/4);
            ch_num = mod(ii_ch-1,4)+1;
            file_IN = fullfile(exp.path.LFP,['LFP_' exp_ID '_TT' num2str(TT) '_ch' num2str(ch_num) '.ncs'])
            file_OUT = fullfile(exp.path.LFP_bands,bands_names(ii_band), ['LFP_' exp_ID '_TT' num2str(TT) '_ch' num2str(ch_num) '.ncs'])
            file_OUT  = str2mat(file_OUT);

            Nlx_filter_CSC2(file_IN, file_OUT, t_start_end, filter_params)
        end
    end
    
end

end
