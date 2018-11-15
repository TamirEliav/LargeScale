function PRE_filter_CSCs(exp_ID)

%%
% file_IN = 'L:\DATA\0148_Boson\20170608\nlx\CSC14.ncs';
% file_OUT = 'L:\Analysis\pre_proc\0148\20170608\Spikes\CSC14_filtered_600-6000Hz.ncs';
% t_start_end = [];
% filter_params.passband  = [600 6000];
% filter_params.fwin      = 2*60;
% tic
% Nlx_filter_CSC(file_IN, file_OUT, t_start_end, filter_params)
% toc

%% params
fwin = 2;         % we will run over the data in 2-min windows but save           % Same as we use for filtering ripples
passband_spikes   = [600 6000];        % Filter for spikes
passband_LFP      = [0.5 400];         % Filter for LFPs
LFP_resamlpe_fs     = 2000;

%% get exp info
[exp_path,exp_info] = DS_get_path(exp_ID);
active_TT_channels = reshape(eval(exp_info.activeChannels)',[],1)';

%% extract LFPs and save them
% skip CSC extraction and exit function if we have already extracted the data
run_LFP_filtering = 1;
% if ~exist(csc_dir_LFP,'dir')
%     mkdir(csc_dir_LFP);
%     run_LFP_filtering = 1;
% end
% if forcecalc
%     % TODO: delete/rename old files!!!
%     run_LFP_filtering = 1;
% end
% if ~param.spikes.generate_LFP_files
%     run_LFP_filtering = 0;
% end
if run_LFP_filtering
    
    t_start_end = [];
    clear filter_params
    filter_params.passband  = passband_LFP;
    filter_params.fwin      = fwin;
    filter_params.resample_fs = LFP_resamlpe_fs;
    
    parfor ii_ch = 1:length(active_TT_channels)
%         if ~active_TT_channels(ii_ch)
%             continue;
%         end
        TT = ceil(ii_ch/4);
        ch_num = mod(ii_ch-1,4)+1;
        file_IN = fullfile(exp_path.nlx,['CSC' num2str(ii_ch-1) '.ncs']);
        file_OUT = fullfile(exp_path.LFP, ['LFP_' exp_info.exp_ID '_TT' num2str(TT) '_ch' num2str(ch_num) '.ncs'])
        
        Nlx_filter_CSC2(file_IN, file_OUT, t_start_end, filter_params)
    end
    
else
    disp('======================');
    disp (['Skipping LFP extraction. Already extracted in ' exp_path.LFP]);
    disp('======================');
end

%% exctract high-pass for spike detection

run_SPIKES_filtering = 1;
% % skip CSC extraction and exit function if we have already extracted the data
% if isempty (dir(csc_dir_spikes))
%     mkdir(csc_dir_spikes);
%     run_SPIKES_filtering = 1;
% end
% if forcecalc
%     % TODO: delete/rename old files!!!
%     run_SPIKES_filtering = 1;
% end

if run_SPIKES_filtering
    
    t_start_end = [];
        clear filter_params
        filter_params.passband  = passband_spikes;
        filter_params.fwin      = fwin;
    parfor ii_ch = 1:length(active_TT_channels)
%         if ~active_TT_channels(ii_ch)
%             continue;
%         end
        TT = ceil(ii_ch/4);
        ch_num = mod(ii_ch-1,4)+1;
        file_IN = fullfile(exp_path.nlx,['CSC' num2str(ii_ch-1) '.ncs']);
        file_OUT = fullfile(exp_path.spikes_raw, ['spikes_' exp_info.exp_ID '_TT' num2str(TT) '_ch' num2str(ch_num) '.ncs'])
        
        Nlx_filter_CSC2(file_IN, file_OUT, t_start_end, filter_params)
    end
    
else
    disp('======================');
    disp (['Skipping high-pass extraction. Already extracted in ' csc_dir_spikes]);
    disp('======================');
    return;
end;






end