function exp_calc_CSC_RAW_stats(exp_ID)

%% load exp details
exp = exp_load_data(exp_ID, 'details', 'path');
dir_IN = exp.path.spikes_raw;
% limits_ts = exp.details.session_ts(1,:);
limits_ts = [];

%% load ref channel
if ~isempty(exp.details.refCh)
    TT = exp.details.refCh(1);
    ch = exp.details.refCh(2);
    csc_file = sprintf('spikes_%s_TT%d_ch%d.ncs',exp_ID,TT,ch);
    csc_file = fullfile(dir_IN,csc_file);
    [csc_ref_ch,~] = Nlx_csc_read(csc_file, limits_ts);
end

%%
total_time = tic;
csc_std = nan(size(exp.details.activeChannels));
csc_abs_median = nan(size(exp.details.activeChannels));
csc_reref_std = nan(size(exp.details.activeChannels));
csc_reref_abs_median = nan(size(exp.details.activeChannels));
parfor TT = exp.details.TT_to_use
    for ch = 1:4%size(exp.details.activeChannels,2)
%         fprintf('TT%d_ch%d\n',TT,ch);
        % load ch data
        csc_file = sprintf('spikes_%s_TT%d_ch%d.ncs',exp_ID,TT,ch);
        csc_file = fullfile(dir_IN,csc_file);
        if ~exist(csc_file,'file')
            continue;
        end
        [csc, timestamps, fs] = Nlx_csc_read(csc_file,limits_ts);
        % calc ch data stats
        csc_std(TT,ch) = std(csc);
        csc_abs_median(TT,ch) = median(abs(csc));
        if ~isempty(exp.details.refCh)
            csc_reref_std(TT,ch) = std(csc-csc_ref_ch);
            csc_reref_abs_median(TT,ch) = median(abs(csc-csc_ref_ch));
        end
    end
end
toc(total_time)

%% arrange results
csc_raw_stats = struct();
csc_raw_stats.csc_std = csc_std;
csc_raw_stats.csc_abs_median = csc_abs_median;
csc_raw_stats.csc_reref_std = csc_reref_std;
csc_raw_stats.csc_reref_abs_median = csc_reref_abs_median;
csc_raw_stats.refCh = exp.details.refCh;

%% save strcut
dir_OUT = 'L:\Analysis\Results\exp\csc_raw_stats';
if ~exist(dir_OUT,'dir')
    mkdir(dir_OUT);
end
file_out = fullfile(dir_OUT, [exp_ID '_exp_csc_raw_stats']);
save(file_out, 'csc_raw_stats');

end





%%
