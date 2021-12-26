function ripples_detect(exp_ID)

%% get exp info
exp = exp_load_data(exp_ID,'details','path','rest');
prm = PARAMS_GetAll();
active_channels = exp.details.activeChannels;
TT_to_use = find(contains(exp.details.TT_loc,{'CA1','CA3'}));
nTT = size(active_channels,1);
nCh = size(active_channels,2);

%% load LFP data
% TODO: load data from all TTs and run detection by TT for all TTs. For 
% detection by pooling tetrodes, use only TT_to_use.
[ripple, ts, fs, params, ch_valid] = LFP_load(exp_ID,'band','ripple');
[gamma , ts, fs, params, ch_valid] = LFP_load(exp_ID,'band','high_gamma');

%% get relevant session times
sleep_ti = exp_get_sessions_ti(exp_ID, 'Sleep1','Sleep2');
sleep_ti(any(isnan(sleep_ti),2),:) = []; % remove nan in case of missing sessions
immobility_ti = [sleep_ti; exp.rest.ti];
is_sleep = any(ts>sleep_ti(:,1)&ts<sleep_ti(:,2),1);
is_immobility = any(ts>immobility_ti(:,1)&ts<immobility_ti(:,2),1);
is_rest = any(ts>exp.rest.ti(:,1)&ts<exp.rest.ti(:,2),1);

%% bands power by hilbert envelope (per channel)
pripple = abs(hilbert(ripple)).^2;
pgamma = abs(hilbert(gamma)).^2;

%% average over ch/TT + smoothing
pripple_TT = squeeze(nanmean(pripple,3)); % sum over channels per TT
pripple_all = squeeze(nanmean(pripple(:,TT_to_use,:),[2 3])); % sum over channels from RELEVANT TTs
pripple_TT = smoothdata(pripple_TT,1, 'gaussian',5*round(fs*prm.ripples.smooth_ker*1e-3), 'includenan');
pripple_all = smoothdata(pripple_all, 'gaussian',5*round(fs*prm.ripples.smooth_ker*1e-3), 'includenan');
pripple_TT  = sqrt(pripple_TT);
pripple_all = sqrt(pripple_all);

pgamma_TT = squeeze(nanmean(pgamma,3)); % sum over channels per TT
pgamma_all = squeeze(nanmean(pgamma(:,TT_to_use,:),[2 3])); % sum over channels from RELEVANT TTs
pgamma_TT = smoothdata(pgamma_TT,1, 'gaussian',5*round(fs*prm.ripples.smooth_ker*1e-3), 'includenan');
pgamma_all = smoothdata(pgamma_all, 'gaussian',5*round(fs*prm.ripples.smooth_ker*1e-3), 'includenan');
pgamma_TT  = sqrt(pgamma_TT);
pgamma_all = sqrt(pgamma_all);

%% zscore (using mean and std from sleep only)
pripple_mean_TT = mean(pripple_TT(is_sleep,:));
pripple_std_TT = std(pripple_TT(is_sleep,:));
pripple_mean_all = mean(pripple_all(is_sleep,:));
pripple_std_all = std(pripple_all(is_sleep,:));
zpripple_TT = (pripple_TT - pripple_mean_TT) ./ pripple_std_TT;
zpripple_all = (pripple_all - pripple_mean_all) ./ pripple_std_all;

%% ripple/gamma power ratio
ripple_gamma_ratio_TT = pripple_TT ./ pgamma_TT;
ripple_gamma_ratio_all = pripple_all ./ pgamma_all;

%% detect!
% TODO: speed up!
ripples_TT = {};
for TT = 1:nTT
    zpripple = zpripple_TT(:,TT);
    zpripple(~is_immobility) = nan;
    ripples_TT{TT} = detect_ripples(zpripple, ts, ripple_gamma_ratio_TT(:,TT), ...
        prm.ripples.high_thr_std, prm.ripples.low_thr_std, prm.ripples.min_width_msec, prm.ripples.merge_thr_msec, prm.ripples.ripple_gamma_power_ratio_thr);
end
zpripple = zpripple_all;
zpripple(~is_immobility) = nan;
ripples_all = detect_ripples(zpripple, ts, ripple_gamma_ratio_all, ...
        prm.ripples.high_thr_std, prm.ripples.low_thr_std, prm.ripples.min_width_msec, prm.ripples.merge_thr_msec, prm.ripples.ripple_gamma_power_ratio_thr);

%% select representative TT (for display etc...)
contrib_events_per_TT = zeros(length(ripples_all),nTT);
t1 = [ripples_all.peak_ts];
for TT=TT_to_use
    events_TT = ripples_TT{TT};
    t2 = [events_TT.peak_ts];
    tdiff = abs(t1-t2');
    thr = 50e3; % 50ms
    contrib_events_per_TT(:,TT) = any(tdiff < thr);
end
num_contrib_events_per_TT = sum(contrib_events_per_TT);
prc_contrib_events_per_TT = 100 .* num_contrib_events_per_TT ./ length(ripples_all);

IX = [ripples_all.peak_IX];
mean_pripple_per_TT = mean(pripple_TT(IX,:));
[~,best_TT] = max(num_contrib_events_per_TT);

%% get LFP traces using selected TT
[LFP, ~, ~, ~, ~] = LFP_load(exp_ID,best_TT);
LFP_raw = nanmean(LFP(:,best_TT,:),3);
LFP_filt = nanmean(ripple(:,best_TT,:),3);

%%
% figure
% hold on
% y = LFP_raw;
% % y = LFP_filt;
% plot(ts, y)
% plot([events_TT.peak_ts],y([events_TT.peak_IX]), 'r*')
% rescale_plot_data('x',[1e-6 ts(1)]);

%% save ripples detection results to mat file
ripples = struct();
ripples.params = prm.ripples;
ripples.by_TT = ripples_TT;
ripples.all = ripples_all;
ripples.zpripple_all = zpripple_all;
ripples.t = ts;
ripples.fs = params.SamplingFrequency;
ripples.LFP_raw = LFP_raw;
ripples.LFP_filt = LFP_filt;
ripples.stats.mean_pripple_per_TT = mean_pripple_per_TT;
ripples.stats.prc_contrib_events_per_TT = prc_contrib_events_per_TT;
ripples.stats.best_TT = best_TT;
ripples.stats.pripple_mean_TT = pripple_mean_TT;
ripples.stats.pripple_std_TT = pripple_std_TT;
ripples.stats.pripple_mean_all = pripple_mean_all;
ripples.stats.pripple_std_all = pripple_std_all;
file_name = fullfile('L:\Analysis\Results\exp\ripples',[exp_ID '_exp_ripples']);
save(file_name,'ripples');

end


%%
function ripples = detect_ripples(zpripple, ts, ripple_gamma_ratio, high_thr, low_thr, min_duration, merge_thr, ratio_thr)
    fs = 1e6 / median(diff(ts));
    zpripple = zpripple(:)';
    ripple_gamma_ratio = ripple_gamma_ratio(:)';
    events  = struct();
    xthr_IX = find(zpripple > high_thr);
    xthr_ts = ts(xthr_IX);
    start_IX = [xthr_IX(1) xthr_IX(find(diff(xthr_ts) > merge_thr*1e3)+1 )              ];
    end_IX   = [           xthr_IX(find(diff(xthr_ts) > merge_thr*1e3)   )  xthr_IX(end)];
    xthr1 = false(1,length(ts));
    for ii = 1:length(start_IX)
        xthr1(start_IX(ii):end_IX(ii)) = true;
    end
    xthr2 = zpripple > low_thr;
    xthr2(xthr1) = true;
    xthr12 = extend_intervals(xthr1,xthr2);
    cc = bwconncomp(xthr12);
    % create ripples events struct
    events=struct();
    events.duration = 1e3/fs .* cellfun(@length,cc.PixelIdxList); % in ms
    events.start_IX = cellfun(@min, cc.PixelIdxList);
    events.end_IX = cellfun(@max, cc.PixelIdxList);
    events.start_ts = ts(events.start_IX);
    events.end_ts = ts(events.end_IX);
    g = bwlabel(xthr12);
    g(g==0)=nan;
    [~,max_IX] = splitapply(@max,zpripple,g); % index relative to event
    events.peak_IX = cellfun(@(IX,max_IX)(IX(max_IX)), cc.PixelIdxList, num2cell(max_IX));
    events.peak_ts = ts(events.peak_IX);
    events.peak_zpripple = zpripple(events.peak_IX);
    events.ripple_gamma_power_ratio_mean = splitapply(@mean,ripple_gamma_ratio,g);
    events.ripple_gamma_power_ratio_at_peak = splitapply(@max,ripple_gamma_ratio,g);
    
    events.valid = events.duration > min_duration &...
                    events.ripple_gamma_power_ratio_at_peak > ratio_thr;
    
    ripples = soa2aos(events);
%     invalid_ripples = ripples(~[ripples.valid]);
    ripples(~[ripples.valid])=[];
end





%%


