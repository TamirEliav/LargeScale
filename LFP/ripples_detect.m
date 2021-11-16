function ripples_detect(exp_ID)

%% TODO: add somewhere here a section that selects the best tetrode i.e. with the strongest 
% ripples and saves the average LFP across valid channel from this tetrode.
% alternatievly, that the best SINGLE CHANNEL instead of TETRODE.

%% get exp info
exp = exp_load_data(exp_ID,'details','path','rest');
prm = PARAMS_GetAll();
active_channels = exp.details.activeChannels;
sleep_ti = exp_get_sessions_ti(exp_ID, 'Sleep1','Sleep2');
sleep_ti(any(isnan(sleep_ti),2),:) = []; % remove nan in case of missing sessions
TT_to_use = find(contains(exp.details.TT_loc,{'CA1','CA3'}));
nTT = size(active_channels,1);
nCh = size(active_channels,2);

%% load data
% TODO: load data from all TTs and run detection by TT for all TTs. For 
% detection by pooling tetrodes, use only TT_to_use
[ripple, ts, fs, params, ch_valid] = LFP_load(exp_ID,'ripple',TT_to_use);
[gamma , ts, fs, params, ch_valid] = LFP_load(exp_ID,'high_gamma',TT_to_use);

%% bands power by hilbert envelope (per channel)
pripple = abs(hilbert(ripple)).^2;
pgamma = abs(hilbert(gamma)).^2;

%%
pripple_TT = squeeze(nansum(pripple,3)); % sum over channels per TT
pripple_all = squeeze(nansum(pripple,[2 3])); % sum over channels and TTs
pripple_TT = smoothdata(pripple_TT,1, 'gaussian',5*round(fs*prm.ripples.smooth_ker*1e-3), 'includenan');
pripple_all = smoothdata(pripple_all, 'gaussian',5*round(fs*prm.ripples.smooth_ker*1e-3), 'includenan');
pripple_TT  = sqrt(pripple_TT);
pripple_all = sqrt(pripple_all);

% TODO: calc STD for each TT and save it to the final strcut

pgamma_TT = squeeze(nansum(pgamma,3)); % sum over channels per TT
pgamma_all = squeeze(nansum(pgamma,[2 3])); % sum over channels and TTs
pgamma_TT = smoothdata(pgamma_TT,1, 'gaussian',5*round(fs*prm.ripples.smooth_ker*1e-3), 'includenan');
pgamma_all = smoothdata(pgamma_all, 'gaussian',5*round(fs*prm.ripples.smooth_ker*1e-3), 'includenan');
pgamma_TT  = sqrt(pgamma_TT);
pgamma_all = sqrt(pgamma_all);

%% zscore
switch 2
    case 1
        is_sleep = any(ts>sleep_ti(:,1)&ts<sleep_ti(:,2),1);
        pripple_TT(~is_sleep,:) = nan;
        pripple_all(~is_sleep,:) = nan;
        pgamma_TT(~is_sleep,:) = nan;
        pgamma_all(~is_sleep,:) = nan;

        zpripple_TT = nanzscore(pripple_TT,0,1);
        zpripple_all = nanzscore(pripple_all,0,1);
        zpgamma_TT = nanzscore(pgamma_TT,0,1);
        zpgamma_all = nanzscore(pgamma_all,0,1);

    case 2
        %% zscore (using mean and std from sleep only, but applied also to rest)
        is_sleep = any(ts>sleep_ti(:,1)&ts<sleep_ti(:,2),1);
        immobility_ti = [sleep_ti; exp.rest.ti];
        is_immobility = any(ts>immobility_ti(:,1)&ts<immobility_ti(:,2),1);
        is_rest = any(ts>exp.rest.ti(:,1)&ts<exp.rest.ti(:,2),1);
%         my_zscore = @(x,IX)( (x-mean(x(IX)))./std(x(IX)) );
        zpripple_TT = my_zscore(pripple_TT, is_sleep, is_immobility);
        zpripple_all = my_zscore(pripple_all, is_sleep, is_immobility);
        zpgamma_TT = my_zscore(pgamma_TT, is_sleep, is_immobility);
        zpgamma_all = my_zscore(pgamma_all, is_sleep, is_immobility);
end

%% ripple/gamma power ratio
ripple_gamma_ratio_TT = pripple_TT ./ pgamma_TT;
ripple_gamma_ratio_all = pripple_all ./ pgamma_all;

%% detect!
% TODO: speed up!
clear ripples_TT ripples_TT_invalid
for TT = 1:nTT
    if ismember(TT,TT_to_use)
        [ripples_TT{TT}, ripples_TT_invalid{TT}] = detect_ripples(zpripple_TT(:,TT), ts, ripple_gamma_ratio_TT(:,TT), ...
            prm.ripples.high_thr_std, prm.ripples.low_thr_std, prm.ripples.min_width_msec, prm.ripples.merge_thr_msec, prm.ripples.ripple_gamma_power_ratio_thr);
    else
        ripples_TT{TT} = [];
        ripples_TT_invalid{TT} = [];
    end
end
[ripples_all, ripples_all_invalid] = detect_ripples(zpripple_all, ts, ripple_gamma_ratio_all, ...
        prm.ripples.high_thr_std, prm.ripples.low_thr_std, prm.ripples.min_width_msec, prm.ripples.merge_thr_msec, prm.ripples.ripple_gamma_power_ratio_thr);

%% save ripples detection results to mat file
ripples = struct();
ripples.params = prm.ripples;
ripples.by_TT = ripples_TT;
ripples.by_TT_invalid = ripples_TT_invalid;
ripples.all = ripples_all;
ripples.all_invalid = ripples_all_invalid;
ripples.zpripple_all = zpripple_all;
ripples.t = ts;
ripples.fs = params.SamplingFrequency;
file_name = fullfile('L:\Analysis\Results\exp\ripples',[exp_ID '_exp_ripples']);
save(file_name,'ripples');

end


%%
function [ripples, invalid_ripples] = detect_ripples(zpripple, ts, ripple_gamma_ratio, high_thr, low_thr, min_duration, merge_thr, ratio_thr)
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
    invalid_ripples = ripples(~[ripples.valid]);
    ripples(~[ripples.valid])=[];
end


%%
function z = my_zscore(x,IX1,IX2)
    z = (x-mean(x(IX1))) ./ std(x(IX1));
    z(~IX2) = nan;
end




%%


