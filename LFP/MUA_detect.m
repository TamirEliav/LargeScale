function MUA_detect(exp_ID)

%% get exp info
exp = exp_load_data(exp_ID,'details','path');
prm = PARAMS_GetAll();
active_channels = exp.details.activeChannels;
nTT = exp.details.numTT;
sleep_ti = exp_get_sessions_ti(exp_ID, 'Sleep1','Sleep2');
sleep_ti(any(isnan(sleep_ti),2),:) = []; % remove nan in case of missing sessions
TT_to_use = find(contains(exp.details.TT_loc,{'CA1','CA3'}));

%% binning params
tstart = min(exp.details.session_ts(:));
tend = max(exp.details.session_ts(:));
bin_size = prm.MUA.bin_size * 1e3; % usec
fs=1e6/bin_size;
edges = tstart:bin_size:tend;
t = edges2centers(edges);

%% load spike detection + bin
N = nan(nTT,length(t));
for TT=1:nTT
    if ~ismember(TT,TT_to_use)
        continue;
    end
    TT_ID = str2mat(exp.details.exp_ID+"_TT"+TT);
    NTT_file = fullfile(exp.path.decoding_spikes_detection, ['spikes_' TT_ID '.NTT']);
	[Timestamps, CellNumbers, Samples, Header] = Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 1, [] );
    
    %% parse header (convert bits to uVolts)
    ADBitVolts = sscanf(Header{contains(Header,'ADBitVolts')},'-ADBitVolts %f %f %f %f');
    Samples = Samples .* ADBitVolts' .* 1e6;

    %% binning
    N(TT,:) = histcounts(Timestamps,edges);
end

%% calc smoothed total FR
ker_SD_n_samples = prm.MUA.smooth_ker*1e3/bin_size;
ker_win_n_samples = ker_SD_n_samples * 5; % smoothdata uses 1/5 of the win size to be the SD
FR = fs.*smoothdata(nanmean(N),'gaussian',ker_win_n_samples);

%% detect! (only during immobility)
is_sleep = any(t>sleep_ti(:,1)&t<sleep_ti(:,2),1);
zFR = FR;
zFR(~is_sleep) = nan;
zFR = nanzscore(zFR);
xthr1 = zFR > prm.MUA.high_thr_std;
xthr2 = zFR > prm.MUA.low_thr_std;
xthr12 = extend_intervals(xthr1,xthr2);
cc = bwconncomp(xthr12);
% create MUA events struct
events = struct();
events.duration = 1e3/fs .* cellfun(@length,cc.PixelIdxList); % in ms
events.start_IX = cellfun(@min, cc.PixelIdxList);
events.end_IX = cellfun(@max, cc.PixelIdxList);
events.start_ts = t(events.start_IX);
events.end_ts = t(events.end_IX);
g = bwlabel(xthr12);
g(g==0)=nan;
[~,max_IX] = splitapply(@max,zFR,g); % index relative to event
events.peak_IX = cellfun(@(IX,max_IX)(IX(max_IX)), cc.PixelIdxList, num2cell(max_IX));
events.peak_ts = t(events.peak_IX);
events.peak_FR = FR(events.peak_IX);
events.peak_zFR = zFR(events.peak_IX);

%% old code for extending the high to low thr
% % % % % extend the high thr crossing to the low thr crossing
% % % % cc1 = bwconncomp(xthr1);
% % % % cc2 = bwconncomp(xthr2);
% % % % thr2_ti =   [cellfun(@min, cc2.PixelIdxList)
% % % %              cellfun(@max, cc2.PixelIdxList)]';
% % % % [~, IX_per_ti] = get_data_in_ti(find(xthr1),thr2_ti);
% % % % cc2_1_TF = cellfun(@length,IX_per_ti)~=0;
% % % % xthr12 = zeros(size(t));
% % % % xthr12_IX = cat(1,cc2.PixelIdxList{cc2_1_TF});
% % % % xthr12(xthr12_IX) = 1;
% % % % cc12 = bwconncomp(xthr12);
% % % % % create MUA events struct
% % % % events=struct();
% % % % events.duration = 1e3/fs .* cellfun(@length,cc12.PixelIdxList); % in ms
% % % % events.start_IX = cellfun(@min, cc12.PixelIdxList);
% % % % events.end_IX = cellfun(@max, cc12.PixelIdxList);
% % % % events.start_ts = t(events.start_IX);
% % % % events.end_ts = t(events.end_IX);
% % % % g = bwlabel(xthr12);
% % % % g(g==0)=nan;
% % % % [~,max_IX] = splitapply(@max,zFR,g); % index relative to event
% % % % events.peak_IX = cellfun(@(IX,max_IX)(IX(max_IX)), cc12.PixelIdxList, num2cell(max_IX));
% % % % events.peak_ts = t(events.peak_IX);
% % % % events.peak_FR = FR(events.peak_IX);
% % % % events.peak_zFR = zFR(events.peak_IX);

%% MUA event triggered FR
trig_win = 1e3;
trig_FR = trigger_signal_by_IX(FR, events.peak_IX, trig_win);
trig_zFR = trigger_signal_by_IX(zFR, events.peak_IX, trig_win);
trig_FR_mean = nanmean(trig_FR,1);
trig_zFR_mean = nanmean(trig_zFR,1);
trig_FR_sem = nanstd(trig_FR,0,1) ./ sqrt(size(trig_FR,1));
trig_zFR_sem = nanstd(trig_zFR,0,1) ./ sqrt(size(trig_zFR,1));
trig_t = linspace(-trig_win,trig_win,size(trig_FR,2));

%% save MUA event detection results
MUA = struct();
MUA.t  = t;
MUA.fs = fs;
MUA.FR = FR;
MUA.zFR = zFR;
MUA.events = soa2aos(events);
MUA.trig.t = trig_t;
MUA.trig.FR_mean = trig_FR_mean;
MUA.trig.FR_sem  = trig_FR_sem;
MUA.trig.zFR_mean = trig_zFR_mean;
MUA.trig.zFR_sem  = trig_zFR_sem;
file_name = fullfile('L:\Analysis\Results\exp\MUA',[exp_ID '_exp_MUA']);
save(file_name,'MUA');

end

%%
function xthr_high_low = extend_intervals(xthr_high,xthr_low)

%% testing
% xthr_high = [0 0 1 0 0 0 1 0 0 0];
% xthr_low  = [0 1 1 1 0 0 1 1 0 1];

%%
xthr_low_lbl = bwlabel(xthr_low);
xthr_low_high_lbl_valid = xthr_low_lbl .* xthr_high;
valid_lbl = unique(xthr_low_high_lbl_valid(xthr_low_high_lbl_valid>0));
xthr_high_low = xthr_low;
xthr_high_low(~ismember(xthr_low_lbl,valid_lbl))=0;

end










%%
