%%
close all
clear
clc

%% load data
behave_dir_IN = 'L:\2bats\bat2299_behavioral_data_from_shaked\';
out_dir = 'L:\2bats\processed_data_adjusted\';
files = dir(fullfile(behave_dir_IN,'*2299*.mat'));
varNames = {'exp_ID','session_names','session_ts'};
varTypes = {'string','string','string'};
T = table('Size',[length(files) length(varNames)], 'VariableNames',varNames,'VariableTypes',varTypes);
for ii_file = 1:length(files)
    %% load data
    file = files(ii_file);
    str = regexp(file.name,'\d+','match');
    bat = str{1};
    day = str{2};
    exp_ID = sprintf('b%s_d%s',bat,day(3:end));
    load(fullfile(behave_dir_IN,sprintf('behav_bat_%s_day_%s.mat',bat,day)));

    %% arrange data
    pos = struct();
    pos.raw=[];
    pos.calib_tunnel=[];
    pos.proj=[];
    pos.proc_1D=struct();
    pos.balls=[];

    bsp = b.bsp(find([b.bsp.self]));
    pos.proc_1D.pos = bsp.pos_upsampled(1,:);
    pos.proc_1D.ts = bsp.ts_synced;
    pos.proc_1D.fs = bsp.fs_upsampled;

    % resample pos with continuous ts of fixed intervals (and place nans in
    % missing ts gaps)
    pos.proc_1D = pos_resample(pos.proc_1D);

    % Calc velocity (using csaps smoothing)
    prm = PARAMS_GetAll();
    pos.proc_1D.vel = [0 diff(pos.proc_1D.pos).*prm.pos.resample_fs];
    pos.proc_1D.pos_csaps_pp = csaps(1:length(pos.proc_1D.ts), pos.proc_1D.pos', prm.pos.csaps_p);
    pos.proc_1D.pos_csaps_p = prm.pos.csaps_p;
    pos.proc_1D.pos_csaps = fnval(pos.proc_1D.pos_csaps_pp, 1:length(pos.proc_1D.ts) );
    pos.proc_1D.vel_csaps = fnval( fnder(pos.proc_1D.pos_csaps_pp), 1:length(pos.proc_1D.ts)) .* prm.pos.resample_fs;
    pos.proc_1D.pos_csaps(isnan(pos.proc_1D.pos)) = nan;
    pos.proc_1D.vel_csaps(isnan(pos.proc_1D.pos)) = nan;

    % cross-overs
    co_ts = cat(1,b.co.ts);
    [co_ts, co_sort_IX] = sort(co_ts,'ascend');
    co_pos = interp1(b.bsp(1).ts_synced, b.bsp(1).pos_upsampled(1,:), co_ts, 'linear','extrap');
    gDir = arrayfun(@(co,d)(d.*ones(1,length(co.ts))),b.co,1:length(b.co),'UniformOutput',false);
    gDir = [gDir{:}];
    gDir = gDir(co_sort_IX);
    if ~isfield(b.co,'click_rate_attention_range')
        for ii_dir = 1:length(b.co)
            b.co(ii_dir).click_rate_attention_range = nan(size(b.co(ii_dir).ts));
            b.co(ii_dir).excluded_from_audio_analysis = ones(size(b.co(ii_dir).ts));
        end
    end
    co_click_rate = cat(1,b.co.click_rate_attention_range);
    co_click_rate = co_click_rate(co_sort_IX);
    co_excluded_from_audio_analysis = cat(1,b.co.excluded_from_audio_analysis);
    co_excluded_from_audio_analysis = co_excluded_from_audio_analysis(co_sort_IX);
    pos.proc_1D.co.ts = co_ts;
    pos.proc_1D.co.pos = co_pos;
    pos.proc_1D.co.direction = gDir';
    pos.proc_1D.co.click_rate = co_click_rate;
    pos.proc_1D.co.excluded_from_audio_analysis = co_excluded_from_audio_analysis;
    
    bsp = b.bsp(find(~[b.bsp.self]));
    pos.proc_1D.other.pos = bsp.pos_upsampled;
    pos.proc_1D.other.ts = bsp.ts_synced;

    %% save position data into mat file
    pos_file_out = fullfile(out_dir,sprintf('%s_exp_pos.mat',exp_ID));
    save(pos_file_out,'pos');

    %% add details to the table
    sessions_name = {b.events.session};
    sessions_ti = [b.events.start_time;b.events.end_time];
    T(ii_file,"exp_ID") = {exp_ID};
    T(ii_file,"session_names") = {sprintf('%s; ',sessions_name{:})};
    T(ii_file,"session_ts") = {['[' sprintf('%d %d; ',sessions_ti) ']']};
    %% save data
end
% write exp details table into an excel file
writetable(T,fullfile(out_dir,'2bats_exp_details.xlsx'));

return

%% check data structure
clear 
clc
load('G:\Tamir\work\PROJECTS\LargeScale\DATA\2299_2299\bat2299_behavioral_data_from_shaked\behav_bat_2299_day_20191201.mat')
co_ts = cat(1,b.co{:});
[co_ts co_sort_IX] = sort(co_ts,'ascend');
co_pos = interp1(b.bsp(1).ts_synced, b.bsp(1).pos_upsampled(1,:), co_ts, 'linear','extrap');
gDir = arrayfun(@(x,d)(d.*ones(1,length(x{:}))),b.co,1:length(b.co),'UniformOutput',false);
gDir = [gDir{:}];
gDir = gDir(co_sort_IX);
co_clrs = interp1([1 2],[1 0 0; 0 0 1],gDir,'linear');
sessions_name = {b.events.session};
sessions_ti = [b.events.start_time;b.events.end_time]';

close all
figure
hold on
for ii = 1:length(b.bsp)
    sdf = b.bsp(ii);
    plot(sdf.ts_synced, sdf.pos_upsampled(1,:),'.');
end
% plot(co_ts,co_pos,'*r')
scatter(co_ts,co_pos,100,co_clrs,'*')

%%
for ii_session = 1:length(sessions_name)
    area(sessions_ti(ii_session,:),[140 140],'FaceAlpha',0.2);
end
legend("self="+[b.bsp.self])
xline(56952012172)
rescale_plot_data('x',[1e-6/60 min([b.bsp.ts_synced])])











%% LFP stuff...
exp = exp_load_data(exp_ID,'details','path');
TT=1;
LFP = struct();
[LFP.signal, LFP.ts, LFP.fs, LFP.params, LFP.ch_valid] = LFP_load(exp_ID,[TT],'band','delta');
LFP.phase = angle(hilbert(LFP.signal(:,TT,1)));
figure
tiledlayout(4,4);
m=[];
for TT=exp.details.TT_to_use
    filename = fullfile(exp.path.decoding_spikes_detection,sprintf('spikes_%s_TT%d.NTT',exp_ID,TT));
    ti = sessions_ti(1,:);
    NTT = nlx_ntt_read(filename,ti);
    
    spikes_phases = interp1(LFP.ts,LFP.phase,NTT.ts,'linear','extrap');
    
    nexttile;
    polarhistogram(spikes_phases);
    hold on
    m(TT) = circ_mean(spikes_phases');
    polarplot([m(TT) m(TT)],rlim,'r');
    title("TT"+TT)
end
figure
polarhistogram(m)













%% 
function pos = pos_resample(pos)

dt = 1e6/pos.fs;
pos2 = struct();
pos2.ts = pos.ts(1):dt:pos.ts(end);
pos2.pos = interp1(pos.ts,pos.pos,pos2.ts,'linear','extrap');
nearest_ts_IX = interp1(pos.ts, 1:length(pos.ts), pos2.ts, 'nearest');
tdiff = pos2.ts - pos.ts(nearest_ts_IX);
valid_IX = abs(tdiff)<2*dt; % max two samples away sounds good enough
pos2.pos(~valid_IX)=nan;
pos2.fs = pos.fs;
pos = pos2;

end

