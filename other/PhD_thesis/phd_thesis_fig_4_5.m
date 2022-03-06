%% PhD thesis figure 4.5 - replay vs ripples/MUA/PE

%%
close all
clear 
clc

%% options
win1 = .5; % in seconds
win2 = 2; % in seconds
xcorr_hist_num_bins = 25;
% use_bat = 34;
% use_bat = 148;
% use_bat = 9861;
% use_bat = 2289;
% use_bat = 194;
% use_bat = 184;
% use_bat = 2382;
%     bat_num    GroupCount
%     _______    __________
% 
%      2289           1    
%        34           2    
%      9861           4    
%       148          15    
%       184          17    
%       194          17    
%      2382          26    

%% define output files
res_dir = 'E:\Tamir\PhD\Thesis\resources\ch_4_seq';
mkdir(res_dir)
fig_name_str = 'Fig_4_5_replay_vs_ripples';
fig_caption_str = ' ';
log_name_str = [fig_name_str '_log_file' '.txt'];
log_name_str = strrep(log_name_str , ':', '-');
log_name_str = strrep(log_name_str , ' ', '_');
log_name_out = fullfile(res_dir, log_name_str);

%% open log file
diary off
diary(log_name_out)
diary on
disp('Log file');
disp(['created: ', datestr(clock)]);
disp('======================================================');
disp([fig_name_str ':' fig_caption_str]);   
disp('======================================================');
disp('');

%% create figure
% figure_size_cm = [21.0 29.7]; % ~A4
figure_size_cm = [21.6 27.9]; % ~US letter
figure ;
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',7);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'DefaultAxesUnits','centimeters');
set(gcf,'PaperType','usletter')
% set(gcf,'PaperType','<custom>');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 figure_size_cm]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]); % position on screen...
set(gcf, 'Renderer', 'painters');
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

% create panels
panels_size = [4 4];
panels(1) = axes('position', [2 19 panels_size]);
panels(2) = axes('position', [7.5 19 panels_size]);
% panels(3) = axes('position', [13 19 panels_size]);
panels(4) = axes('position', [2 13 panels_size]);
panels(5) = axes('position', [7.5 13 panels_size]);
panels(6) = axes('position', [13 13 panels_size]);
panels(7) = axes('position', [2 7 panels_size]);
panels(8) = axes('position', [7.5 7 panels_size]);
panels(9) = axes('position', [13 7 panels_size]);
panels(10) = axes('position', [2 1 panels_size]);

%% arrange data
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
clear exp_list
groupsummary(T,'bat_num')
if exist('use_bat','var')
    T = groupfilter(T,"bat_num",@(x)x==use_bat,'bat_num');
end
bats = unique(T.bat_num)

triggered_zpripple_all = {};
triggered_MUA_zFR_all = {};
xcorr_posterior_ripples_tdiff = {};
xcorr_posterior_MUA_tdiff = {};
xcorr_posterior_PE_tdiff = {};
xcorr_posterior_ripples_fraction = {};
xcorr_posterior_MUA_fraction = {};
xcorr_posterior_PE_fraction = {};
n_ripples_per_replay = {};
seqs_all = [];
for ii_exp = 1:height(T)
    %% load exp data
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details','path','ripples','MUA','PE');
    epoch_type = 'sleep';
    params_opt = 11;
    [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    seqs = [events.seq_model];
    [seqs.ts]=disperse([events.peak_ts]);

    %% apply inclusion criteria
    seqs([seqs.score]<0.5)=[];
    seqs([seqs.distance]<3)=[];
%     [seqs.ts] = disperse(mean([seqs.start_ts; seqs.end_ts]));
    if isempty(seqs)
        continue;
    end
    seqs_all = [seqs_all seqs];

    %% trigger signals (ripples/MUA) around posterior events
    win_samples = round(win1*exp.ripples.fs);
    trigger_IX = interp1(exp.ripples.t, 1:length(exp.ripples.t), [seqs.ts], 'nearest');
    triggered_zpripple_all{ii_exp} = trigger_signal_by_IX(exp.ripples.zpripple', trigger_IX, win_samples);
    win_samples = round(win1*exp.MUA.fs);
    trigger_IX = interp1(exp.MUA.t, 1:length(exp.MUA.t), [seqs.ts], 'nearest');
    triggered_MUA_zFR_all{ii_exp} = trigger_signal_by_IX(exp.MUA.zFR, trigger_IX, win_samples);

    %% events xcorr (posterior vs ripples)
    t1 = [seqs.ts];
    t2 = [exp.ripples.events.peak_ts];
    tdiff = (t2'-t1).*1e-6;
    fraction = sum(any(abs(tdiff)<win1)) / size(tdiff,2); % fraction of the posterior events which had ripple event nearby
    tdiff = tdiff(:);
    tdiff(abs(tdiff)>win1) = [];
    xcorr_posterior_ripples_tdiff{ii_exp} = tdiff;
    xcorr_posterior_ripples_fraction{ii_exp} = fraction;

    %% events xcorr (posterior vs MUA)
    t1 = [seqs.ts];
    t2 = [exp.MUA.events.peak_ts];
    tdiff = (t2'-t1).*1e-6;
    fraction = sum(any(abs(tdiff)<win1)) / size(tdiff,2); % fraction of the posterior events which had MUA event nearby
    tdiff = tdiff(:);
    tdiff(abs(tdiff)>win1) = [];
    xcorr_posterior_MUA_tdiff{ii_exp} = tdiff;
    xcorr_posterior_MUA_fraction{ii_exp} = fraction;

    %% events xcorr (posterior vs PE)
    t1 = [seqs.ts];
    t2 = [exp.PE.thr.peak_ts];
    tdiff = (t2'-t1).*1e-6;
    fraction = sum(any(abs(tdiff)<win1)) / size(tdiff,2); % fraction of the posterior events which had population event nearby
    tdiff = tdiff(:);
    tdiff(abs(tdiff)>win1) = [];
    xcorr_posterior_PE_tdiff{ii_exp} = tdiff;
    xcorr_posterior_PE_fraction{ii_exp} = fraction;

    %% count num ripples in each replay
    ti = [seqs.start_ts; seqs.end_ts]';
    ts = [exp.ripples.events.peak_ts];
    [~,IX_per_ti] = get_data_in_ti(ts,ti);
    n_ripples_per_replay{ii_exp} = cellfun(@length,IX_per_ti);
end

%% concatanated data across sessions
triggered_zpripple_cat = cat(1,triggered_zpripple_all{:});
triggered_MUA_zFR_cat = cat(1,triggered_MUA_zFR_all{:});

%% event-triggered ripple power 
axes(panels(1))
cla
hold on
t = linspace(-win1, win1, size(triggered_zpripple_cat,2));
shadedErrorBar(t, triggered_zpripple_cat, {@nanmean,@nansem});
xlabel('Time from replay (s)')
ylabel('Ripple power (z)')

%% event-triggered MUA
axes(panels(2))
cla
hold on
t = linspace(-win1, win1, size(triggered_MUA_zFR_cat,2));
shadedErrorBar(t, triggered_MUA_zFR_cat, {@nanmean,@nansem});
xlabel('Time from replay (s)')
ylabel('MUA firing rate (z)')


%%
bins_edges = linspace(-win1,win1,xcorr_hist_num_bins);

%% xcorr - posterior vs ripples
axes(panels(4))
cla
hold on
xlabel('Time from replay (s)')
ylabel('Counts')
histogram(cat(1,xcorr_posterior_ripples_tdiff{:}), bins_edges);
prc = [xcorr_posterior_ripples_fraction{:}] .* 100;
disp('posterior vs ripples:');
fprintf('%.2g %s %.2g %% (mean+-SD)\n',mean(prc),'\pm',std(prc));
fprintf('%.2g %s %.2g %% (mean+-SEM)\n',mean(prc),'\pm',nansem(prc));
% fraction_str = sprintf('%.2g %s %.2g %%',mean(prc),'\pm',std(prc));
% text(1,0.9,fraction_str,'Units','normalized','HorizontalAlignment','right','FontSize',8,'Interpreter','tex');
h=title('Replay events vs. ripple events','Units','normalized','Position',[0.5 1.05]);

%% xcorr - posterior vs MUA
axes(panels(5))
cla
hold on
xlabel('Time from replay (s)')
ylabel('Counts')
histogram(cat(1,xcorr_posterior_MUA_tdiff{:}), bins_edges);
prc = [xcorr_posterior_MUA_fraction{:}] .* 100;
disp('posterior vs MUA:');
fprintf('%.2g %s %.2g %% (mean+-SD)\n',mean(prc),'\pm',std(prc));
fprintf('%.2g %s %.2g %% (mean+-SEM)\n',mean(prc),'\pm',nansem(prc));
% fraction_str = sprintf('%.2g %s %.2g %%',mean(prc),'\pm',std(prc));
% text(1,0.9,fraction_str,'Units','normalized','HorizontalAlignment','right','FontSize',8,'Interpreter','tex');
h=title('Replay events vs. MUA events','Units','normalized','Position',[0.5 1.05]);

%% xcorr - posterior vs ripples
axes(panels(6))
cla
hold on
xlabel('Time from replay (s)')
ylabel('Counts')
histogram(cat(1,xcorr_posterior_PE_tdiff{:}), bins_edges);
prc = [xcorr_posterior_PE_fraction{:}] .* 100;
disp('posterior vs PE:');
fprintf('%.2g %s %.2g %% (mean+-SD)\n',mean(prc),'\pm',std(prc));
fprintf('%.2g %s %.2g %% (mean+-SEM)\n',mean(prc),'\pm',nansem(prc));
% fraction_str = sprintf('%.2g %s %.2g %%',mean(prc),'\pm',std(prc));
% text(1,0.9,fraction_str,'Units','normalized','HorizontalAlignment','right','FontSize',8,'Interpreter','tex');
h=title('Replay events vs. population events','Units','normalized','Position',[0.5 1.05]);

%% num ripples per replay event vs duration
axes(panels(7))
x = [seqs_all.duration];
y = [n_ripples_per_replay{:}];
cla
hold on
% plot(x,y,'o')
violinplot(x,y);
ylabel('Replay duration (s)');
xlabel('No. of ripples');

%% num ripples per replay event vs distance
axes(panels(8))
x = [seqs_all.distance];
y = [n_ripples_per_replay{:}];
cla
hold on
% % plot(x,y,'o')
violinplot(x,y);
ylabel('Distance (m)');
xlabel('No. of ripples');

%% num ripples per replay event vs compression
axes(panels(9))
x = [seqs_all.compression];
y = [n_ripples_per_replay{:}];
cla
hold on
% plot(x,y,'o');
violinplot(x,y);
ylabel('Compression');
xlabel('No. of ripples');

%% num ripples per replay event vs score
axes(panels(10))
x = [seqs_all.score];
y = [n_ripples_per_replay{:}];
cla
hold on
% plot(x,y,'o')
violinplot(x,y);
ylabel('Score');
xlabel('No. of ripples');

%%
axes(panels(1));
text(-0.3,1.1, 'A', 'Units','normalized','FontWeight','bold');
axes(panels(4));
text(-0.3,1.12, 'B', 'Units','normalized','FontWeight','bold');

%% save fig
if exist('bats_to_include','var')
    bats_str = ['_bats_' char(strjoin(""+bats,'_'))];
else
    bats_str = '_bats_all';
end

fig_name_out = fullfile(res_dir, [fig_name_str bats_str]);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
disp('figure was successfully saved to pdf/tiff/fig formats');
diary off




%%







