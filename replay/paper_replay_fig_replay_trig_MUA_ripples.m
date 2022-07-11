%% Replay figure: Replay vs ripples/MUA/PE

%%
close all
clear 
clc

%% options
% bats_to_include = [34 9861];
win = 2; % in seconds

%% define output files
res_dir =  'L:\paper_replay\figures\';
mkdir(res_dir)
fig_name_str = 'Fig_replay_trig_MUA_ripples';
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
panels(1) = axes('position', [4 15 panels_size]);
panels(2) = axes('position', [12 15 panels_size]);

%% arrange data
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
clear exp_list
groupsummary(T,'bat_num')
if exist('bats_to_include','var')
    T = groupfilter(T,"bat_num",@(x)ismember(x,bats_to_include),'bat_num');
end
bats = unique(T.bat_num)

%% arrange data
epoch_type = 'sleep'; params_opt = 11;
data_sleep = arrange_data(T, epoch_type, params_opt, win);
epoch_type = 'rest'; params_opt = 11;
data_rest = arrange_data(T, epoch_type, params_opt, win);

%% event-triggered ripple power 
axes(panels(1))
cla
hold on
t = linspace(-win, win, size(data_sleep.triggered_zpripple_cat,2));
shadedErrorBar(t, data_sleep.triggered_zpripple_cat, {@nanmean,@nansem},'lineprops',{'r-'});
shadedErrorBar(t, data_rest.triggered_zpripple_cat, {@nanmean,@nansem},'lineprops',{'b-'});
xlabel('Time from replay (s)')
ylabel('Ripple power (z)')

%% event-triggered MUA
axes(panels(2))
cla
hold on
t = linspace(-win, win, size(data_sleep.triggered_MUA_zFR_cat,2));
shadedErrorBar(t, data_sleep.triggered_MUA_zFR_cat, {@nanmean,@nansem},'lineprops',{'r-'});
shadedErrorBar(t, data_rest.triggered_MUA_zFR_cat, {@nanmean,@nansem},'lineprops',{'b-'});
xlabel('Time from replay (s)')
ylabel('MUA firing rate (z)')

%%
% axes(panel_legend)
panel_legend = axes('position', [9 21 0.3 0.3]);
cla
hold on
plot([0 1],[1 1],'-r','LineWidth',2)
plot([0 1],[0 0],'-b','LineWidth',2)
text(1.2,1, "Sleep",FontSize=7)
text(1.2,0, "Rest",FontSize=7)
axis off

%%
% %%
% bins_edges = linspace(-win,win,xcorr_hist_num_bins);
% 
% %% xcorr - posterior vs ripples
% axes(panels(4))
% cla
% hold on
% xlabel('Time from replay (s)')
% ylabel('Counts')
% histogram(cat(1,xcorr_posterior_ripples_tdiff{:}), bins_edges);
% prc = [xcorr_posterior_ripples_fraction{:}] .* 100;
% disp('posterior vs ripples:');
% fprintf('%.2g %s %.2g %% (mean+-SD)\n',mean(prc),'\pm',std(prc));
% fprintf('%.2g %s %.2g %% (mean+-SEM)\n',mean(prc),'\pm',nansem(prc));
% % fraction_str = sprintf('%.2g %s %.2g %%',mean(prc),'\pm',std(prc));
% % text(1,0.9,fraction_str,'Units','normalized','HorizontalAlignment','right','FontSize',8,'Interpreter','tex');
% h=title('Replay events vs. ripple events','Units','normalized','Position',[0.5 1.05]);
% 
% %% xcorr - posterior vs MUA
% axes(panels(5))
% cla
% hold on
% xlabel('Time from replay (s)')
% ylabel('Counts')
% histogram(cat(1,xcorr_posterior_MUA_tdiff{:}), bins_edges);
% prc = [xcorr_posterior_MUA_fraction{:}] .* 100;
% disp('posterior vs MUA:');
% fprintf('%.2g %s %.2g %% (mean+-SD)\n',mean(prc),'\pm',std(prc));
% fprintf('%.2g %s %.2g %% (mean+-SEM)\n',mean(prc),'\pm',nansem(prc));
% % fraction_str = sprintf('%.2g %s %.2g %%',mean(prc),'\pm',std(prc));
% % text(1,0.9,fraction_str,'Units','normalized','HorizontalAlignment','right','FontSize',8,'Interpreter','tex');
% h=title('Replay events vs. MUA events','Units','normalized','Position',[0.5 1.05]);
% 
% %% xcorr - posterior vs ripples
% axes(panels(6))
% cla
% hold on
% xlabel('Time from replay (s)')
% ylabel('Counts')
% histogram(cat(1,xcorr_posterior_PE_tdiff{:}), bins_edges);
% prc = [xcorr_posterior_PE_fraction{:}] .* 100;
% disp('posterior vs PE:');
% fprintf('%.2g %s %.2g %% (mean+-SD)\n',mean(prc),'\pm',std(prc));
% fprintf('%.2g %s %.2g %% (mean+-SEM)\n',mean(prc),'\pm',nansem(prc));
% % fraction_str = sprintf('%.2g %s %.2g %%',mean(prc),'\pm',std(prc));
% % text(1,0.9,fraction_str,'Units','normalized','HorizontalAlignment','right','FontSize',8,'Interpreter','tex');
% h=title('Replay events vs. population events','Units','normalized','Position',[0.5 1.05]);
% 
% %% num ripples per replay event vs duration
% axes(panels(7))
% x = [seqs_all.duration];
% y = [n_ripples_per_replay{:}];
% cla
% hold on
% % plot(x,y,'o')
% violinplot(x,y);
% ylabel('Replay duration (s)');
% xlabel('No. of ripples');
% 
% %% num ripples per replay event vs distance
% axes(panels(8))
% x = [seqs_all.distance];
% y = [n_ripples_per_replay{:}];
% cla
% hold on
% % % plot(x,y,'o')
% violinplot(x,y);
% ylabel('Distance (m)');
% xlabel('No. of ripples');
% 
% %% num ripples per replay event vs compression
% axes(panels(9))
% x = [seqs_all.compression];
% y = [n_ripples_per_replay{:}];
% cla
% hold on
% % plot(x,y,'o');
% violinplot(x,y);
% ylabel('Compression');
% xlabel('No. of ripples');
% 
% %% num ripples per replay event vs score
% axes(panels(10))
% x = [seqs_all.score];
% y = [n_ripples_per_replay{:}];
% cla
% hold on
% % plot(x,y,'o')
% violinplot(x,y);
% ylabel('Score');
% xlabel('No. of ripples');

% %%
% axes(panels(1));
% text(-0.3,1.1, 'A', 'Units','normalized','FontWeight','bold');
% axes(panels(4));
% text(-0.3,1.12, 'B', 'Units','normalized','FontWeight','bold');

%% save fig
if exist('bats_to_include','var')
    bats_str = ['bats_' char(strjoin(""+bats,'_'))];
else
    bats_str = 'bats_all';
end

fig_name_out = fullfile(res_dir, sprintf('%s_%s_%s',fig_name_str, epoch_type, bats_str));
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
disp('figure was successfully saved to pdf/tiff/fig formats');
diary off




%%
function data = arrange_data(T, epoch_type, params_opt, win)

%%
% epoch_type = 'sleep';
% epoch_type = 'rest';
% params_opt = 11;
% win = 2; % in seconds
% xcorr_hist_num_bins = 25;

%%
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
    [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    if isempty(events)
        continue;
    end
    seqs = [events.seq_model];

    %% apply inclusion criteria
    [seqs, TF] = decoding_apply_seq_inclusion_criteria(seqs);
    if isempty(seqs)
        continue;
    end
    [seqs.ts] = disperse(mean([seqs.start_ts; seqs.end_ts]));
    seqs_all = [seqs_all seqs];

    %% trigger signals (ripples/MUA) around posterior events
    win_samples = round(win*exp.ripples.fs);
    trigger_IX = interp1(exp.ripples.t, 1:length(exp.ripples.t), [seqs.ts], 'nearest');
    triggered_zpripple_all{ii_exp} = trigger_signal_by_IX(exp.ripples.zpripple', trigger_IX, win_samples);
    win_samples = round(win*exp.MUA.fs);
    trigger_IX = interp1(exp.MUA.t, 1:length(exp.MUA.t), [seqs.ts], 'nearest');
    triggered_MUA_zFR_all{ii_exp} = trigger_signal_by_IX(exp.MUA.zFR, trigger_IX, win_samples);

    %% events xcorr (posterior vs ripples)
    t1 = [seqs.ts];
    t2 = [exp.ripples.events.peak_ts];
    tdiff = (t2'-t1).*1e-6;
    fraction = sum(any(abs(tdiff)<win)) / size(tdiff,2); % fraction of the posterior events which had ripple event nearby
    tdiff = tdiff(:);
    tdiff(abs(tdiff)>win) = [];
    xcorr_posterior_ripples_tdiff{ii_exp} = tdiff;
    xcorr_posterior_ripples_fraction{ii_exp} = fraction;

    %% events xcorr (posterior vs MUA)
    t1 = [seqs.ts];
    t2 = [exp.MUA.events.peak_ts];
    tdiff = (t2'-t1).*1e-6;
    fraction = sum(any(abs(tdiff)<win)) / size(tdiff,2); % fraction of the posterior events which had MUA event nearby
    tdiff = tdiff(:);
    tdiff(abs(tdiff)>win) = [];
    xcorr_posterior_MUA_tdiff{ii_exp} = tdiff;
    xcorr_posterior_MUA_fraction{ii_exp} = fraction;

    %% events xcorr (posterior vs PE)
    t1 = [seqs.ts];
    t2 = [exp.PE.thr.peak_ts];
    tdiff = (t2'-t1).*1e-6;
    fraction = sum(any(abs(tdiff)<win)) / size(tdiff,2); % fraction of the posterior events which had population event nearby
    tdiff = tdiff(:);
    tdiff(abs(tdiff)>win) = [];
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

%% arrange output in struct
data.triggered_zpripple_cat = triggered_zpripple_cat;
data.triggered_MUA_zFR_cat = triggered_MUA_zFR_cat;

end



