%% plot median/std of neural raw channels for bats across days

%%
clear;clc

%% load exp summary and choose exps
exp_t = DS_get_exp_summary();
bat_num = 148;
exp_t(~contains(exp_t.recordingArena, '200m'),:) = [];
exp_t(exp_t.position_data_exist==0,:) = [];
exp_t(exp_t.neural_data_exist==0,:) = [];
exp_t(~ismember(exp_t.batNum, [bat_num] ),:) = [];
exp_t

exps = cellfun(@(x)(exp_load_data(x,'details','csc_raw_stats')), exp_t.exp_ID, 'UniformOutput', false);
exps = exps(cellfun(@(x)(isfield(x,'csc_raw_stats')),exps));
exps = [exps{:}];

%% arrange data
dates = cellfun(@(x)(x.date), {exps.details});
stats = [exps.csc_raw_stats];
abs_median = cat(3,stats.csc_abs_median);
reref_abs_median = cat(3,stats.csc_reref_abs_median);

%% plot to figure
figure('Units','normalized','Position',[0 0 1 1])
pnl = panel();
pnl.pack('v',2,'h',4);
pnl.margin=25;
h=pnl.title(sprintf('Raw (filtered) neural data noise levels, bat %04d',bat_num));
h.FontSize = 20;
h.Position = [0.5 1.05];
h=pnl(1).title('No re-ref');
h.FontSize = 16;
for TT=1:4
    pnl(1,TT).select();
    title(sprintf('TT %d',TT))
    plot(dates, squeeze(abs_median(TT,:,:)),'.-', 'LineWidth',2);
    xlabel('Date')
    ylabel('median(abs(raw data)) [uVolt]')
    ylim([0 10])
    legend({'ch1','ch2','ch3','ch4'})
end
h=pnl(2).title('with re-ref');
h.FontSize = 16;
for TT=1:4
    pnl(2,TT).select();
    title(sprintf('TT %d',TT))
    plot(dates, squeeze(reref_abs_median (TT,:,:)),'.-', 'LineWidth',2);
    xlabel('Date')
    ylabel('median(abs(raw data)) [uVolt]')
    ylim([0 10])
    legend({'ch1','ch2','ch3','ch4'})
end

%% save figure
dir_OUT = 'L:\Analysis\Results\pre_proc';
if ~exist(dir_OUT,'dir') ;mkdir(dir_OUT); end
fig_filename = fullfile(dir_OUT, sprintf('csc_raw_stats_bat_%04d',bat_num));
saveas(gcf,fig_filename, 'tif')
saveas(gcf,fig_filename, 'fig')



%%
