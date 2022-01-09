%% plot median/std of neural raw channels for bats across days

%%
clear;clc

%% load exp summary and choose exps
exp_t = DS_get_exp_summary();
bat_num = 2382;
% exp_t(~contains(exp_t.recordingArena, '200m'),:) = [];
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
numTT = exps(1).details.numTT;

%% plot to figure (absolute values)
fig=figure;
fig.WindowState = 'maximized';
nrows = sqrt(numTT);
ncols = sqrt(numTT);
tiledlayout(nrows,ncols);
h=sgtitle(sprintf('Raw (filtered in spikes freq) neural data noise levels, bat %04d',bat_num));
for TT=1:numTT
    nexttile
    hold on
    title(sprintf('TT %d',TT))
    plot(dates, squeeze(abs_median(TT,:,:)),'o-', 'LineWidth',2);
    plot(dates, squeeze(reref_abs_median(TT,:,:)),'o--', 'LineWidth',2);
    title("TT"+TT)
    xlabel('Date')
    ylabel({'median(abs(data))';'[uVolt]'})
end
hl=legend(["ch"+[1:4] "ch"+[1:4]+" (re-ref)"],'Location','bestoutside');
dir_OUT = 'L:\Analysis\Results\pre_proc';
mkdir(dir_OUT);
fig_filename = fullfile(dir_OUT, sprintf('csc_raw_stats_bat_%04d',bat_num));
saveas(fig,fig_filename, 'tif')
saveas(fig,fig_filename, 'fig')

%% plot to figure (relative values)
fig=figure;
fig.WindowState = 'maximized';
nrows = sqrt(numTT);
ncols = sqrt(numTT);
tiledlayout(nrows,ncols);
h=sgtitle(sprintf('Raw (filtered in spikes freq) neural data relative noise levels, bat %04d',bat_num));
for TT=1:numTT
    nexttile
    hold on
    title(sprintf('TT %d',TT));
    v1 = squeeze(abs_median(TT,:,:));
    v2 = squeeze(reref_abs_median(TT,:,:));
    v = v1-v2; % if we get positive values it means that re-reference worked!
    plot(dates, v, 'o-', 'LineWidth',2);
    yline(0)
    title("TT"+TT)
    xlabel('Date')
    ylabel({'median(|(data)|)-median(|(data reref)|)';'[uVolt]'})
end
hl=legend(["ch"+[1:4]],'Location','bestoutside');
dir_OUT = 'L:\Analysis\Results\pre_proc';
mkdir(dir_OUT);
fig_filename = fullfile(dir_OUT, sprintf('csc_raw_stats_bat_%04d_relative',bat_num));
saveas(fig,fig_filename, 'tif')
saveas(fig,fig_filename, 'fig')


%%
