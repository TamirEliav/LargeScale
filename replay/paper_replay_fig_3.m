%% Replay - Fig 2 - replay population 
%%
clear 
clc
close all

%% plotting options
behavior_ex_opt = 5

%% graphics params
% timediff_max = inf;
timediff_max = 100;

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Fig_3';
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
fig = figure;
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
% set(gcf, 'color', 'none');
set(groot, 'defaultAxesColor','None')
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

% create panels
panels{1}(1) = axes('position', [3 20 4 3]);
panels{2}(1) = axes('position', [8 20 5 3]);
panels{3}(1) = axes('position', [3 15 3 3]);
panels{4}(1) = axes('position', [8 15 3 3]);
panels{5}(1) = axes('position', [13 15 3 3]);
% panels{5}(2) = axes('position', [13 10 3 3]);

%% panels A - experimental setup
axes(panels{1}(1))
hold on
tunnel_2bats_image = imread('E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\resources\tunnel_2bats.jpg');
imshow(tunnel_2bats_image);

%% panels B - Behavior example
axes(panels{2}(1))
% figure
hold on
switch behavior_ex_opt
    case 1
        exp_ID = 'b2299_d191205';
        ti_seconds_in_session = [5200 5600];
    case 2
        exp_ID = 'b2299_d191213';
        ti_seconds_in_session = [5850 6200];
    case 3
        exp_ID = 'b2299_d191213';
        ti_seconds_in_session = [2900 3600];
    case 4
        exp_ID = 'b2299_d191203';
        ti_seconds_in_session = [500 900];
    case 5
        exp_ID = 'b2299_d191203';
        ti_seconds_in_session = [2700 3200];
end
exp = exp_load_data(exp_ID, 'details', 'pos');
epoch_type = 'rest';
params_opt = 11;
event_type = 'posterior';
[events, params]= decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
[seqs, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
events(~TF)=[];
session_ti = exp_get_sessions_ti(exp_ID,'Behave');
t0 = session_ti(1);
ti = ti_seconds_in_session.*1e6+t0;
lw = 2;
plot(exp.pos.proc_1D.ts, interp_nans(exp.pos.proc_1D.pos),'LineWidth',lw);
plot(exp.pos.proc_1D.other.ts, interp_nans(exp.pos.proc_1D.other.pos(1,:)),'LineWidth',1);
plot(exp.pos.proc_1D.co.ts,exp.pos.proc_1D.co.pos,'xk','MarkerSize',8)
plot([seqs.start_ts;seqs.end_ts],[seqs.start_pos; seqs.end_pos],'-m','LineWidth',1.3);
% plot([seqs.start_ts],[seqs.start_pos],'.m','MarkerSize',10)
plot([seqs([seqs.direction]== 1).end_ts],[seqs([seqs.direction]== 1).end_pos],'^m','MarkerSize',2,'MarkerFaceColor','m');
plot([seqs([seqs.direction]==-1).end_ts],[seqs([seqs.direction]==-1).end_pos],'vm','MarkerSize',2,'MarkerFaceColor','m');
xlim(ti)
rescale_plot_data('x',[1e-6/60 ti(1)]);
ylim([0 135])
% xticks(linspace(0,135,4))
yticks(linspace(0,135,4))
xlabel('Time (min)', 'Units','normalized', 'Position',[0.5 -0.16]);
ylabel('Position (m)', 'Units','normalized', 'Position',[-0.13 .5]);
hax=gca;
hax.XRuler.TickLabelGapOffset = -.5;
hax.YRuler.TickLabelGapOffset = 1;

%% legend
if exist('panels_behavior_legend','var')
    delete(panels_behavior_legend);
end
panels_behavior_legend = axes('position', [8 23.35 5 0.4]);
cla
hold on
plot([0 0.1],[0.8 0.8],'-','LineWidth',lw);
plot([0 0.1],[0.0 0.0],'-','LineWidth',1);
plot(0.6,0.8,'xk','MarkerSize',8);
plot([0.5 0.6]+0.05,[0 0],'-m','LineWidth',lw);
text(.15, .8, 'Recorded bat','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
text(.15, .0, 'Other bat','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
text(.7, .8, 'Cross-overs','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
text(.7, .0, 'Replay','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
xlim([0 1])
ylim([0 1])
axis off

%% panels C - scatter plot
axes(panels{3}(1))
hold on
data = load('F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & replay distance_gt_10.mat');
scatter(data.x(data.TF),data.y(data.TF),5,'k','filled');
axis equal
xlim([0 135])
ylim([0 135])
xticks(linspace(0,135,4))
yticks(linspace(0,135,4))
xlabel('Previous cross-over position (m)', 'Units','normalized', 'Position',[0.5 -0.16]);
ylabel('Replay position (m)', 'Units','normalized', 'Position',[-0.2 .5]);
% text(0.8,0.3,"{\itn} = "+ data.stats.n,'Units','normalized','FontSize',9);
% text(0.6,0.2,"{\itr} = "+ sprintf('%.2g',data.stats.Pearson.r), 'Units','normalized','FontSize',7);
% text(0.6,0.1,"{\itP} = "+ sprintf('%.2g',data.stats.Pearson.p), 'Units','normalized','FontSize',7);
text(.3,1.2,"{\it\rho} = "+ sprintf('%.2f',data.stats.Spearman.r), 'Units','normalized','FontSize',7);
text(.3,1.1,"{\itP} = "+ sprintf('%.2f',data.stats.Spearman.p), 'Units','normalized','FontSize',7);
h=refline(1,0);
h.Color = .8.*[1 1 1];
hax=gca;
hax.XRuler.TickLength(1) = 0.035;
hax.YRuler.TickLength(1) = 0.024;
hax.XRuler.TickLabelGapOffset = -.5;
hax.YRuler.TickLabelGapOffset = 1;

%% panels D - scatter plot (control - next crossover)
axes(panels{4}(1))
cla
hold on
data = load('F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & replay distance_gt_10.mat');
X = [data.seqs_all.next_co_pos];
Y = data.y;
X = X(data.TF)';
Y = Y(data.TF)';
[stats.Pearson.r, stats.Pearson.p] = corr(X,Y,'type','Pearson',rows='pairwise',tail='right');
[stats.Spearman.r, stats.Spearman.p] = corr(X,Y,'type','Spearman',rows='pairwise',tail='right');
scatter(X,Y,5,'k','filled');
axis equal
xlim([0 135])
ylim([0 135])
xticks(linspace(0,135,4))
yticks(linspace(0,135,4))
xlabel('Next cross-over position (m)', 'Units','normalized', 'Position',[0.5 -0.16]);
ylabel('Replay position (m)', 'Units','normalized', 'Position',[-0.2 .5]);
% text(0.8,0.3,"{\itn} = "+ data.stats.n,'Units','normalized','FontSize',9);
% text(0.8,0.2,"{\itr} = "+ sprintf('%.2g',stats.Pearson.r), 'Units','normalized','FontSize',7);
% text(0.8,0.1,"{\itP} = "+ sprintf('%.2g',stats.Pearson.p), 'Units','normalized','FontSize',7);
text(.3,1.2,"{\it\rho} = "+ sprintf('%.2f',stats.Spearman.r), 'Units','normalized','FontSize',7);
text(.3,1.1,"{\itP} = "+ sprintf('%.2f',stats.Spearman.p), 'Units','normalized','FontSize',7);
h=refline(1,0);
h.Color = .8.*[1 1 1];
hax=gca;
hax.XRuler.TickLength(1) = 0.035;
hax.YRuler.TickLength(1) = 0.024;
hax.XRuler.TickLabelGapOffset = -.5;
hax.YRuler.TickLabelGapOffset = 1;

%% panels E - time diff vs pos diff
axes(panels{5}(1))
cla
hold on
data = load('F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & replay distance_gt_10.mat');
xxx = [data.seqs_all.prev_co_time_diff];
yyy = abs(data.x-data.y);
ccc = [data.seqs_all.prev_co_same_map];
valid = data.TF & xxx < timediff_max;
xxx = xxx(valid)';
yyy = yyy(valid)';
ccc = ccc(valid)';
[stats.Pearson.r, stats.Pearson.p] = corr(xxx,yyy,'type','Pearson',rows='pairwise',tail='right');
[stats.Spearman.r, stats.Spearman.p] = corr(xxx,yyy,'type','Spearman',rows='pairwise',tail='right');
% scatter(xxx,yyy,10,ccc,'filled');
scatter(xxx,yyy,5,'k','filled');
xlabel('{\Delta} Time (s)', 'Units','normalized', 'Position',[0.5 -0.16]);
ylabel('{\Delta} Position (m)', 'Units','normalized', 'Position',[-0.2 .5]);
% text(0.8,0.3,"{\itn} = "+ data.stats.n,'Units','normalized','FontSize',9);
% text(1.1,.2,"{\itr} = "+ sprintf('%.2g',stats.Pearson.r), 'Units','normalized','FontSize',7);
% text(1.1,.1,"{\itP} = "+ sprintf('%.2g',stats.Pearson.p), 'Units','normalized','FontSize',7);
text(.3,1.2,"{\it\rho} = "+ sprintf('%.2g',stats.Spearman.r), 'Units','normalized','FontSize',7);
text(.3,1.1,"{\itP} = "+ sprintf('%.2g',stats.Spearman.p), 'Units','normalized','FontSize',7);
hax=gca;
hax.XRuler.TickLength(1) = 0.035;
hax.YRuler.TickLength(1) = 0.024;
hax.XRuler.TickLabelGapOffset = -.5;
hax.YRuler.TickLabelGapOffset = 1;

%% panels E - time diff vs pos diff (box plot option)
if length(panels{5}) == 2
axes(panels{5}(2))
cla
hold on
data = load('F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & replay distance_gt_10.mat');
xxx = [data.seqs_all.prev_co_time_diff];
yyy = abs(data.x-data.y);
ccc = [data.seqs_all.prev_co_same_map];
valid = data.TF & xxx < timediff_max;
xxx = xxx(valid)';
yyy = yyy(valid)';
ccc = ccc(valid)';

thr = 25;
IX1 = xxx<thr;
IX2 = xxx>=thr;
G = zeros(size(xxx));
G(IX1) = 1;
G(IX2) = 2;
groupcounts(G)
boxplot(yyy,G)
pval = ranksum(yyy(IX1),yyy(IX2),"tail","left");
text(0.5,1.1,"{\itP} = "+ sprintf('%.2g',pval), 'Units','normalized','FontSize',7,'HorizontalAlignment','center');
hax=gca;
hax.XTick = [1 2];
hax.XTickLabel = {"< "+thr,">= "+thr};
hax.XLim = [.5 2.5];
box off
xlabel('{\Delta} Time (s)', 'Units','normalized', 'Position',[0.5 -0.16]);
ylabel('{\Delta} Position (m)', 'Units','normalized', 'Position',[-0.2 .5]);
hax.XRuler.TickLength(1) = 0.025;
hax.YRuler.TickLength(1) = 0.024;
hax.XRuler.TickLabelGapOffset = -.5;
hax.YRuler.TickLabelGapOffset = 1;
end

%% add panel letters
font_size = 11;
axes(panels{1}(1))
text(-0.2,1.12, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{2}(1))
text(-0.2,1.12, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{3}(1))
text(-0.35,1.12, 'c', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{4}(1))
text(-0.35,1.12, 'd', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{5}(1))
text(-0.35,1.12, 'e', 'Units','normalized','FontWeight','bold','FontSize',font_size);


%%
fig_name = sprintf('%s_panel_B_opt_%d',fig_name_str,behavior_ex_opt);
fig_name + "_" + timediff_max 
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%
