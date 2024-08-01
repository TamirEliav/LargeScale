%% Replay - Fig 3 - 2 bats crossovers
%%
clear 
clc
close all

%% plotting options
behavior_ex_opt = 2 % number 2 was chosen

replay_ex_opt = 1;
replay_examples_list = {
    {'rest',11,'b2299_d191213',20}
    };
replay_examples_list = cellfun(@(c)cell2struct(c,{'epoch_type','params_opt','exp_ID','event_num'},2), replay_examples_list)

%% graphics params
% timediff_max = inf;
timediff_max = 100;
epoch_types = {'sleep','rest'};
exp_types = {'1 bat','2 bats'};
panels_xlim = [-0.0950    1.9950; -1.9500   40.9500; -0.4500   56.5500; -0.0210    0.4410];
panels_xticks = {[0 0.5 1 1.5],[0 20 40],[0:10:50],[0 0.2 0.4]};
smooth_method = 'movmean';
smooth_n_points = 19;

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
panels{1}(1) = axes('position', [2.7 20.3 4 3]);
panels{2}(1) = axes('position', [8 20 5 3]);
panels{3}(1) = axes('position', [14.8 20 3 3]);
panels{3}(2) = axes('position', [14.8 23 3 .5]);
x = linspace(3,10,4);
x = x+[0 0.15 -0.1 0.1];
% panels{4}(1) = axes('position', [x(1) 15 1 3]);
% panels{4}(2) = axes('position', [x(2) 15 1 3]);
% panels{4}(3) = axes('position', [x(3) 15 1 3]);
% panels{4}(4) = axes('position', [x(4) 15 1 3]);
panels{5}(1,1) = axes('position', [3 15 3 3]);
% panels{5}(1,2) = axes('position', [3 12 3 1.5]);
% panels{5}(1,3) = axes('position', [3  7 3 3]);
% panels{5}(1,4) = axes('position', [3  3 3 3]);
panels{5}(2,1) = axes('position', [7.5 15 3 3]);
% panels{5}(2,2) = axes('position', [7.5 12 3 1.5]);
% panels{5}(2,3) = axes('position', [7.5  7 3 3]);
% panels{5}(2,4) = axes('position', [7.5  3 3 3]);
% panels{6}(1) = axes('position', [16 6 3 2]);
% panels{6}(2) = axes('position', [16 3 3 2]);

%% panels A - experimental setup
axes(panels{1}(1))
hold on
tunnel_2bats_image = imread('E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\resources\tunnel_2bats.jpg');
imshow(tunnel_2bats_image);

% delete(panels{1}(2))
panels{1}(2) = axes('position', [2.7 22.9 .2 .3]);
axes(panels{1}(2))
cla reset
hold on
axis off
plot([0 1],[1 1],'LineWidth',2)
plot([0 1],[0 0],'LineWidth',1)
text(1.5,1,'Recorded bat','Units','normalized','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
text(1.5,0,'Other bat','Units','normalized','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
xlim([0 1])
ylim([0 1])

%% panels B - Behavior example
axes(panels{2}(1))
cla
hold on
switch behavior_ex_opt
    case 1
        exp_ID = 'b2299_d191205';
        ti_seconds_in_session = [5200 5600];
    case 2
        exp_ID = 'b2299_d191213';
        ti_seconds_in_session = 5850 +[0 4*60]+[0.4 -0.1]; %[5850 6200];
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
plot([seqs.start_ts;seqs.end_ts],[seqs.start_pos; seqs.end_pos],'-k','LineWidth',1.);
plot([seqs.start_ts],[seqs.start_pos],'.k','MarkerSize',7)
% plot([seqs([seqs.direction]== 1).end_ts],[seqs([seqs.direction]== 1).end_pos],'^m','MarkerSize',2,'MarkerFaceColor','m');
% plot([seqs([seqs.direction]==-1).end_ts],[seqs([seqs.direction]==-1).end_pos],'vm','MarkerSize',2,'MarkerFaceColor','m');
xline(35969648929.0)
xlim(ti)
rescale_plot_data('x',[1e-6/60 ti(1)]);
ylim([0 135])
% xticks(linspace(0,135,4))
yticks(linspace(0,135,4))
xlabel('Time (min)', 'Units','normalized', 'Position',[0.5 -0.11]);
ylabel('Position (m)', 'Units','normalized', 'Position',[-0.13 .5]);
hax=gca;
hax.XRuler.TickLabelGapOffset = -1.8;
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
plot([0.5 0.6]+0.05,[0 0],'-m','LineWidth',1.3);
plot(0.65,0,'>m','MarkerSize',2,'MarkerFaceColor','m');
text(.15, .8, 'Recorded bat','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
text(.15, .0, 'Other bat','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
text(.7, .8, 'Cross-overs','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
text(.7, .0, 'Replay','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
xlim([0 1])
ylim([0 1])
axis off

%% replay example
if 1
%% load data
win_s = 1;
cmap = bone;
cmap = flipud(cmap);
ex = replay_examples_list(replay_ex_opt);
addFieldsToWorkspace(ex);
decode = decoding_load_data(exp_ID, epoch_type, params_opt);
exp = exp_load_data(exp_ID,'details','path','MUA','ripples');
events = decoding_load_events_quantification(exp_ID,epoch_type,params_opt,"posterior");
event = events(event_num);
seq = event.seq_model;
seq_ti = [event.start_ts event.end_ts];
t0 = mean(seq_ti);
ti = t0 + [-1 1].*win_s*1e6;

TT = exp.ripples.stats.best_TT;
[LFP.signal, LFP.ts, LFP.fs, LFP.params] = LFP_load(exp_ID,TT,'band','ripple','limits_ts',ti);
LFP.avg_signal = nanmean(LFP.signal,[2 3]);

%% plot LFP
% axes(panels{3}(3));
% cla
% hold on
% plot(LFP.ts, LFP.avg_signal,'k');
% xlim(seq_ti+[-1 1].*0.2*range(seq_ti))
% xticks([])
% yticks([])
% rescale_plot_data('x',[1e-6 seq_ti(1)]);
% axis off
% title(sprintf('%s_%s_%d',epoch_type,exp_ID,event_num),'Interpreter','none','FontWeight','normal','FontSize',6);

%% plot posterior (state)
axes(panels{3}(2));
cla
hold on
IX = get_data_in_ti(decode.time, ti);
prob_t = decode.time(IX);
prob_state = squeeze(decode.posterior_state(event.state_num,IX));
plot(prob_t, prob_state, 'k','LineWidth',2);
hax = gca;
hax.XLim = prob_t([1 end]);
hax.YLim = [0 1];
box on
colormap(cmap);
hax.XLim = seq_ti+[-1 1].*0.2*range(seq_ti);
hax.TickDir = 'out';
hax.TickLength = [0.02 0.02];
hax.XRuler.TickLabelGapOffset = -4;
ylabel('Prob.','Units','normalized','Position',[-0.22 0.5]);
rescale_plot_data('x',[1e-6 seq_ti(1)]);
    
%% plot posterior (position)
axes(panels{3}(1));
cla
hold on
IX = get_data_in_ti(decode.time, ti);
prob_t = decode.time(IX);
prob_pos = squeeze(decode.posterior(:,event.state_num,IX));
imagesc(prob_t, decode.pos, prob_pos);
plot([seq.start_ts seq.end_ts],[seq.start_pos seq.end_pos],'-r','LineWidth',0.8);
hax = gca;
clim_prctiles = [1 99];
hax.CLim = prctile(prob_pos(:),clim_prctiles);
% hax.CLim = [0 max(prob_pos(:))];
hax.XLim = prob_t([1 end]);
hax.YLim = [min(decode.pos) max(decode.pos)] + [-1 1].*median(diff(decode.pos))*1;
box on
colormap(cmap);
hax.XLim = seq_ti+[-1 1].*0.2*range(seq_ti);
hax.TickDir = 'out';
hax.TickLength = [0.02 0.02];
hax.XRuler.TickLabelGapOffset = -2;
xlabel('Time (s)','Units','normalized','Position',[0.5 -0.105]);
ylabel('Position (m)','Units','normalized','Position',[-0.22 0.5]);
rescale_plot_data('x',[1e-6 seq_ti(1)]);

%% link x axes
linkaxes(panels{3},'x');

%% add colorbar
hcb = colorbar('east');
hcb.Units = 'centimeters';
cb_offset = 1.12;
hcb.Position(1) = panels{3}(1).Position(1) + panels{3}(1).Position(3)*cb_offset;
cb_offset_middle = (hcb.Position(1)+hcb.Position(3)/2-panels{3}(1).Position(1))/panels{3}(1).Position(3);
hcb.Label.Rotation = -90;
hcb.Label.Position(1) = 1.5;
hcb.Label.String = 'Probability';
hcb.Ticks = [];
text(cb_offset_middle,1,'Max','Units','normalized','FontSize',7,'HorizontalAlignment','center');
text(cb_offset_middle,0,'0','Units','normalized','FontSize',7,'HorizontalAlignment','center');
% text(cb_offset_middle,1,clim_prctiles(2)+"%",'Units','normalized','FontSize',7,'HorizontalAlignment','center');
% text(cb_offset_middle,0,clim_prctiles(1)+"%",'Units','normalized','FontSize',7,'HorizontalAlignment','center');
end







%% panels C+D - scatters (load data)
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & same map.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & same map & forward.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & same map & reverse.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & replay distance_gt_5.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & replay distance_gt_10.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & forward.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & reverse.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay pos between 25-115, dec acc_gt_65%.mat';
data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay pos between 25-115, dec acc_gt_70%.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay pos between 25-115 & same map, dec acc_gt_65%.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay pos between 25-115 & same map, dec acc_gt_70%.mat';

%% panels C - main scatter plot
axes(panels{5}(1,1))
cla reset
hold on
data = load(data_filename);
X = data.x(data.TF);
Y = data.y(data.TF);
scatter(X,Y,5,'k','filled');
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
text(.3,1.1,"{\itP} = "+ sprintf('%.2g',data.stats.Spearman.p), 'Units','normalized','FontSize',7);
% text(.3,1.3,"{\itn} = "+ sprintf('%d',sum(data.TF)), 'Units','normalized','FontSize',7);
% text(0,-.4,data.msg_str, 'Units','normalized','FontSize',10);
h=refline(1,0);
h.Color = .8.*[1 1 1];
hax=gca;
hax.XRuler.TickLength(1) = 0.035;
hax.YRuler.TickLength(1) = 0.024;
hax.XRuler.TickLabelGapOffset = -.5;
hax.YRuler.TickLabelGapOffset = 1;

% axes(panels{5}(1,3))
% cla reset
% hold on
% [~,IX_sorted] = sort(X,'ascend');
% X2 = smoothdata(X(IX_sorted),smooth_method,smooth_n_points);
% Y2 = smoothdata(Y(IX_sorted),smooth_method,smooth_n_points);
% scatter(X,Y,5,'k','filled');
% plot(X2,Y2,'-m')
% axis equal
% xlim([0 135])
% ylim([0 135])
% xticks(linspace(0,135,4))
% yticks(linspace(0,135,4))
% 
% axes(panels{5}(1,2))
% cla
% hold on
% diff_data_prev = abs(X-Y);
% nShuffles = 10000;
% diff_shuffles_prev = zeros(1,nShuffles);
% rng(0);
% for ii_shuffle = 1:nShuffles
%     IX =randperm(length(X));    
%     diff_shuffles_prev(ii_shuffle) = mean(abs(X(IX)-Y));
% end
% diff_data_z = (mean(diff_data_prev)-mean(diff_shuffles_prev))/std(diff_shuffles_prev);
% pval = normcdf(diff_data_z,0,1);
% hh=histogram(diff_shuffles_prev);
% hh.FaceColor = .5*[1 1 1];
% xline(mean(diff_data_prev),'r');
% text(0.7,0.85,"{\itP} = "+sprintf('%.2f',pval),'units','normalized','FontSize',7);
% xlabel({'|\Deltaposition| (m)';'Replay vs. previous cross-over'});
% ylabel('Counts');
% 
% axes(panels{5}(1,4))
% cla reset
% hold on
% g = nan(size(X));
% g(X<45) = 1;
% g(X>=45 & X<=90) = 2;
% g(X>90) = 3;
% boxplot(Y,g)
% xticklabels({'<45','>45 & < 90','>90'})

%% panels D - scatter plot (control - next crossover)
axes(panels{5}(2,1))
cla reset
hold on
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

% axes(panels{5}(2,3))
% cla reset
% hold on
% [~,IX_sorted] = sort(X,'ascend');
% X2 = smoothdata(X(IX_sorted),smooth_method,smooth_n_points);
% Y2 = smoothdata(Y(IX_sorted),smooth_method,smooth_n_points);
% scatter(X,Y,5,'k','filled');
% plot(X2,Y2,'-m')
% axis equal
% xlim([0 135])
% ylim([0 135])
% xticks(linspace(0,135,4))
% yticks(linspace(0,135,4))

% axes(panels{5}(2,2))
% cla
% hold on
% diff_data_next = abs(X-Y);
% nShuffles = 10000;
% diff_shuffles_next = zeros(1,nShuffles);
% rng(0);
% for ii_shuffle = 1:nShuffles
%     IX =randperm(length(X));    
%     diff_shuffles_next(ii_shuffle) = mean(abs(X(IX)-Y));
% end
% diff_data_z = (mean(diff_data_next)-mean(diff_shuffles_next))/std(diff_shuffles_next);
% pval = normcdf(diff_data_z,0,1);
% hh=histogram(diff_shuffles_next);
% hh.FaceColor = .5*[1 1 1];
% xline(mean(diff_data_next),'r');
% text(0.7,0.85,"{\itP} = "+sprintf('%.2f',pval),'units','normalized','FontSize',7);
% xlabel({'|\Deltaposition| (m)';'Replay vs. next cross-over'});
% ylabel('Counts');
% 
% axes(panels{5}(2,4))
% cla reset
% hold on
% g = nan(size(X));
% g(X<45) = 1;
% g(X>=45 & X<=90) = 2;
% g(X>90) = 3;
% boxplot(Y,g)
% xticklabels({'<45','>45 & < 90','>90'})

%%
% axes(panels{6}(1))
% cla reset
% hold on
% hax=gca
% ecdf(diff_data_prev)
% ecdf(diff_data_next)
% hl=legend({'Next';'Previous'},'Location','none');
% hl.Units = 'centimeters';
% hl.Position = [hax.Position([1 2]) + [2 0.3] .1 .1];
% hl.Box = 'off';
% xlabel('\DeltaPosition')
% ylabel('CDF')
% %%
% axes(panels{6}(2))
% cla reset
% hold on
% hax=gca;
% hh = histogram(diff_data_prev,'DisplayStyle','stairs','DisplayName','Previous cross-over','EdgeColor','r','BinWidth',10);
% hh = histogram(diff_data_next,'DisplayStyle','stairs','DisplayName','Next cross-over','EdgeColor','b','BinWidth',10);
% hl=legend({'Next';'Previous'},'Location','none');
% hl.Units = 'centimeters';
% hl.Position = [hax.Position([1 2]) + [2.4 1.6] .1 .1];
% hl.Box = 'off';
% xlabel('\DeltaPosition')
% ylabel('Counts')

%% panels E - time diff vs pos diff (recency effect)
% axes(panels{6}(1))
% cla
% hold on
% xxx = [data.seqs_all.prev_co_time_diff];
% yyy = abs(data.x-data.y);
% ccc = [data.seqs_all.prev_co_same_map];
% valid = data.TF & xxx < timediff_max;
% xxx = xxx(valid)';
% yyy = yyy(valid)';
% ccc = ccc(valid)';
% [stats.Pearson.r, stats.Pearson.p] = corr(xxx,yyy,'type','Pearson',rows='pairwise',tail='right');
% [stats.Spearman.r, stats.Spearman.p] = corr(xxx,yyy,'type','Spearman',rows='pairwise',tail='right');
% % scatter(xxx,yyy,10,ccc,'filled');
% scatter(xxx,yyy,5,'k','filled');
% xlabel('{\Delta} Time (s)', 'Units','normalized', 'Position',[0.5 -0.16]);
% ylabel('{\Delta} Position (m)', 'Units','normalized', 'Position',[-0.2 .5]);
% % text(0.8,0.3,"{\itn} = "+ data.stats.n,'Units','normalized','FontSize',9);
% % text(1.1,.2,"{\itr} = "+ sprintf('%.2g',stats.Pearson.r), 'Units','normalized','FontSize',7);
% % text(1.1,.1,"{\itP} = "+ sprintf('%.2g',stats.Pearson.p), 'Units','normalized','FontSize',7);
% text(.3,1.2,"{\it\rho} = "+ sprintf('%.2g',stats.Spearman.r), 'Units','normalized','FontSize',7);
% text(.3,1.1,"{\itP} = "+ sprintf('%.2g',stats.Spearman.p), 'Units','normalized','FontSize',7);
% hax=gca;
% hax.XRuler.TickLength(1) = 0.035;
% hax.YRuler.TickLength(1) = 0.024;
% hax.XRuler.TickLabelGapOffset = -.5;
% hax.YRuler.TickLabelGapOffset = 1;

%% panels E - time diff vs pos diff (box plot option)
% if length(panels{6}) == 2
% axes(panels{6}(2))
% cla
% hold on
% xxx = [data.seqs_all.prev_co_time_diff];
% yyy = abs(data.x-data.y);
% ccc = [data.seqs_all.prev_co_same_map];
% valid = data.TF & xxx < timediff_max;
% xxx = xxx(valid)';
% yyy = yyy(valid)';
% ccc = ccc(valid)';
% 
% thr = 25;
% IX1 = xxx<thr;
% IX2 = xxx>=thr;
% G = zeros(size(xxx));
% G(IX1) = 1;
% G(IX2) = 2;
% groupcounts(G)
% boxplot(yyy,G)
% pval = ranksum(yyy(IX1),yyy(IX2),"tail","left");
% text(0.5,1.1,"{\itP} = "+ sprintf('%.2g',pval), 'Units','normalized','FontSize',7,'HorizontalAlignment','center');
% hax=gca;
% hax.XTick = [1 2];
% hax.XTickLabel = {"< "+thr,">= "+thr};
% hax.XLim = [.5 2.5];
% box off
% xlabel('{\Delta} Time (s)', 'Units','normalized', 'Position',[0.5 -0.16]);
% ylabel('{\Delta} Position (m)', 'Units','normalized', 'Position',[-0.2 .5]);
% hax.XRuler.TickLength(1) = 0.025;
% hax.YRuler.TickLength(1) = 0.024;
% hax.XRuler.TickLabelGapOffset = -.5;
% hax.YRuler.TickLabelGapOffset = 1;
% end

%% ------------------------------------------------------------------------
% % % % compare 1 bat vs 2 bats experiments
% % % if 0
% % % %% load data - 1 bat / 2 bats experiments
% % % epoch_types = {'sleep','rest'};
% % % exp_list_1bat = decoding_get_inclusion_list();
% % % exp_list_2bats = {
% % % %     'b2299_d191202', % <70% accuracy, excluded
% % % %     'b2299_d191203', % <70% accuracy, excluded
% % %     'b2299_d191204',
% % %     'b2299_d191205',
% % %     'b2299_d191208',
% % %     'b2299_d191209',
% % % %     'b2299_d191210', % <70% accuracy, excluded
% % %     'b2299_d191213',
% % %     };
% % % exp_lists = {exp_list_1bat,exp_list_2bats};
% % % 
% % % events_all_per_session = {};
% % % for ii_exp_type = 1:length(exp_types)
% % %     exp_list = exp_lists{ii_exp_type};
% % %     for ii_epoch_type = 1:length(epoch_types)
% % %         for ii_exp = 1:length(exp_list)
% % %             exp_ID = exp_list{ii_exp};
% % %             exp = exp_load_data(exp_ID,'details');
% % %             epoch_type = epoch_types{ii_epoch_type};
% % %             params_opt = 11;
% % %             event_type = 'posterior';
% % %             [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
% % %             [~, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
% % %             events(~TF) = [];
% % %             if isempty(events)
% % %                 continue;
% % %             end
% % %             [events.recordingArena] = deal(exp.details.recordingArena);
% % %             event_struct = events(1,[]);
% % %             [events(2:end).prev_event] = disperse(events(1:end-1));
% % %             events(1).prev_event = event_struct;
% % %             events_all_per_session{ii_exp_type,ii_epoch_type}{ii_exp} = events;
% % %         end
% % %     end
% % % end
% % % 
% % % %%
% % % features_names = {'duration';'compression';'distance';'distance_norm';};
% % % xlable_strs = {
% % %     'Replay duration (s)';
% % %     {'Compression ratio';'(replay speed / flight speed)'};
% % %     'Replay distance (m)';
% % %     {'Replay distance','(norm. to environment size)'};
% % %     };
% % % 
% % % %% plot boxplots - pool sleep/rest
% % % boxplot_panels_ylimits = [0 0.8; 0 25; 3 23; .02 .18];
% % % pvals = [];
% % % for ii_fn = 1:length(features_names)
% % %     axes(panels{4}(ii_fn));
% % %     cla
% % %     hold on
% % %     fn = features_names{ii_fn};
% % %     X={};
% % %     G={};
% % %     for ii_exp_type = 1:length(exp_types)
% % %         events1 = [events_all_per_session{ii_exp_type,1}{:}];
% % %         events2 = [events_all_per_session{ii_exp_type,2}{:}];
% % %         seqs1 = [events1.seq_model];
% % %         seqs2 = [events2.seq_model];
% % %         seqs = [seqs1 seqs2];
% % %         x = [seqs.(fn)];
% % %         g = ii_exp_type;
% % %         X{ii_exp_type} = x;
% % %         G{ii_exp_type} = ones(size(seqs)).*g;
% % %         m1 = prctile(x,50);
% % %         m2 = prctile(x,[25 75]);
% % %         m3 = prctile(x,[5 95]);
% % %         w = 0.2;
% % %         lw = 1.3;
% % %         clr = 'k';
% % %         plot([g-w g+w],[m1 m1],'Color',clr,'LineWidth',lw);
% % %         plot([g g],m3,'Color',clr,'LineWidth',lw);
% % %         rectangle('Position',[g-w m2(1) 2*w diff(m2)],'EdgeColor',clr,'FaceColor','none','LineWidth',lw);
% % %     end
% % %     pval = ranksum(X{1},X{2});
% % %     str = genSignifStrAstricks(pval);
% % %     xx = [g g-1];
% % %     yy = boxplot_panels_ylimits(ii_fn,[2 2]);
% % % %         hax= gca;
% % % %         yy = hax.YLim([2 2]);
% % %     plot(xx,yy,'k-');
% % %     font_size = 10;
% % %     if strcmp(str,'n.s.')
% % %         font_size = 8;
% % %         yy = yy.*1.02;
% % %     end
% % %     text(mean(xx),mean(yy),str,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',font_size);
% % %     pvals(ii_fn,ii_epoch_type) = pval;
% % %     xlim([0.2 2.8])
% % %     ylim(boxplot_panels_ylimits(ii_fn,:))
% % %     xlabel('');
% % %     ylabel_x_pos = [-0.65 -0.61 -0.55 -0.7];
% % %     ylabel(xlable_strs{ii_fn},'units','normalized','position',[ylabel_x_pos(ii_fn) 0.5]);
% % %     xticks(1:length(exp_types));
% % %     xticklabels(exp_types);
% % %     xtickangle(50);    
% % % end
% % % 
% % % end % if 0

%% add panel letters
font_size = 11;
axes(panels{1}(1))
text(-0.3,1.2, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(2))
axes(panels{2}(1))
text(-0.2,1.3, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{3}(1))
text(-0.38,1.3, 'c', 'Units','normalized','FontWeight','bold','FontSize',font_size);
% axes(panels{4}(1))
% text(-1.1,1.15, 'd', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{5}(1,1))
text(-0.35,1.15, 'e', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{5}(2,1))
text(-0.35,1.15, 'f', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = sprintf('%s',fig_name_str);
% fig_name = sprintf('%s_panel_B_opt_%d',fig_name,behavior_ex_opt);
fig_name = sprintf('%s_max_tdiff_%ds',fig_name,timediff_max);
[~,data_str,~] = fileparts(data_filename);
fig_name = sprintf('%s__%s',fig_name,data_str);
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%

%%
function str = genSignifStrAstricks(pval)
if pval < 0.001
    str = '***';
elseif pval < 0.01
    str = '**';
elseif pval < 0.05
    str = '*';
else
    str = 'n.s.';
end
end
