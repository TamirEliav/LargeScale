%% Replay - Fig 3 - replay population 
%%
clear 
clc
close all

%% plotting options
% replay_coverage_examples_IX = [19 21 22 23];
% replay_coverage_examples_IX = [24 25 27 29];
% replay_coverage_examples_IX = [41 53 54 55];
% replay_coverage_examples_IX = [61 65 68 69];
% replay_coverage_examples_IX = [70 71 73 74];
% replay_coverage_examples_IX = [75 77 78 81];
% replay_coverage_examples_IX = [21 19 70 78];
% replay_coverage_examples_IX = [21 19 70 29];
% replay_coverage_examples_IX = [21 19 70 29];
replay_coverage_examples_IX = [21 19 29 70];

% 5	b0184_d191202 21
% 1	b0184_d191130 19
% 31	b2382_d190801 70
% 30	b0184_d191212 29

single_session_seq_ex_IX = [21];

%% graphics params
epoch_type_clrs = {[.6 .1 .8],[.1 .8 .1]};
PastFutureClrs = [0.1.*[1 1 1]; 1 1 1];
directions_clrs = {[0    0.3843    0.7451];[ 0.5216    0.2471         0]};
panel_H_bin_size_options = [0.1 0.15 0.2 0.125]; % we chose opt 3 (0.2)
panel_H_bin_size_opt = 3;
panel_H_bin_size = panel_H_bin_size_options(panel_H_bin_size_opt)

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
close all
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
panels{1}(1) = axes('position', [3 23.3 3 2.5]);
panels{1}(2) = axes('position', [7 23.3 3 2.5]);
panels{1}(3) = axes('position', [11 23.3 3 2.5]);
panels{1}(4) = axes('position', [15 23.3 3 2.5]);

panels{2}(1) = axes('position', [3 11 3.5 10]);
panels{2}(2) = axes('position', [7 11 3.5 10]);

x = 12.3;
panels{3}(1) = axes('position', [x 19.5 7 1.5]);
panels{3}(2) = axes('position', [x 17.5 7 1.5]);
panels{3}(3) = axes('position', [x 15.5 7 1.5]);
panels{3}(4) = axes('position', [x 13.5 7 1.5]);
panels{4}(1) = axes('position', [x 11 7 1.5]);

panels{5}(1) = axes('position', [3     7.5 2 2]);
panels{5}(2) = axes('position', [4.1     9 .35 .35]);
panels{6}(1) = axes('position', [7.3    7.5 2 2]);
panels{7}(1,1) = axes('position', [10.9 7.5 2 2]);
panels{7}(2,1) = axes('position', [13.5 7.5 2 2]);
panels{8}(1) = axes('position', [17.5   7.5 2 2]);
% panels{8}(2) = axes('position', [17.5  4 2 2]);
% panels{8}(3) = axes('position', [17.5  1.5 3 2]);

total_offset = [0 0.5]+[0 -3];
for ii = 1:length(panels)
    subpanels = panels{ii};
    subpanels = subpanels(:);
    for jj = 1:length(subpanels)
        subpanels(jj).Position([1 2]) = subpanels(jj).Position([1 2]) + total_offset;
    end
end


%% load data
epoch_types = {'sleep','rest'};
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
events_all_per_session = {};
for ii_epoch_type = 1:length(epoch_types)
    for ii_exp = 1:length(exp_list)
        exp_ID = exp_list{ii_exp};
        exp = exp_load_data(exp_ID,'details','rest');
        epoch_type = epoch_types{ii_epoch_type};
        params_opt = 11;
        event_type = 'posterior';
        [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
        [~, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
        events(~TF) = [];
        if isempty(events)
            continue;
        end
        if isfield(events,'rest_ball_num')
            [events.rest_ball_loc] = disperse(exp.rest.balls_loc([events.rest_ball_num]));
        else
            [events.rest_ball_num] = disperse(nan(size(events)));
            [events.rest_ball_loc] = disperse(nan(size(events)));
        end
        
        [events.recordingArena] = deal(exp.details.recordingArena);
        event_struct = events(1,[]);
        [events(2:end).prev_event] = disperse(events(1:end-1));
        events(1).prev_event = event_struct;
        events_all_per_session{ii_epoch_type}{ii_exp} = events;
    end
end

%%
events_all = {[events_all_per_session{1}{:}], [events_all_per_session{2}{:}]};
seqs_all = {[events_all{1}.seq_model], [events_all{2}.seq_model]};
gArena = {categorical({events_all{1}.recordingArena}), categorical({events_all{2}.recordingArena}) };
seqs_pooled = [seqs_all{:}];

%%
features_names = {'duration';'compression';'distance';'distance_norm';};
xlable_strs = {
    'Replay duration (s)';
    {'Compression ratio';'(replay speed / flight speed)'};
    'Replay distance (m)';
    {'Replay distance','(norm. to environment size)'};
    };

%% panels A-D
clc
line_styles = {'-',':'};
panels_xlim = [-0.0950    1.9950; -1.9500   40.9500; -0.4500   56.5500; -0.0210    0.4410];
panels_xticks = {[0 0.5 1 1.5],[0 20 40],[0:10:50],[0 0.2 0.4]};
disp('panels a-d:')
for ii_fn = 1:length(features_names)
    fn = features_names{ii_fn};
    axes(panels{1}(ii_fn));
    cla
    hold on
    lw = 1.1;
    for ii_epoch_type = 1:length(epoch_types)
        X = [seqs_all{ii_epoch_type}.(fn)];
        g = gArena{ii_epoch_type};
%         histogram(X,'DisplayStyle','stairs','Normalization','pdf','LineStyle',line_styles{ii_epoch_type}, 'EdgeColor','k','LineWidth',lw);
        h1=histogram(X(g=='200m'),'DisplayStyle','stairs','Normalization','pdf','LineStyle',line_styles{1}, 'EdgeColor',epoch_type_clrs{ii_epoch_type},'LineWidth',lw);
        h2=histogram(X(g=='120m'),'DisplayStyle','stairs','Normalization','pdf','LineStyle',line_styles{2}, 'EdgeColor',epoch_type_clrs{ii_epoch_type},'LineWidth',lw);
        fprintf('%s: median=%.2g  mean=%.2g (%s,200m)\n', fn, median(X(g=='200m')), mean(X(g=='200m')), epoch_types{ii_epoch_type} )
        fprintf('%s: median=%.2g  mean=%.2g (%s,120m)\n', fn, median(X(g=='120m')), mean(X(g=='120m')), epoch_types{ii_epoch_type} )
    end
%     xlabel_pos = [];
    xlabel(xlable_strs{ii_fn}, 'Units','normalized', 'Position',[0.5 -0.14]);
%     xlabel(xlable_strs{ii_fn});
    if ii_fn == 1
        ylabel('Probability density');
    end
    hax=gca;
    hax.XLim = panels_xlim(ii_fn,:);
    hax.XTick = panels_xticks{ii_fn};
    hax.XRuler.TickLength(1) = 0.03;
    hax.XRuler.TickLabelGapOffset = -1;
    hax.YRuler.TickLabelGapOffset = 1;
    fprintf('%s: median = %.2g, mean = %.2g (pooled)\n',fn,median([seqs_pooled.(fn)]),mean([seqs_pooled.(fn)]))
end

%% legend (hists)
if exist('panels_hist_legend','var')
    delete(panels_hist_legend);
end
panels_hist_legend = axes('position', [4 22 0.5 0.5]);
cla
hold on
t = linspace(0,1,100);
x = pulstran(t,linspace(0,1,3),'rectpuls',1/6);
x(x>0) = nan;x(~isnan(x)) = 0;
clear h
% plot(  t  ,   x  ,       'color', epoch_type_clrs{2}, 'LineWidth',lw,'Clipping','off');
% plot(  t  ,   x  +1.5,   'color', epoch_type_clrs{1}, 'LineWidth',lw,'Clipping','off');
plot([0 1], [0 0],     ':', 'color', epoch_type_clrs{2}, 'LineWidth',lw,'Clipping','off');
plot([0 1], [0 0]+1.5, ':', 'color', epoch_type_clrs{1}, 'LineWidth',lw,'Clipping','off');
plot([0 1], [.5 .5],     'color', epoch_type_clrs{2}, 'LineWidth',lw,'Clipping','off');
plot([0 1], [.5 .5]+1.5, 'color', epoch_type_clrs{1}, 'LineWidth',lw,'Clipping','off');
text(1.3, 2.0, 'Sleep (200m)','FontSize',7,'HorizontalAlignment','left');
text(1.3, 1.5, 'Sleep (130m)','FontSize',7,'HorizontalAlignment','left');
text(1.3, 0.5, 'Awake (200m)','FontSize',7,'HorizontalAlignment','left');
text(1.3, 0.0, 'Awake (130m)','FontSize',7,'HorizontalAlignment','left');
xlim([0 1])
ylim([0 1])
axis off


%% panel F - coverage per session
load('E:\Tamir\work\PROJECTS\LargeScale\paper_replay\data_prepared_for_figures\replay_coverage.mat');
for ii_ex = 1:length(panels{3})
    axes(panels{3}(ii_ex))
    cla
    hold on
    ii_exp = replay_coverage_examples_IX(ii_ex);
    c = squeeze(coverage_all(ii_exp,:,:,:));
    x = linspace(0,1,size(c,3));
    lw = 1.5;
    plot(x, squeeze(c(1,1,:)),'-','LineWidth',lw,'Color',directions_clrs{1});
    plot(x, squeeze(c(2,1,:)),'--','LineWidth',lw,'Color',directions_clrs{1});
    plot(x, squeeze(c(1,2,:)),'-','LineWidth',lw,'Color',directions_clrs{2});
    plot(x, squeeze(c(2,2,:)),'--','LineWidth',lw,'Color',directions_clrs{2});
    hax=gca;
    hax.XTick = [0 1];
    hax.XLim = [0 1];
    hax.XRuler.TickLabelGapOffset = -1.5;
    hax.YRuler.TickLabelGapOffset = 1;
    if ii_ex == 1
        title('Four example sessions: replay spatial coverage','Units','normalized','Position',[0.5 1.14],'FontWeight','normal','HorizontalAlignment','center');
%         hl=legend({'Sleep dir 1','Rest dir 1','Sleep dir 2','Rest dir 2'},'NumColumns',2);
%         hl.Units = 'normalized';
%         hl.Position([1 2]) = [0.22 0.715];
%         hl.Position([3 4]) = [0.1 0.025];
%         hl.Box = 'off';
    end
    if ii_ex == length(panels{3})
        xlabel('Position (norm.)', 'Units','normalized', 'Position',[0.5 -0.06]);
        ylabel({'Replay coverage';'(counts)'}, 'Units','normalized', 'Position',[-0.1 2.4]);
    end
end

%% panel F legend
if exist('panels_coverage_legend','var')
    delete(panels_coverage_legend);
end
panels_coverage_legend(1) = axes('position', [panels{3}(1).Position([1 2])+[.5 1.25] 0.5 0.5]);
panels_coverage_legend(2) = axes('position', [panels{3}(1).Position([1 2])+[2.6 1.25] 0.5 0.5]);
arrow_str = {'\rightarrow','\leftarrow'};
for ii_dir = 1:2
    axes(panels_coverage_legend(ii_dir))
    cla
    hold on
    t = linspace(0,1,100);
    x = pulstran(t,linspace(0,1,3),'rectpuls',1/6);
    x(x>0) = nan;x(~isnan(x)) = 0;
    clear h
    plot(  t  ,   x,     'color', directions_clrs{ii_dir}, 'LineWidth',lw,'Clipping','off');
    plot([0 1], [.5 .5], 'color', directions_clrs{ii_dir}, 'LineWidth',lw,'Clipping','off');
    text(1.3, 0.5, "Sleep dir "+ii_dir,'FontSize',7,'HorizontalAlignment','left');
    text(1.3, 0, "Awake dir "+ii_dir,'FontSize',7,'HorizontalAlignment','left');
    text(0.2, -0.3, arrow_str{ii_dir},'Color',directions_clrs{ii_dir},'FontWeight','bold', 'FontSize',10);
    xlim([0 1])
    ylim([0 1])
    axis off
end

%% panel G - coverage total
axes(panels{4}(1))
cla
hold on
c = squeeze(sum(coverage_all,1));
x = linspace(0,1,size(c,3));
lw = 2;
plot(x, squeeze(c(1,1,:)),'-','LineWidth',lw,'Color',directions_clrs{1});
plot(x, squeeze(c(2,1,:)),'--','LineWidth',lw,'Color',directions_clrs{1});
plot(x, squeeze(c(1,2,:)),'-','LineWidth',lw,'Color',directions_clrs{2});
plot(x, squeeze(c(2,2,:)),'--','LineWidth',lw,'Color',directions_clrs{2});
hax=gca;
hax.XTick = [0 1];
hax.XLim = [0 1];
hax.XRuler.TickLabelGapOffset = -1.5;
hax.YRuler.TickLabelGapOffset = 1;
xlabel('Position (norm.)', 'Units','normalized', 'Position',[0.5 -0.06]);
ylabel({'Replay coverage';'(counts)'}, 'Units','normalized', 'Position',[-0.1 .5]);
title("Population: {\itn} = "+size(coverage_all,1)+" sessions",'Units','normalized','Position',[0.5 0.9],'FontWeight','normal');

%% panel E  - all seqs in a session (sleep)
axes(panels{2}(1))
hax=gca;
cla
hold on
exp_ID = T(single_session_seq_ex_IX,:).exp_ID{1};
epoch_type ='sleep';
params_opt = 11;
events = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
[seqs, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
events(~TF) = [];
epoch_sep = find(diff([events.epoch_num])~=0)+0.5;
seqs_edges = [seqs.start_pos_norm; seqs.end_pos_norm];
seqs_IX = 1:length(seqs);
plot([0 1],repmat(epoch_sep,2,1),':','LineWidth',1.5,'Color',0.5*[1 1 1]);
h=plot(seqs_edges,[seqs_IX;seqs_IX],'-','LineWidth',.55);
dir_1_IX = [seqs.state_direction]==1;
dir_2_IX = [seqs.state_direction]==-1;
[h(dir_1_IX).Color] = disperse(repelem(directions_clrs(1),length(dir_1_IX)));
[h(dir_2_IX).Color] = disperse(repelem(directions_clrs(2),length(dir_2_IX)));
scatter(seqs_edges(1,:),seqs_IX, 3, [seqs.state_direction]==-1, "filled");
hax.Colormap = cell2mat(directions_clrs);
hax=gca;
hax.XTick = [0 1];
hax.YTick = [1 10*ceil(length(seqs)/10)];
hax.XLim = [0 1];
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = 1;
xlabel('Position (norm.)', 'Units','normalized', 'Position',[0.5 -0.017]);
ylabel('Replay event no.', 'Units','normalized', 'Position',[-0.06 .5]);
title('Sleep replays', 'Units','normalized', 'Position',[0.47 0.99],'FontWeight','normal');
text(1.043, 1.05, 'Example session: all individual replays','FontSize',7,'HorizontalAlignment','center','Units','normalized')

%% Panel E - all seqs in a session (rest)
axes(panels{2}(2))
hax=gca;
cla
hold on
exp_ID = T(single_session_seq_ex_IX,:).exp_ID{1};
epoch_type ='rest';
params_opt = 11;
events = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
[seqs, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
events(~TF)=[];
epoch_sep = find(diff([events.epoch_num])~=0)+0.5;
seqs_edges = [seqs.start_pos_norm; seqs.end_pos_norm];
seqs_IX = 1:length(seqs);
% plot([0 1],repmat(epoch_sep,2,1),'-','LineWidth',.2,'Color',0.9*[1 1 1]);
h=plot(seqs_edges,[seqs_IX;seqs_IX],'-','LineWidth',.55);
dir_1_IX = [seqs.state_direction]==1;
dir_2_IX = [seqs.state_direction]==-1;
[h(dir_1_IX).Color] = disperse(repelem(directions_clrs(1),length(dir_1_IX)));
[h(dir_2_IX).Color] = disperse(repelem(directions_clrs(2),length(dir_2_IX)));
scatter(seqs_edges(1,:),seqs_IX, 3, [seqs.state_direction]==-1, "filled");
hax.Colormap = cell2mat(directions_clrs);
hax.XTick = [0 1];
hax.YTick = [1 10*ceil(length(seqs)/10)];
hax.XLim = [0 1];
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = 1;
xlabel('Position (norm.)', 'Units','normalized', 'Position',[0.5 -0.017]);
% ylabel('Replay event no.', 'Units','normalized', 'Position',[-0.09 .5]);
title('Awake replays', 'Units','normalized', 'Position',[0.47 0.99],'FontWeight','normal');

%% Coverage correlations (pdf)
axes(panels{5}(1))
cla reset
hold on
sleep_vs_rest_coverage_corr = ccc_all(:,:,1,2);
sleep_vs_rest_coverage_corr = sleep_vs_rest_coverage_corr(:);
h=histogram(sleep_vs_rest_coverage_corr,'FaceColor',0.15*[1 1 1],'Normalization','pdf');
h.BinWidth = panel_H_bin_size;
h.BinLimits = [-1 1];
nbins = 50;
% h=histogram(cat(1,ccc_shuffled{:}),'FaceColor','k','Normalization','pdf');
% h.DisplayStyle = 'stairs';
% h.BinLimits = [-1 1];
% h.EdgeColor = 'm';
% h.LineWidth = 1.5;
% h.NumBins = nbins;
xi = linspace(-1,1,nbins);
f = ksdensity(cat(1,ccc_shuffled{:}),xi,'Function','pdf');
plot(xi,f,'k','LineWidth',1)
[~,P_KS, KS_stats] = kstest2(sleep_vs_rest_coverage_corr, cat(1,ccc_shuffled{:}));
text(0.5,1.15,['{\itP}_{KS} = ' sprintf('%.2g',P_KS)],'units','normalized','FontSize',7,'HorizontalAlignment','center');
hax=gca;
hax.XLim = [-1 1];
hax.XTick = -1:0.5:1;
hax.YTick = [0 hax.YTick(end)];
hax.TickDir = 'out';
hax.XRuler.TickLength = [0.03 0];
hax.YRuler.TickLength = [0.02 0];
hax.XRuler.TickLabelGapOffset = -0;
hax.YRuler.TickLabelGapOffset = 1;
hax.XTickLabelRotation = 0;
xlabel({'Corr. of replay positions (r)';'sleep vs. awake'}, 'Units','normalized', 'Position',[0.5 -0.25]);
ylabel({'Probability','density'}, 'Units','normalized', 'Position',[-0.08 .5]);

%% add legend
axes(panels{5}(2))
cla reset
hold on
axis off
axis equal
patch([0 0 1 1],[0 1 1 0], 0.15*[1 1 1],'EdgeColor','k');
plot([0 1],[1 1].*2, 'k-','LineWidth',1.5)
text(1.5,0.5, 'Data','FontSize',7)
text(1.5,2, 'Shuffle','FontSize',7)
axis ij


%% replay rate - PRE vs POST sleep
axes(panels{6}(1))
cla reset
hold on
plot(replay_rate(:,3),replay_rate(:,4),'ok','MarkerSize',2)
% plot(replay_rate(:,3),replay_rate(:,4),'.k')
plot([1e-3 1e0],[1e-3 1e0],'Color',[1 1 1].*0.5);
pval = signrank(replay_rate(:,3),replay_rate(:,4));
disp('panel i stats:')
fprintf('replay-rate sleep-after vs before:\n');
fprintf('median ratio = %.2g\n', nanmedian(replay_rate(:,4)./replay_rate(:,3)) );
fprintf('pval (signrank test)= %.2g\n', pval );
% replay_rate_diff_sleep = diff(replay_rate(:,[3 4]),1,2);
% h=histogram(replay_rate_diff_sleep,'FaceColor','k','Normalization','pdf');
% pval = signrank(replay_rate_diff_sleep);
text(0.5,1.15,['{\itP}_{signrank} = ' sprintf('%.1g',pval)],'units','normalized','FontSize',7,'HorizontalAlignment','center');
axis square
axis equal
hax=gca;
hax.XScale = 'log';
hax.YScale = 'log';
hax.XLim = [5e-4 5e-1];
hax.YLim = [5e-4 5e-1];
hax.XTick = [1e-3 1e-2 1e-1 1e0];
hax.YTick = [1e-3 1e-2 1e-1 1e0];
hax.TickDir = 'out';
hax.XRuler.TickLength = [0.04 0];
hax.YRuler.TickLength = [0.04 0];
hax.XRuler.TickLabelGapOffset = -0;
hax.YRuler.TickLabelGapOffset = 1;
hax.XTickLabelRotation = 0;
xlabel({'Replay rate (events/s)';'sleep before'}, 'Units','normalized', 'Position',[0.5 -0.25]);
ylabel({'Replay rate (events/s)';'sleep after'}, 'Units','normalized', 'Position',[-0.38 0.5]);

%% replay i vs i+1 (no corr) - OLD!!
if 0
bins_size = 20; % in meters
for ii_epoch_type = 1:length(epoch_types)
    events2 = events_all{ii_epoch_type};
    invalid = cellfun(@isempty,{events2.prev_event});
    events2(invalid)=[];
    events1 = [events2.prev_event];
    seqs2 = [events2.seq_model];
    seqs1 = [events1.seq_model];
    xxx = [seqs1.middle_pos];
    yyy = [seqs2.middle_pos];

    % scatter 
    axes(panels{7}(ii_epoch_type,1))
    cla reset
    hold on
    axis equal
    plot(xxx,yyy,'k.','MarkerSize',2,'Color',epoch_type_clrs{ii_epoch_type});
    h=refline(1,0);
    h.Color = 0.5*[1 1 1];
%     xlabel({'Position of';'replay {\iti} (m)'},'Units','normalized','Position',[0.5 -0.3]);
    if ii_epoch_type == 1
        xlabel({'Position of replay {\iti} (m)'},'Units','normalized','Position',[1.15 -0.3]);
%         xlabel({'Position of';'replay {\iti} (m)'},'Units','normalized','Position',[1 -0.3]); 
        ylabel({'Position of';'replay {\iti+1} (m)'},'Units','normalized','Position',[-0.35 0.5])
    end
    hax=gca;
    hax.XRuler.TickLength(1) = 0.05;
    hax.YRuler.TickLength(1) = 0.03;
    hax.XRuler.TickLabelGapOffset = -0;
    hax.YRuler.TickLabelGapOffset = 1;

    % density
    axes(panels{7}(ii_epoch_type,2))
    cla reset
    hold on
    axis equal
    edges = 0:bins_size:200;
    xedges = 0:bins_size:200;
    yedges = 0:bins_size:201;
    xedges = edges;
    yedges = edges;
    N1 = histcounts(xxx,xedges,'Normalization','pdf');
    N2 = histcounts(yyy,yedges,'Normalization','pdf');
    N12 = histcounts2(xxx,yyy,xedges,yedges,'Normalization','pdf');
    N12 = N12 ./ (N1'.*N2);
    imagesc(xedges,yedges,log(N12));
%     imagesc(xedges,yedges,N12);
%     colormap bone
    cmap = jet;
%     cmap = 1-cmap;
    colormap(cmap)
    hax=gca;
    hax.CLim(1) = 0;

    axes(panels{7}(ii_epoch_type,1))
    switch epoch_types{ii_epoch_type}
        case 'sleep'
            title('Sleep','FontWeight','normal');
        case 'rest'
            title('Awake','FontWeight','normal');
    end
end
end


%% Replay i vs i+1 (no corr)
% sigma = 1.5;
% bins_size = 20;
% edges = 0:bins_size:200;
% xedges = 0:bins_size:200;
% yedges = 0:bins_size:201;
edges = linspace(0,1,11);
xedges = edges;
yedges = edges;
nbinsx = length(xedges)-1;
nbinsy = length(yedges)-1;
nEpochTypes = length(events_all_per_session);
nSessions = length(events_all_per_session{1});
M = nan(nEpochTypes,nSessions,nbinsx,nbinsy);
events_count = zeros(nEpochTypes,nSessions);
for ii_epoch_type = 1:nEpochTypes
    for ii_session = 1:nSessions
        events2 = events_all_per_session{ii_epoch_type}{ii_session};
        if length(events2) == 0
            continue
        end
        invalid = cellfun(@isempty,{events2.prev_event});
        events2(invalid)=[];
        if length(events2) == 0
            continue
        end
        events1 = [events2.prev_event];
        seqs2 = [events2.seq_model];
        seqs1 = [events1.seq_model];
        xxx = [seqs1.middle_pos_norm];
        yyy = [seqs2.middle_pos_norm];
        N1 = histcounts(xxx,xedges,'Normalization','pdf');
        N2 = histcounts(yyy,yedges,'Normalization','pdf');
        N12 = histcounts2(xxx,yyy,xedges,yedges,'Normalization','pdf');
        N12 = N12 ./ (N1'.*N2);
        M(ii_epoch_type,ii_session,:,:) = N12;
        events_count(ii_epoch_type,ii_session) = length(xxx);
    end
end
W = events_count./sum(events_count')';
M_weighted = M.*W;
M_weighted_avg = squeeze(nanmean(M_weighted,[2]));

for ii_epoch_type = 1:nEpochTypes

    M2 = squeeze(M_weighted_avg(ii_epoch_type,:,:));
%     M3 = imgaussfilt(M2,sigma,'FilterDomain','spatial');

    axes(panels{7}(ii_epoch_type,1))
    cla reset
    hold on
    axis equal
%     imagesc(xedges,yedges,log(M2));
    imagesc(xedges,yedges,M2);
    cmap = jet;
    colormap(cmap)
    if ii_epoch_type == 1
        xlabel({'Position of replay {\iti} (norm.)'},'Units','normalized','Position',[1.2 -0.25]);
%         xlabel({'Position of';'replay {\iti} (m)'},'Units','normalized','Position',[1 -0.3]); 
        ylabel({'Position of';'replay {\iti+1} (norm.)'},'Units','normalized','Position',[-0.25 0.5])
    end
    hax=gca;
    hax.CLim(1) = 0;
    hax.XRuler.TickLength(1) = 0.05;
    hax.YRuler.TickLength(1) = 0.03;
    hax.XRuler.TickLabelGapOffset = -1;
    hax.YRuler.TickLabelGapOffset = 1;
    switch epoch_types{ii_epoch_type}
        case 'sleep'
            title('Sleep','FontWeight','normal')
        case 'rest'
            title('Awake','FontWeight','normal')
    end
    if ii_epoch_type==2
        hcb = colorbar('east');
        hcb.Units = 'centimeters';
        hcb.Position(1) = hax.Position(1) + hax.Position(3)*1.07;
        hcb.Ticks = [];
        hcb.Label.String = 'Density (norm.)';
        hcb.Label.Rotation = -90;
        hcb.Label.Units = 'normalized';
        hcb.Label.Position = [1.3 0.5];
        text(1.14,1+0.03,'Max','Units','normalized','FontSize',7,'HorizontalAlignment','center');
        text(1.14,0-0.03,'0','Units','normalized','FontSize',7,'HorizontalAlignment','center');
    end

%     axes(panels{6}(ii_epoch_type,2))
%     cla reset
%     hold on
%     axis equal
% %     imagesc(xedges,yedges,log(M3));
%     imagesc(xedges,yedges,M3);
%     cmap = jet;
%     colormap(cmap)
end

%%
% % % % n_shuffles = 1000;
% % % % rng(0)
% % % % % TODO: shuffle must be done per session (because over-representation is
% % % % % different per sessions)
% % % % shuffles = arrayfun(@randperm,repelem(length(xxx),n_shuffles),'UniformOutput',false);
% % % % shuffles = cat(1,shuffles{:});
% % % % corr_shuffled = corr(shuffles',yyy','tail','right');

%% Takeoff/Landing X Past/Future
axes(panels{8}(1))
cla reset
hold on
% take only rest (awake)
seqs = seqs_all{2};
events = events_all{2};
TakeLand_thr = 0.05;
gTakeLand = classify_replay_landing_takeoff_other(seqs, TakeLand_thr);
gPastFuture = categorical( ([events.rest_ball_num] == 1 & [seqs.state_direction] == 1) | ...
                           ([events.rest_ball_num] == 2 & [seqs.state_direction] == -1), ...
                           [false true],["Past","Future"])';
gForRev = categorical([seqs.forward],[true false],["Forward","Reverse"])';
g = gTakeLand .* gPastFuture;
[N,G] = histcounts(g);
N = reshape(N,2,3)';
G = reshape(G,2,3)';
% do some stats
TakeoffsFuturePast_binom_pval = myBinomTest(sum(gTakeLand=='Takeoff' & gPastFuture=='Future'),sum(gTakeLand=='Takeoff'),0.5,'two');
LandingFuturePast_binom_pval = myBinomTest(sum(gTakeLand=='Landing' & gPastFuture=='Past'),sum(gTakeLand=='Landing'),0.5,'two');
disp('panel 2k stats:')
fprintf('Takeoff: past vs future, pval=%.2g (binom test)\n',TakeoffsFuturePast_binom_pval);
fprintf('Landing: past vs future, pval=%.2g (binom test)\n',LandingFuturePast_binom_pval);
takeoff_IX = find(all(contains(G,'Takeoff'),2));
midair_IX = find(all(contains(G,'Mid-air'),2));
landing_IX = find(all(contains(G,'Landing'),2));
new_order = [takeoff_IX midair_IX landing_IX];
N = N(new_order,:);
G = G(new_order,:);
N([1 3],:) = N([1 3],:)./TakeLand_thr;
N([2],:) = N([2],:)./(1-2*TakeLand_thr);
N = N ./ sum(N,"all");
hb=bar(N);
[hb.FaceColor] = disperse(PastFutureClrs');
% hb(1).FaceColor = 1.0.*[1 1 1];
% hb(1).FaceColor = 0.1.*[1 1 1];
hax=gca;
hax.XTick = 1:3;
hax.XTickLabel = {'Takeoff zone','Mid-air zone','Landing zone'};
hax.XRuler.TickLabelGapOffset = -1;
ylabel('Fraction')
% create legend
w = 0.08*range(hax.XLim);
h = 0.10*range(hax.YLim);
x = hax.XLim(1)+0.1*range(hax.XLim);
y1 = hax.YLim(1)+0.9*range(hax.YLim);
y2 = hax.YLim(1)+0.75*range(hax.YLim);
rectangle(Position=[x y1 w h],FaceColor=hb(1).FaceColor);
rectangle(Position=[x y2 w h],FaceColor=hb(2).FaceColor);
text(x+1.3*w,y1+h/2,'Past','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',7);
text(x+1.3*w,y2+h/2,'Future','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',7);

% hl=legend(["past","Future"],'box','off','location','none')
% hl.Units = hax.Units;
% hl.Position = [hax.Position([1 2])+[0.5 0.75].*hax.Position([3 4]) .5 .5];
% hl


%% add panel letters
font_size = 11;
axes(panels{1}(1))
text(-0.2,1.12, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(2))
text(-0.2,1.12, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(3))
text(-0.2,1.12, 'c', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(4))
text(-0.2,1.12, 'd', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{2}(1))
text(-0.25,1.035, 'e', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{3}(1))
text(-0.1,1.3, 'f', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{4}(1))
text(-0.1,1.18, 'g', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{5}(1))
text(-0.3,1.2, 'h', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{6}(1))
text(-0.85,1.2, 'i', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{7}(1))
text(-0.6,1.2, 'j', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{8}(1))
text(-0.45,1.2, 'k', 'Units','normalized','FontWeight','bold','FontSize',font_size);



%%
fig_name = sprintf('%s_panel_E_%d_%d_%d_%d_panel_G_%d_panel_H_bin_opt_%d', ...
    fig_name_str, ...
    replay_coverage_examples_IX, ...
    single_session_seq_ex_IX,...
    panel_H_bin_size_opt)
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%
