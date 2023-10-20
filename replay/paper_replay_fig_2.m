%% Replay - Fig 2 - replay population 
%%
clear 
clc

%% plotting options
% replay_coverage_examples_IX = [19 21 22 23];
% replay_coverage_examples_IX = [24 25 27 29];
% replay_coverage_examples_IX = [41 53 54 55];
% replay_coverage_examples_IX = [61 65 68 69];
% replay_coverage_examples_IX = [70 71 73 74];
% replay_coverage_examples_IX = [75 77 78 81];
replay_coverage_examples_IX = [21 19 70 78]

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
fig_name_str = 'Fig_2';
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
panels{1}(1) = axes('position', [3 23 3 3]);
panels{1}(2) = axes('position', [7 23 3 3]);
panels{1}(3) = axes('position', [11 23 3 3]);
panels{1}(4) = axes('position', [15 23 3 3]);

panels{2}(1) = axes('position', [3 19.5 7 1.5]);
panels{2}(2) = axes('position', [3 17.5 7 1.5]);
panels{2}(3) = axes('position', [3 15.5 7 1.5]);
panels{2}(4) = axes('position', [3 13.5 7 1.5]);

panels{3}(1) = axes('position', [3 11 7 1.5]);

panels{4}(1) = axes('position', [11.0 11 3.5 10]);
panels{4}(2) = axes('position', [15.0 11 3.5 10]);

panels{5}(1) = axes('position', [3   7.5 2.5 2]);
panels{6}(1,1) = axes('position', [7.5 7.5 2 2]);
panels{6}(2,1) = axes('position', [10.2 7.5 2 2]);
% panels{6}(1,2) = axes('position', [7.5 4.3 2 2]);
% panels{6}(2,2) = axes('position', [10.2 4.3 2 2]);
% panels{6}(1,3) = axes('position', [7.5 1.5 2 2]);
% panels{6}(2,3) = axes('position', [10.2 1.5 2 2]);
panels{7}(1) = axes('position', [13.5  7.5 3 2]);

%% load data
epoch_types = {'sleep','rest'};
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
events_all_per_session = {};
for ii_epoch_type = 1:length(epoch_types)
    for ii_exp = 1:length(exp_list)
        exp_ID = exp_list{ii_exp};
        exp = exp_load_data(exp_ID,'details');
        epoch_type = epoch_types{ii_epoch_type};
        params_opt = 11;
        event_type = 'posterior';
        [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
        [~, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
        events(~TF) = [];
        if isempty(events)
            continue;
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

%%
features_names = {'duration';'compression';'distance';'distance_norm';};
xlable_strs = {
    'Replay duration (s)';
    {'Compression ratio';'(replay speed / flight speed)'};
    'Replay distance (m)';
    {'Replay distance','(norm. to environment size)'};
    };

%% panels A-D
axes(panels{1}(1))
hold on
line_styles = {'-','--'};
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
        histogram(X(g=='200m'),'DisplayStyle','stairs','Normalization','pdf','LineStyle',line_styles{1}, 'EdgeColor',epoch_type_clrs{ii_epoch_type},'LineWidth',lw);
        histogram(X(g=='120m'),'DisplayStyle','stairs','Normalization','pdf','LineStyle',line_styles{2}, 'EdgeColor',epoch_type_clrs{ii_epoch_type},'LineWidth',lw);
    end
%     xlabel_pos = [];
%     xlabel(xlable_strs{ii_fn}, 'Units','normalized', 'Position',[0.5 -0.11]));
    xlabel(xlable_strs{ii_fn});
    if ii_fn == 1
        ylabel('Probability density');
    end
    hax=gca;
%     hax.TickLength(1) = [0.025];
    hax.XRuler.TickLength(1) = 0.03;
    hax.XRuler.TickLabelGapOffset = -1.2;
    hax.YRuler.TickLabelGapOffset = 1;
end

%% legend (hists)
if exist('panels_hist_legend','var')
    delete(panels_hist_legend);
end
panels_hist_legend = axes('position', [4 25 0.5 0.5]);
cla
hold on
t = linspace(0,1,100);
x = pulstran(t,linspace(0,1,3),'rectpuls',1/6);
x(x>0) = nan;x(~isnan(x)) = 0;
clear h
plot(  t  ,   x  ,       'color', epoch_type_clrs{2}, 'LineWidth',lw,'Clipping','off');
plot(  t  ,   x  +1.5,   'color', epoch_type_clrs{1}, 'LineWidth',lw,'Clipping','off');
plot([0 1], [.5 .5],     'color', epoch_type_clrs{2}, 'LineWidth',lw,'Clipping','off');
plot([0 1], [.5 .5]+1.5, 'color', epoch_type_clrs{1}, 'LineWidth',lw,'Clipping','off');
text(1.3, 2.0, 'Sleep (200m)','FontSize',7,'HorizontalAlignment','left');
text(1.3, 1.5, 'Sleep (130m)','FontSize',7,'HorizontalAlignment','left');
text(1.3, 0.5, 'Rest (200m)','FontSize',7,'HorizontalAlignment','left');
text(1.3, 0.0, 'Rest (130m)','FontSize',7,'HorizontalAlignment','left');
xlim([0 1])
ylim([0 1])
axis off


%% panel E - coverage per session
load('E:\Tamir\work\PROJECTS\LargeScale\paper_replay\data_prepared_for_figures\replay_coverage.mat');
for ii_ex = 1:length(panels{2})
    axes(panels{2}(ii_ex))
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
        title('Example sessions','Units','normalized','Position',[0.2 1.1],'FontWeight','normal');
%         hl=legend({'Sleep dir 1','Rest dir 1','Sleep dir 2','Rest dir 2'},'NumColumns',2);
%         hl.Units = 'normalized';
%         hl.Position([1 2]) = [0.22 0.715];
%         hl.Position([3 4]) = [0.1 0.025];
%         hl.Box = 'off';
    end
    if ii_ex == length(panels{2})
        xlabel('Position (norm.)', 'Units','normalized', 'Position',[0.5 -0.11]);
        ylabel({'Replay coverage';'(counts)'}, 'Units','normalized', 'Position',[-0.1 2.4]);
    end
end

%% panel E legend
if exist('panels_coverage_legend','var')
    delete(panels_coverage_legend);
end
panels_coverage_legend(1) = axes('position', [panels{2}(1).Position([1 2])+[.5 1.25] 0.5 0.5]);
panels_coverage_legend(2) = axes('position', [panels{2}(1).Position([1 2])+[2.5 1.25] 0.5 0.5]);
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
    text(1.3, 0, "Rest dir "+ii_dir,'FontSize',7,'HorizontalAlignment','left');
    text(0.2, -0.3, arrow_str{ii_dir},'Color',directions_clrs{ii_dir},'FontWeight','bold', 'FontSize',10);
    xlim([0 1])
    ylim([0 1])
    axis off
end


%% panel F - coverage total
axes(panels{3}(1))
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
xlabel('Position (norm.)', 'Units','normalized', 'Position',[0.5 -0.11]);
ylabel({'Replay coverage';'(counts)'}, 'Units','normalized', 'Position',[-0.1 .5]);
title("{\itn} = "+size(coverage_all,1)+" sessions",'Units','normalized','Position',[0.2 1.05],'FontWeight','normal');

%% panel G  - all seqs in a session (sleep)
axes(panels{4}(1))
hax=gca;
cla
hold on
exp_ID = T(single_session_seq_ex_IX,:).exp_ID{1};
epoch_type ='sleep';
params_opt = 11;
events = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
[seqs, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
seqs_edges = [seqs.start_pos_norm; seqs.end_pos_norm];
seqs_IX = 1:length(seqs);
h=plot(seqs_edges,[seqs_IX;seqs_IX],'-','LineWidth',.1);
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
xlabel('Position (norm.)', 'Units','normalized', 'Position',[0.5 -0.02]);
ylabel('Replay event no.', 'Units','normalized', 'Position',[-0.1 .5]);
title('Example sleep session', 'Units','normalized', 'Position',[0.5 1],'FontWeight','normal');

%% Panel G - all seqs in a session (rest)
axes(panels{4}(2))
hax=gca;
cla
hold on
exp_ID = T(single_session_seq_ex_IX,:).exp_ID{1};
epoch_type ='rest';
params_opt = 11;
events = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
[seqs, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
events(~TF)=[];
seqs_edges = [seqs.start_pos_norm; seqs.end_pos_norm];
seqs_IX = 1:length(seqs);
h=plot(seqs_edges,[seqs_IX;seqs_IX],'-');
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
xlabel('Position (norm.)', 'Units','normalized', 'Position',[0.5 -0.02]);
% ylabel('Replay event no.', 'Units','normalized', 'Position',[-0.09 .5]);
title('Example rest session', 'Units','normalized', 'Position',[0.5 1],'FontWeight','normal');

%% Panel H - coverage correlations
axes(panels{5}(1))
cla reset
hold on
h=histogram(ccc_all(:),'FaceColor','k');
h.BinWidth = panel_H_bin_size;
h.BinLimits = [-1 1];;
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
xlabel({'Coverage correlation';'sleep vs. rest'}, 'Units','normalized', 'Position',[0.5 -0.25]);
ylabel('Counts', 'Units','normalized', 'Position',[-0.12 .5]);

%% panel I - replay i vs i+1 (no corr)
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
    axes(panels{6}(ii_epoch_type,1))
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
    axes(panels{6}(ii_epoch_type,2))
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

    axes(panels{6}(ii_epoch_type,1))
    title(epoch_types{ii_epoch_type})
end
end


%% panel I - replay i vs i+1 (no corr)
bins_size = 20;
edges = 0:bins_size:200;
xedges = 0:bins_size:200;
yedges = 0:bins_size:201;
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
        xxx = [seqs1.middle_pos];
        yyy = [seqs2.middle_pos];
        N1 = histcounts(xxx,xedges,'Normalization','pdf');
        N2 = histcounts(yyy,yedges,'Normalization','pdf');
        N12 = histcounts2(xxx,yyy,xedges,yedges,'Normalization','pdf');
        N12 = N12 ./ (N1'.*N2);
        M(ii_epoch_type,ii_session,:,:) = N12;
        events_count(ii_epoch_type,ii_session) = length(xxx);
    end
end
W = events_count./sum(events_count')';
M_weighted = M./W;
M_weighted_avg = squeeze(nanmean(M_weighted,[2]));

for ii_epoch_type = 1:nEpochTypes

    sigma = 1.5;
    M2 = squeeze(M_weighted_avg(ii_epoch_type,:,:));
    M3 = imgaussfilt(M2,sigma);

    axes(panels{6}(ii_epoch_type,1))
    cla reset
    hold on
    axis equal
%     imagesc(xedges,yedges,log(M2));
    imagesc(xedges,yedges,M2);
    cmap = jet;
    colormap(cmap)
    if ii_epoch_type == 1
        xlabel({'Position of replay {\iti} (m)'},'Units','normalized','Position',[1.2 -0.25]);
%         xlabel({'Position of';'replay {\iti} (m)'},'Units','normalized','Position',[1 -0.3]); 
        ylabel({'Position of';'replay {\iti+1} (m)'},'Units','normalized','Position',[-0.35 0.5])
    end
    hax=gca;
    hax.XRuler.TickLength(1) = 0.05;
    hax.YRuler.TickLength(1) = 0.03;
    hax.XRuler.TickLabelGapOffset = -1;
    hax.YRuler.TickLabelGapOffset = 1;
    title(epoch_types{ii_epoch_type})

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

%% Panel J - Takeoff/Landing X Past/Future
axes(panels{7}(1))
cla reset
hold on
% take only rest 
seqs = seqs_all{2};
events = events_all{2};
TakeLand_thr = 0.05;
gTakeLand = classify_replay_landing_takeoff_other(seqs, TakeLand_thr);
gPastFuture = categorical( ([events.rest_ball_num] == 1 & [seqs.state_direction] == 1) | ...
                           ([events.rest_ball_num] == 2 & [seqs.state_direction] == -1), ...
                           [false true],["Past","Future"])';
g = gTakeLand .* gPastFuture;
[N,G] = histcounts(g);
N = reshape(N,2,3)';
G = reshape(G,2,3)';
takeoff_IX = find(all(contains(G,'Takeoff'),2));
midair_IX = find(all(contains(G,'Mid-air'),2));
landing_IX = find(all(contains(G,'Landing'),2));
new_order = [takeoff_IX midair_IX landing_IX];
N = N(new_order,:);
G = G(new_order,:);
N([1 3],:) = N([1 3],:)./TakeLand_thr;
N([2],:) = N([2],:)./(1-2*TakeLand_thr);
hb=bar(N);
[hb.FaceColor] = disperse(PastFutureClrs');
% hb(1).FaceColor = 1.0.*[1 1 1];
% hb(1).FaceColor = 0.1.*[1 1 1];
hax=gca;
hax.XTick = 1:3;
hax.XTickLabel = {'Takeoff','Mid-air','Landing'};
ylabel('Counts (norm.)')
% create legend
w = 0.08*range(hax.XLim);
h = 0.10*range(hax.YLim);
x = hax.XLim(1)+0.25*range(hax.XLim);
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
text(-0.1,1.2, 'e', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{3}(1))
text(-0.1,1.12, 'f', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{4}(1))
text(-0.25,1.035, 'g', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{5}(1))
text(-0.3,1.2, 'h', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{6}(1))
text(-0.7,1.2, 'i', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{7}(1))
text(-0.2,1.2, 'j', 'Units','normalized','FontWeight','bold','FontSize',font_size);



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
