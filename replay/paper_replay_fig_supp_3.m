%% Replay - Fig supp 3 - replay vs flight speed
clear 
clc
close all

%% data options 
params_opt = 11; % decoding opt 
% params_opt = 21; % decoding opt (random walk = fixed speed)

%% plotting options

%% graphics params
line_styles = {'-',':'};
line_widths = [1.1 2];
clrs = {[.6 .1 .8],[.1 .8 .1]};
lw = 1.1;

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Extended_Data_Fig_3';
fig_caption_str = 'replay and flight speed histograms';
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
panels{1}(1) = axes('position', [3 15 5 4]);
panels{1}(2) = axes('position', [5.3 17.8 0.5 0.5]);

panels{2}(1) = axes('position', [10 15 5 4]);
panels{2}(2) = axes('position', [15.7 18 0.5 0.5]);

panels{3}(1) = axes('position', [10   9 5 4]);
panels{3}(2) = axes('position', [10.8 11.8 0.5 0.5]);



%% load data
if ~exist('events_all','var')
    disp('loading data!!!')
epoch_types = {'sleep','rest'};
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
events_all_per_session = {};
speed_mean_per_session = nan(length(epoch_types),height(T));
speed_median_per_session = nan(length(epoch_types),height(T));
speed_std_per_session = nan(length(epoch_types),height(T));
flight_speed_all = [];
exp_arena_all = {};
for ii_exp = 1:height(T)
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details','flight');
    flight_speed_all(ii_exp,:) = [exp.flight.speed_traj.vel_median_median];
    exp_arena_all{ii_exp} = exp.details.recordingArena;
    for ii_epoch_type = 1:length(epoch_types)
        epoch_type = epoch_types{ii_epoch_type};
        event_type = 'posterior';
        [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
        [seqs, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
        events(~TF) = [];
        if isempty(events)
            continue;
        end
        seqs_map_IX = [seqs.state_direction];
        seqs_map_IX(seqs_map_IX == -1) = 2;
        [seqs.session_flight_speed] = disperse(flight_speed_all(ii_exp,seqs_map_IX));
        [events.seq_model] = disperse(seqs);
        [events.recordingArena] = deal(exp.details.recordingArena);
        event_struct = events(1,[]);
        [events(2:end).prev_event] = disperse(events(1:end-1));
        events(1).prev_event = event_struct;
        events_all_per_session{ii_epoch_type}{ii_exp} = events;
        speed_mean_per_session(ii_epoch_type,ii_exp) = mean([seqs.speed]);
        speed_median_per_session(ii_epoch_type,ii_exp) = median([seqs.speed]);
        speed_std_per_session(ii_epoch_type,ii_exp) = std([seqs.speed]);
    end
end
flight_speed_all = mean(abs(flight_speed_all),2);
events_all = {[events_all_per_session{1}{:}], [events_all_per_session{2}{:}]};
seqs_all = {[events_all{1}.seq_model], [events_all{2}.seq_model]};
gArena = {categorical({events_all{1}.recordingArena}), categorical({events_all{2}.recordingArena}) };
end

%% tests (speed correlations: replay vs flight)
% figure
% tiledlayout('flow')
% nexttile
% hold on
% plot(flight_speed_all, speed_mean_per_session(1,:), 'ob');
% plot(flight_speed_all, speed_mean_per_session(2,:), 'or');
% xlim([3 10])
% nexttile
% hold on
% plot(flight_speed_all, speed_median_per_session(1,:), 'ob');
% plot(flight_speed_all, speed_median_per_session(2,:), 'or');
% xlim([3 10])
% 
% clc
% [r p] = corr(flight_speed_all, speed_mean_per_session(1,:)','rows','pairwise','tail','right','type','Pearson')
% [r p] = corr(flight_speed_all, speed_median_per_session(1,:)','rows','pairwise','tail','right','type','Pearson')
% [r p] = corr(flight_speed_all, speed_mean_per_session(2,:)','rows','pairwise','tail','right','type','Pearson')
% [r p] = corr(flight_speed_all, speed_median_per_session(2,:)','rows','pairwise','tail','right','type','Pearson')
% 
% xxx = [flight_speed_all flight_speed_all]';
% yyy = speed_mean_per_session;
% [r p] = corr(xxx(:), yyy(:),'rows','pairwise','tail','right','type','Pearson')
% 
% figure
% plot(xxx(:), yyy(:), 'o')
% 
% figure
% plot(flight_speed_all, speed_mean_per_session(2,:), 'or');


%% panel A - replay speed
axes(panels{1}(1))
cla
hold on
for ii_epoch_type = 1:length(epoch_types)
    X = [seqs_all{ii_epoch_type}.speed];
    g = gArena{ii_epoch_type};
    histogram(X(g=='200m'),'DisplayStyle','stairs','Normalization','pdf','LineStyle',line_styles{1}, 'EdgeColor',clrs{ii_epoch_type},'LineWidth',line_widths(1));
    histogram(X(g=='120m'),'DisplayStyle','stairs','Normalization','pdf','LineStyle',line_styles{2}, 'EdgeColor',clrs{ii_epoch_type},'LineWidth',line_widths(2));
    [~,P,KSSTAT] = kstest2(X(g=='200m'),X(g=='120m'))
    median(X(g=='200m'))
    median(X(g=='120m'))
end
xlabel('Replay speed (m/s)', 'Units','normalized', 'Position',[0.5 -0.12]);
ylabel('Probability density', 'Units','normalized', 'Position',[-0.2 .5]);
hax=gca;
hax.TickLength(1) = [0.016];
hax.XRuler.TickLabelGapOffset = -1;

%% legend (replay speed hists)
axes(panels{1}(2))
cla
hold on
t = linspace(0,1,100);
x = pulstran(t,linspace(0,1,3),'rectpuls',1/6);
x(x>0) = nan;x(~isnan(x)) = 0;
clear h
% plot(  t  ,   x  ,       'color', clrs{2}, 'LineWidth',lw,'Clipping','off');
% plot(  t  ,   x  +1.5,   'color', clrs{1}, 'LineWidth',lw,'Clipping','off');
plot([0 1], [0 0],    ':', 'color', clrs{2}, 'LineWidth',line_widths(2),'Clipping','off');
plot([0 1], [0 0]+1.5,':', 'color', clrs{1}, 'LineWidth',line_widths(2),'Clipping','off');
plot([0 1], [.5 .5],     'color', clrs{2}, 'LineWidth',line_widths(1),'Clipping','off');
plot([0 1], [.5 .5]+1.5, 'color', clrs{1}, 'LineWidth',line_widths(1),'Clipping','off');
text(1.3, 2.0, 'Sleep (200m)','FontSize',7,'HorizontalAlignment','left');
text(1.3, 1.5, 'Sleep (130m)','FontSize',7,'HorizontalAlignment','left');
text(1.3, 0.5, 'Awake (200m)','FontSize',7,'HorizontalAlignment','left');
text(1.3, 0.0, 'Awake (130m)','FontSize',7,'HorizontalAlignment','left');
xlim([0 1])
ylim([0 1])
axis off


%% panel C - Flight speed
axes(panels{3}(1))
cla
hold on
X = flight_speed_all;
g=categorical(exp_arena_all);
histogram(X(g=='200m'),'DisplayStyle','stairs','Normalization','pdf','LineStyle',line_styles{1}, 'EdgeColor','k','LineWidth',line_widths(1),'NumBins',7);
histogram(X(g=='120m'),'DisplayStyle','stairs','Normalization','pdf','LineStyle',line_styles{2}, 'EdgeColor','k','LineWidth',line_widths(2),'NumBins',7);
xlabel('Flight speed (m/s)', 'Units','normalized', 'Position',[0.5 -0.12]);
ylabel('Probability density', 'Units','normalized', 'Position',[-0.15 .5]);
hax=gca;
hax.XLim = [0 10];
hax.XTick = [0 5 10];
hax.TickLength(1) = [0.016];
hax.XRuler.TickLabelGapOffset = -1;
[~,P,CI,STATS] = ttest2(X(g=='200m'), X(g=='120m'), 'Tail','both');
text(0.6,1.1,"{\itP} = "+ sprintf('%.2g',P), 'Units','normalized','FontSize',7);

%% legend (flight speed hists)
axes(panels{3}(2))
cla
hold on
t = linspace(0,1,100);
x = pulstran(t,linspace(0,1,3),'rectpuls',1/6);
x(x>0) = nan;x(~isnan(x)) = 0;
clear h
% plot(  t  ,   x  ,       'color', 'k', 'LineWidth',lw,'Clipping','off');
plot([0 1], [0 0],':',   'color', 'k', 'LineWidth',line_widths(2),'Clipping','off');
plot([0 1], [.5 .5],     'color', 'k', 'LineWidth',line_widths(1),'Clipping','off');
text(1.3, 0.5, '200m','FontSize',7,'HorizontalAlignment','left');
text(1.3, 0.0, '130m','FontSize',7,'HorizontalAlignment','left');
xlim([0 1])
ylim([0 1])
axis off

%% panel B - speed correlations (replay vs flight)
axes(panels{2}(1))
cla
hold on
r=[];
p=[];
for ii_epoch_type = 1:length(epoch_types)
    X = flight_speed_all;
    Y = speed_mean_per_session(ii_epoch_type,:)';
    lm = fitlm(X,Y)
    h=plot(lm);
    h(1).Marker = 'o';
    h(1).MarkerSize = 5;
    h(1).Color = clrs{ii_epoch_type};
    h(2).Color = clrs{ii_epoch_type};
    h(3).Color = clrs{ii_epoch_type};
    h(4).Color = clrs{ii_epoch_type};
    legend off
%     plot(X,Y,'o','Color',clrs{ii_epoch_type},'MarkerSize',5);
    [r(ii_epoch_type) p(ii_epoch_type)] = corr(X, Y,'rows','pairwise','tail','right','type','Pearson');
%     text(0.1,0.9, "{\itr} = "+ sprintf('%.2f',r), 'Units','normalized','FontSize',7);
%     text(0.1,0.8, "{\itP} = "+ sprintf('%.3f',p), 'Units','normalized','FontSize',7,'FontWeight','normal');
end
title("")
xlabel('Flight speed (m/s)', 'Units','normalized', 'Position',[0.5 -0.12]);
ylabel('Replay speed (m/s)', 'Units','normalized', 'Position',[-0.15 .5]);
hax=gca;
hax.TickLength(1) = [0.016];
hax.XRuler.TickLabelGapOffset = -1;

%% legend (speed correlations)
axes(panels{2}(2))
cla
hold on
plot(0,1,'o', 'color', clrs{1}, 'MarkerSize',4,'Clipping','off');
plot(0,0,'o', 'color', clrs{2}, 'MarkerSize',4,'Clipping','off');
text(0.4, 1, 'Sleep'+":  {\itr} = " + sprintf('%.2f',r(1)) + ",  {\itP} = "+ sprintf('%.1g',p(1)),'FontSize',7,'HorizontalAlignment','left');
text(0.4, 0, 'Awake'+":   {\itr} = " + sprintf('%.2f',r(2)) + ",  {\itP} = "+ sprintf('%.2f',p(2)),'FontSize',7,'HorizontalAlignment','left');
xlim([0 1])
ylim([0 1])
axis off

%% add panel letters
font_size = 11;
axes(panels{1}(1))
text(-0.3,1.1, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{2}(1))
text(-0.25,1.1, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{3}(1))
text(-0.25,1.1, 'c', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = sprintf('%s_decoding_opt_%d',fig_name_str, params_opt);
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%
