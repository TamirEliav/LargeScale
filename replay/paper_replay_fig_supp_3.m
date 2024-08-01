%% Replay - Fig supp 3 - replay vs flight speed
clear 
clc
close all

%% data options 
params_opt = 11; % decoding opt 
% params_opt = 21; % decoding opt (random walk = fixed speed)
novelty_session_num_thr = 5;
minSeqsThr = 30;

%% plotting options

%% graphics params
line_styles = {'-',':'};
line_widths = [1.1 2];
clrs = {[.6 .1 .8],[.1 .8 .1]};
epoch_type_clrs = {[.6 .1 .8],[.1 .8 .1]};
directions_clrs = {[0    0.3843    0.7451];[ 0.5216    0.2471         0]};
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
y = 20.8;
w = 3.8;
h = 3;
panels{1}(1) = axes('position', [2 y w h]);
panels{1}(2) = axes('position', [3.1 y+2.2 0.5 0.5]);

panels{2}(1) = axes('position', [7 y w h]);
panels{2}(2) = axes('position', [7.3 y+3.2 0.5 0.4]);

panels{3}(1) = axes('position', [12   y w h]);
panels{3}(2) = axes('position', [12.3 y+2.8 0.5 0.5]);

panels{4}(1) = axes('position', [17   y w h]);

x = linspace(2,18,5);
y = [10 6]+6.5;
for ii = 1:length(y)
    for jj = 1:length(x)
        panels{5}(ii,jj) = axes('position', [x(jj) y(ii) 2.8 2.8]);
    end
end

panels{6}(1) = axes('position', [2 2 3.5 8]);
panels{6}(2) = axes('position', [6.5 2 3.5 8]);
panels{7}(1) = axes('position', [12 6.4 4 3.5]);
panels{8}(1) = axes('position', [12 2 4 3.5]);
panels{8}(2) = axes('position', [14 4 .5 .5]);
panels{9}(1) = axes('position', [17.8 6.4 3 3]);

%% load data
if ~exist('events_all','var')
    disp('loading data!!!')
epoch_types = {'sleep','rest'};
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
bats = unique(T.bat_num);
bats_clr_map = containers.Map(num2cell(bats),{'r','g','b','k','m','c',[0.8 0.8 0]});
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
    h(1).MarkerSize = 3;
    h(1).Color = clrs{ii_epoch_type};
    h(2).Color = clrs{ii_epoch_type};
%     h(3).Color = clrs{ii_epoch_type};
%     h(4).Color = clrs{ii_epoch_type};
    delete(h(3))
    delete(h(4))
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


%% plot flight speed per bat
axes(panels{4})
cla reset
hold on
for ii_bat = 1:length(bats)
    bat_num = bats(ii_bat);
    bat_session_IX = find(T.bat_num==bat_num);
    x = flight_speed_all(bat_session_IX);
    xi = linspace(5.5,8.5,100);
    c = bats_clr_map(bat_num);
    f = ksdensity(x,xi,"Bandwidth",0.05);
    plot(xi,f,'DisplayName',"bat"+bat_num,'Color',c,'LineWidth',1.6);
%     plot(xi,f,'DisplayName',sprintf('bat%d(%s)',bat_num,exp.details.recordingArena),'Color',c,'LineWidth',2);
end
xlim([5.5 8.5])
ylim([0 8])
xlabel('Flight speed (m/s)','Units','normalized','Position',[0.5 -0.12])
ylabel('Probability density','Units','normalized','Position',[-0.12 0.5])
hax=gca;
hax.XRuler.TickLength(1) = 0.035;
hax.YRuler.TickLength(1) = 0.035;
hax.XRuler.TickLabelGapOffset = -2;
hax.YRuler.TickLabelGapOffset = 0;
% hl=legend("bat "+bats,'Location','northwest','Box','off')
% hl.Position([1 2]) = hl.Position([1 2]) + [-0.07 0.25].*hl.Position([3 4]);


%% panel e - scatters of compression/duration/distance/speed
scatters_features = {
    'compression','duration';
    'speed','duration';
    'compression','distance';
    'speed','distance';
    'distance','duration';
    };
fn_map = containers.Map();
fn_map('compression') = 'Compression ratio';
fn_map('duration') = 'Replay duration (s)';
fn_map('distance') = 'Replay distance (m)';
fn_map('speed') = 'Replay speed (m/s)';
for ii_epoch_type = 1:length(epoch_types)
    for ii_scatter = 1:size(scatters_features,1)
        axes(panels{5}(ii_epoch_type,ii_scatter))
        cla reset
        hold on
        events = [events_all_per_session{ii_epoch_type}{:}];
        seqs = [events.seq_model];
        fn_x = scatters_features{ii_scatter,2};
        fn_y = scatters_features{ii_scatter,1};
        x = [seqs.(fn_x)];
        y = [seqs.(fn_y)];
        plot(x,y,'.','Color',clrs{ii_epoch_type})
        xlabel(fn_map(fn_x))
        ylabel(fn_map(fn_y))
        hax=gca;
        hax.XRuler.TickLabelGapOffset = -1;
        hax.YRuler.TickLabelGapOffset = 1;
    end
end





%% ========================================================================
%% arrange sessions to load (novelty bats)
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
novelty_exposure_bats = [184 2382];
novelty_exposure_bats_first_session = {'b0184_d191127','b2382_d190623'};
bat_novelty_session_mapping = containers.Map(novelty_exposure_bats,novelty_exposure_bats_first_session);
TF_novelty_bats = ismember(T.bat_num, novelty_exposure_bats);

%% get session number (from exposure)
exp_summary_t = DS_get_exp_summary();
for ii_exp = 1:height(exp_summary_t)
    exp = exp_summary_t(ii_exp,:);
    if any(exp.batNum == novelty_exposure_bats)
        exposure_exp_ID = bat_novelty_session_mapping(exp.batNum);
        exposure_date = exp_summary_t.date(exposure_exp_ID);
        TF = exp_summary_t.batNum ==exp.batNum &...
            exp_summary_t.date >= exposure_date &...
            exp_summary_t.date <= exp.date;
        session_num = sum(TF);
    else
        session_num = nan;
    end
    exp_summary_t{ii_exp,'session_num_from_exposure'} = session_num;
end
T = innerjoin(T,exp_summary_t,'Keys','exp_ID');

%% arrange to early and late sessions
early_days_IX = find(T.session_num_from_exposure <= novelty_session_num_thr);
late_days_IX = find(T.session_num_from_exposure > novelty_session_num_thr);
early_late_IX = {early_days_IX,late_days_IX};
exp_list = T.exp_ID;

%% load data
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
        [events.session_num_from_exposure] = deal(T.session_num_from_exposure(ii_exp));
        [events.recordingArena] = deal(exp.details.recordingArena);
        event_struct = events(1,[]);
        [events(2:end).prev_event] = disperse(events(1:end-1));
        events(1).prev_event = event_struct;
        if ~isfield(events,'rest_ball_num')
            [events.rest_ball_num] = deal(nan);
        end
        events_all_per_session{ii_epoch_type,ii_exp} = events;
    end
end

%% replay directionality bias - calc
directionality_contrast_index = [];
directionality_fraction = [];
directionality_binom_pval = [];
nSeqs = [];
for ii_epoch_type = 1:length(epoch_types)+1
    for ii_exp = 1:length(exp_list)
        switch ii_epoch_type
            case {1,2}
                events = events_all_per_session{ii_epoch_type,ii_exp};
            case 3
                events = [events_all_per_session{:,ii_exp}];
        end
        if isempty(events)
            continue;
        end
        seqs = [events.seq_model];
        n1 = sum([seqs.direction]==1);
        n2 = sum([seqs.direction]==-1);
        directionality_contrast_index(ii_epoch_type,ii_exp) = (n1-n2)/(n1+n2);
        directionality_fraction(ii_epoch_type,ii_exp) = (n1)/(n1+n2);
        binom_pval = myBinomTest(n1,n1+n2,0.5);
        directionality_binom_pval(ii_epoch_type,ii_exp) = binom_pval;
        directionality_binom_surprise(ii_epoch_type,ii_exp) = -log10(binom_pval);
        nSeqs(ii_epoch_type,ii_exp) = length(seqs);
    end
end
directionality_vals = cat(3,directionality_contrast_index,directionality_fraction,directionality_binom_pval,directionality_binom_surprise);
directionality_labels = {'contrast index','fraction','binom pval','binom surprise'};
% bat_sym_map = containers.Map(num2cell([184;2382]),{'+','x'});


%% replay directionality bias - example session
replay_directionality_bias_ex_exp_ID = 'b0184_d191130'
ii_exp = find(strcmp(T.exp_ID,replay_directionality_bias_ex_exp_ID));
example_session_num_form_exposure = T.session_num_from_exposure(replay_directionality_bias_ex_exp_ID);
for ii_epoch_type = 1:length(epoch_types)
    axes(panels{6}(ii_epoch_type))
    cla reset
    hold on
    events = events_all_per_session{ii_epoch_type,ii_exp};
    seqs = [events.seq_model];
    epoch_sep = find(diff([events.epoch_num])~=0)+0.5;
    seqs_edges = [seqs.start_pos_norm; seqs.end_pos_norm];
    seqs_IX = 1:length(seqs);
    if strcmp(epoch_types(ii_epoch_type),'sleep')
        plot([0 1],repmat(epoch_sep,2,1),':','LineWidth',1.5,'Color',0.5*[1 1 1]);
    end
    h=plot(seqs_edges,[seqs_IX;seqs_IX],'-','LineWidth',.55);
    dir_1_IX = [seqs.state_direction]==1;
    dir_2_IX = [seqs.state_direction]==-1;
    [h(dir_1_IX).Color] = disperse(repelem(directions_clrs(1),length(dir_1_IX)));
    [h(dir_2_IX).Color] = disperse(repelem(directions_clrs(2),length(dir_2_IX)));
    scatter(seqs_edges(1,:),seqs_IX, 3, [seqs.state_direction]==-1, "filled");
    hax=gca;
    hax.Colormap = cell2mat(directions_clrs);
    hax.XTick = [0 1];
    hax.YTick = [1 10*ceil(length(seqs)/10)];
    hax.XLim = [0 1];
    hax.YRuler.TickLabelGapOffset = -1;
    hax.XRuler.TickLabelGapOffset = 1;
    xlabel('Position (norm.)', 'Units','normalized', 'Position',[0.5 -0.02]);
    ylabel('Replay event no.', 'Units','normalized', 'Position',[-0.1 0.5]);
    switch epoch_types{ii_epoch_type}
        case 'sleep'
            title_str = 'Sleep replays';
        case 'rest'
            title_str = 'Awake replays';
    end
    title(title_str, 'Units','normalized', 'Position',[0.5 0.99],'FontWeight','normal');
end
text(-.2, 1.07, 'Example session: replay directionality over-representation','FontSize',9,'HorizontalAlignment','center','Units','normalized')


%% replay directionality bias - plot trend over exposure to enviroenment
axes(panels{7}(1));
cla reset
hold on
clrs = [epoch_type_clrs,'k'];
for ii_epoch_type = 1:length(epoch_types)+1
    c = clrs{ii_epoch_type};
%         sym = arrayfun(@(bat_num)bat_sym_map(bat_num),T.bat_num,'UniformOutput',false);
    x = T.session_num_from_exposure;
    y = directionality_contrast_index(ii_epoch_type,:);
    y(nSeqs(ii_epoch_type,:)<minSeqsThr) = nan;
    g = findgroups(T.bat_num);
    hs = gscatter(x,y,g,c,'ox',4,false);
    nPointsSmooth = 3;
    k = (nPointsSmooth-1)/2;
    xi = [-1 1].*k + [(1-k):(max(x)+k)]';
    xx=[];
    yy=[];
    for ii_xi = 1:size(xi,1)
        TF = x>=xi(ii_xi,1) & x<=xi(ii_xi,2);
        xx(ii_xi) = nanmedian(x(TF));
        yy(ii_xi) = nanmedian(y(TF));
    end
    plot(xx,yy,'-','Color',c);
end
xlabel('Session no.');
ylabel('Replay directionality')
ylim([-1 1])
xticks([1 40])
text(5.5,-0.75,"\leftarrow"+"session #"+example_session_num_form_exposure)
% h = annotation('arrow');
% h.Parent = gca;
% h.Units = 'normalized';
% h.X = [.8 .5];
% h.Y = [1 1].*0.3;

%% Inter-replay-interval (pooled across bats)
axes(panels{8}(1));
cla reset
hold on
IRI_bin_size = 0.1;
IRI_bin_limits = [0 20];
IRI = load('E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\internal_review_figures\inter_replay_interval_all.mat');
hax=gca;
% hax.ColorOrder = epoch_type_clrs;
for ii_epoch_type = 1:length(epoch_types)
    x = [IRI.IRI_all{ii_epoch_type,:}];
    hh=histogram(x,'Normalization','pdf','DisplayStyle','Stairs','BinWidth',IRI_bin_size,'EdgeColor',epoch_type_clrs{ii_epoch_type});
    hh.LineWidth = 1;
end
% x = [IRI.IRI_all{:,:}];
% hh=histogram(x,'Normalization','pdf','DisplayStyle','Stairs','BinWidth',IRI_bin_size,'EdgeColor','k');
% hh.LineWidth = 2;
xlim(IRI_bin_limits)
xlabel('Inter-replay-interval (s)')
ylabel('Probability density')

%% add legend
axes(panels{8}(2));
cla reset
hold on
axis off
x = [0 1];
y = [0 0];
plot(x,y+1,'Color',epoch_type_clrs{1})
plot(x,y+2,'Color',epoch_type_clrs{2})
% plot(x,y+3,'Color','k')
x = 1.2;
text(x,1,'Sleep')
text(x,2,'Awake')
% text(x,3,'Pooled')
axis ij

%% scatters: number of ripples vs replay duration
load("L:\processed_data_structs\replay_events.mat");
events = [events_all_per_session{:}];
seqs = [events.seq_model];

% features_fn = {'duration','compression','distance','distance_norm'};
features_fn = {'duration'};
fn_map = containers.Map();
fn_map('duration') = 'Replay duration (s)';
fn_map('compression') = 'Compression ratio';
fn_map('distance') = 'Replay distance (m)';
fn_map('distance_norm') = 'Replay distance (norm)';
nFeatures = length(features_fn);
for ii_fn = 1:nFeatures
    axes(panels{9}(ii_fn));
    cla reset
    hold on
    fn = features_fn{ii_fn};
    fn_label = fn_map(fn);
    X = [seqs.(fn)];
    Y = [events.num_ripples];
    Y(Y==0)=nan;
    jitter_sd = .1;
    rng(0);
    Y_jitter = randn(size(Y)).*jitter_sd;
    plot(X,Y+Y_jitter,'.k')
%     swarmchart(Y,X,'.k','XJitter','density','XJitterWidth',.5);
    [r,pval] = corr(X',Y','type','Spearman','rows','pairwise');
    fprintf('%s: r=%.2g pval=%.2g\n',fn,r,pval);
%     view([90 -90])
%     yticks(1:6)
    xlabel(fn_label)
    ylabel('No. ripples')
end


%% add panel letters
font_size = 11;
axes(panels{1}(1))
text(-0.3,1.1, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{2}(1))
text(-0.3,1.1, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{3}(1))
text(-0.3,1.1, 'c', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{4}(1))
text(-0.3,1.1, 'd', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{5}(1,1))
text(-0.4,1.2, 'e', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{6}(1))
text(-0.3,1.1, 'f', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{7}(1))
text(-0.3,1.1, 'g', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{8}(1))
text(-0.3,1.1, 'h', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{9}(1))
text(-0.3,1.1, 'i', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = sprintf('%s_decoding_opt_%d',fig_name_str, params_opt);
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%
