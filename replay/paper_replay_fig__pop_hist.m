%% Replay Fig: Population histograms
%%
close all
clear 
clc

%% plotting options
epoch_types = {'sleep','rest'};
clrs = [1 0 0; 0 0 1];
PastFutureClrs = [0.1.*[1 1 1]; 1 1 1];

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Fig_replay_pop_hist';
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
panels_size = [3 3];
clear panels_hist
panels_hist(1) = axes('position', [2 20 panels_size]);
panels_hist(2) = axes('position', [7 20 panels_size]);
panels_hist(3) = axes('position', [12 20 panels_size]);
panels_hist(4) = axes('position', [17 20 panels_size]);
panel_over_represent = axes('position', [5 14 10 2.5]);
panel_TakeLandPastFuture = axes('position', [5 10 3 2.5]);
panel_seq_dynamics(1) = axes('position', [5 4 3 3]);
panel_seq_dynamics(2) = axes('position', [11 4 3 3]);

%% load data
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
    'Replay distance (norm)';
    };

%% plot all pop hists
for ii_fn = 1:length(features_names)
    fn = features_names{ii_fn};
    axes(panels_hist(ii_fn));
    cla
    hold on
    for ii_epoch_type = 1:length(epoch_types)
        X = [seqs_all{ii_epoch_type}.(fn)];
        g = gArena{ii_epoch_type};
        histogram(X(g=='200m'),'DisplayStyle','stairs','Normalization','pdf','LineStyle','-', 'EdgeColor',clrs(ii_epoch_type,:));
        histogram(X(g=='120m'),'DisplayStyle','stairs','Normalization','pdf','LineStyle','--', 'EdgeColor',clrs(ii_epoch_type,:));
    end
    xlabel(xlable_strs{ii_fn});
    ylabel('Probability density');
    hax=gca;
    hax.TickLength(1) = [0.025];
end

%% legend (main hists)
panels_hist_legend = axes('position', [10 24 0.5 0.5]);
% axes(panels_hist_legend);
cla
hold on
t = linspace(0,1,100);
x = pulstran(t,linspace(0,1,3),'rectpuls',1/6);
x(x>0) = nan;x(~isnan(x)) = 0;
clear h
plot(  t  ,   x  ,   '-', 'color', clrs(1,:), 'LineWidth',1,'Clipping','off');
plot(  t  ,   x  +1.5, '-', 'color', clrs(2,:), 'LineWidth',1,'Clipping','off');
plot([0 1], [.5 .5],   '-', 'color', clrs(1,:), 'LineWidth',1,'Clipping','off');
plot([0 1], [.5 .5]+1.5, '-', 'color', clrs(2,:), 'LineWidth',1,'Clipping','off');
text(1.3, 0, 'Sleep (120m)','FontSize',7,'HorizontalAlignment','left');
text(1.3, .5, 'Sleep (200m)','FontSize',7,'HorizontalAlignment','left');
text(1.3, 1.5, 'Rest (120m)','FontSize',7,'HorizontalAlignment','left');
text(1.3, 2, 'Rest (200m)','FontSize',7,'HorizontalAlignment','left');
xlim([0 1])
ylim([0 1])
axis off

%% over representations (near balls)
axes(panel_over_represent)
cla
hold on
lw = 1.2;
pos_bins = linspace(0,1,1000);
for ii_epoch_type = 1:length(epoch_types)
    seqs = [seqs_all{ii_epoch_type}];
    xi = [seqs.start_pos_norm; seqs.end_pos_norm]';
    [~,~,~,X] = get_data_in_ti(pos_bins,xi);
    X = [X{:}];
    histogram(X,'DisplayStyle','stairs','EdgeColor',clrs(ii_epoch_type,:),'BinEdges',pos_bins,'LineStyle','-','Normalization','pdf','LineWidth',lw);
end
% X_sleep = [seqs_sleep.middle_pos_norm];
% X_rest = [seqs_rest.middle_pos_norm];
% histogram(X_sleep,'DisplayStyle','stairs','EdgeColor',clrs(1,:),'BinEdges',linspace(0,1,50),'LineStyle','-','Normalization','probability','LineWidth',lw);
% histogram(X_rest,'DisplayStyle','stairs','EdgeColor',clrs(2,:),'BinEdges',linspace(0,1,50),'LineStyle','-','Normalization','probability','LineWidth',lw);
% histogram(X_sleep(gArenaSleep=='200m'),'DisplayStyle','stairs','EdgeColor',clrs(1,:),'BinEdges',linspace(0,1,50),'LineStyle','-','Normalization','pdf');
% histogram(X_sleep(gArenaSleep=='120m'),'DisplayStyle','stairs','EdgeColor',clrs(1,:),'BinEdges',linspace(0,1,50),'LineStyle','--','Normalization','pdf');
% histogram(X_rest(gArenaRest=='200m'),'DisplayStyle','stairs','EdgeColor',clrs(2,:),'BinEdges',linspace(0,1,50),'LineStyle','-','Normalization','pdf');
% histogram(X_rest(gArenaRest=='120m'),'DisplayStyle','stairs','EdgeColor',clrs(2,:),'BinEdges',linspace(0,1,50),'LineStyle','--','Normalization','pdf');
xlabel('Replay position (norm.)');
ylabel('Probability density');
hax = gca;
xx = hax.XLim(1)+[.15 .20].*range(hax.XLim);
yy1 = hax.YLim(1)+[.8 .8].*range(hax.YLim);
yy2 = hax.YLim(1)+[.9 .9].*range(hax.YLim);
plot(xx, yy1, 'LineWidth',lw,'Color',clrs(1,:),'Clipping','off');
plot(xx, yy2, 'LineWidth',lw,'Color',clrs(2,:),'Clipping','off');
text(xx(2)+0.015, yy1(1), 'Sleep', 'HorizontalAlignment','left', 'VerticalAlignment','middle','FontSize',7);
text(xx(2)+0.015, yy2(1), 'Rest', 'HorizontalAlignment','left', 'VerticalAlignment','middle','FontSize',7);

%% Takeoff/Landing X Past/Future
axes(panel_TakeLandPastFuture)
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
clc
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
hax.XTickLabel = {'Takeoff','Midair','Landing'};
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

%% 
for ii_epoch_type = 1:length(epoch_types)
    axes(panel_seq_dynamics(ii_epoch_type));
    cla reset
    hold on
    events2 = events_all{ii_epoch_type};
    invalid = cellfun(@isempty,{events2.prev_event});
    events2(invalid)=[];
    events1 = [events2.prev_event];
    seqs2 = [events2.seq_model];
    seqs1 = [events1.seq_model];
    hold on
    plot([seqs2.middle_pos],[seqs1.middle_pos],'k.','MarkerSize',2,'Color',clrs(ii_epoch_type,:));
    h=refline(1,0);
    h.Color = 0.5*[1 1 1];
    xlabel('Position of replay {\iti} (m)')
    ylabel('Position of replay {\iti+1} (m)')
end

%% print/save the figure
fig_name_out = fullfile(res_dir, sprintf('%s',fig_name_str));

print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
disp('figure was successfully saved to pdf/tiff/fig formats');

%%

