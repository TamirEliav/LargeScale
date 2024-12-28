%% Replay - Fig supp 4 - per bat results
%%
clear 
clc
close all

%% data options 
params_opt = 11; % decoding opt 
% params_opt = 21; % decoding opt (random walk = fixed speed)

%% plotting options
features_names = {'duration';'compression';'distance';'distance_norm';};
xlable_strs = {
    {'Replay';'duration (s)'};
    {'Compression ratio';'(replay speed / flight speed)'};
    {'Replay';'distance (m)'};
    {'Replay';'distance (norm.)'};
    };
nFeatures = length(features_names);
nbins = 50; % spatial bins for coverage

%% graphics params
PastFutureClrs = [0.1.*[1 1 1]; 1 1 1];
directions_clrs = {[0    0.3843    0.7451];[ 0.5216    0.2471         0]};
panels_xlim = [-0.0950    1.9950; -1.9500   40.9500; -0.4500   56.5500; -0.0210    0.4410];
panels_xticks = {[0 0.5 1 1.5],[0 20 40],[0:25:50],[0 0.2 0.4]};

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Figure S4';
fig_caption_str = ' per bat results';
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

%% get number of bats
exp_t = DS_get_exp_summary();
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
T = innerjoin(T,exp_t,'Keys','exp_ID','RightVariables','recordingArena');
clear exp_t
T2 = groupsummary(T,{'bat_num','recordingArena'});
T2 = sortrows(T2,{'recordingArena','bat_num'})
bats = unique(T2.bat_num,'stable');
nBats = length(bats);
for ii=1:height(T2)
    if T2.recordingArena(ii)=="120m"
        T2.recordingArena{ii}='130m';
    end
end

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
clear panels
xoffset = 2;
panels_hist_size = [2 1.5];
panles_hist_pos_x = linspace(1.3, 9.8, nFeatures)+xoffset;
panles_hist_pos_y = linspace(24.5, 12, nBats);
panles_hist_pos_y(end+1) = panles_hist_pos_y(end)-4;
panles_hist_pos_y = panles_hist_pos_y - 4;
for ii_bat=1:length(panles_hist_pos_y)
    panel_pos(2) = panles_hist_pos_y(ii_bat);
    for ii_feature=1:nFeatures
        panel_pos(1) = panles_hist_pos_x(ii_feature);
        panels{1}(ii_bat,ii_feature) = axes('position', [panel_pos panels_hist_size]);
    end
    panels{2}(ii_bat) = axes('position', [12.4+xoffset panel_pos(2) 7 1.5]);
end

%%
% panels{1}(1) = axes('position', [3 23 3 3]);
% panels{1}(2) = axes('position', [7 23 3 3]);
% panels{1}(3) = axes('position', [11 23 3 3]);
% panels{1}(4) = axes('position', [15 23 3 3]);
% 
% panels{2}(1) = axes('position', [3 18.5 8 2]);
% panels{2}(2) = axes('position', [3 16.0 8 2]);
% panels{2}(3) = axes('position', [3 13.5 8 2]);
% panels{2}(4) = axes('position', [3 11.0 8 2]);
% 
% panels{3}(1) = axes('position', [3 8 8 2]);
% 
% % panels{4}(1) = axes('position', [12.5 15 8 6]);
% % panels{4}(2) = axes('position', [12.5  8 8 6]);
% 
% panels{4}(1) = axes('position', [12.5 8 4 13]);
% panels{4}(2) = axes('position', [17   8 4 13]);
% 
% panels{5}(1) = axes('position', [3 4.5 2 2]);
% 
% panels{6}(1) = axes('position', [6.5 4.5 3 2]);

%% load data
if ~exist('events_all','var')
epoch_types = {'sleep','rest'};
events_all_per_session = {};
for ii_epoch_type = 1:length(epoch_types)
    for ii_exp = 1:length(T.exp_ID)
        exp_ID = T.exp_ID{ii_exp};
        exp = exp_load_data(exp_ID,'details');
        epoch_type = epoch_types{ii_epoch_type};
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
events_all = {[events_all_per_session{1}{:}], [events_all_per_session{2}{:}]};
seqs_all = {[events_all{1}.seq_model], [events_all{2}.seq_model]};
gArena = {categorical({events_all{1}.recordingArena}), categorical({events_all{2}.recordingArena}) };
end

%% panels A-D: features hists per bat
axes(panels{1}(1))
hold on
clrs = {[.6 .1 .8],[.1 .8 .1]};
for ii_fn = 1:nFeatures
    for ii_bat = 1:nBats
        fn = features_names{ii_fn};
        axes(panels{1}(ii_bat,ii_fn));
        cla
        hold on
        lw = 1.1;
        for ii_epoch_type = 1:length(epoch_types)
            bat_exp_IX = T.bat_num == T2.bat_num(ii_bat);
            bat_events = [events_all_per_session{ii_epoch_type}{bat_exp_IX}];
            bat_seqs = [bat_events.seq_model];
            X = [bat_seqs.(fn)];
            histogram(X,'DisplayStyle','stairs','Normalization','pdf','LineStyle','-', 'EdgeColor',clrs{ii_epoch_type},'LineWidth',lw);
        end
        hax=gca;
        hax.TickLength(1) = [0.035];
        hax.XRuler.TickLabelGapOffset = -1;
        hax.YRuler.TickLabelGapOffset = 1;
        hax.XLim = panels_xlim(ii_fn,:);
        hax.XTick = panels_xticks{ii_fn};
        hax.XTickLabelRotation = 0;
        if ii_bat == nBats
            xlabel_pos = [0.45 0.5 0.55 0.5];
            xlabel(xlable_strs{ii_fn}, 'Units','normalized', 'Position',[xlabel_pos(ii_fn) -0.3]);
        end
        if ii_fn == 1
            ylabel({'Probability';'density'}, 'Units','normalized', 'Position',[-0.2 .5]);
        end
%         if ii_fn == 1
%             bat_arena = unique(T.recordingArena(bat_exp_IX));
%             text(-.6, .5, {"bat "+T2.bat_num(ii_bat);"("+bat_arena+")"},'FontSize',7,'HorizontalAlignment','left','Units','normalized','FontWeight','bold');
%         end
    end
    linkaxes(panels{1}(:,ii_fn),'x');
end

%% legend (features hists)
if exist('panels_hist_legend','var')
    delete(panels_hist_legend);
end
panels_hist_legend = axes('position', [xoffset+2.2 21.5 0.2 0.3]);
cla
hold on
plot([0 1], [1 1], 'color', clrs{1}, 'LineWidth',lw,'Clipping','off');
plot([0 1], [0 0], 'color', clrs{2}, 'LineWidth',lw,'Clipping','off');
text(1.5, 1.0, 'Sleep','FontSize',7,'HorizontalAlignment','left');
text(1.5, 0.0, 'Awake','FontSize',7,'HorizontalAlignment','left');
xlim([0 1])
ylim([0 1])
axis off

%% panel E - coverage per session per bat
% coverage = load('E:\Tamir\work\PROJECTS\LargeScale\paper_replay\data_prepared_for_figures\replay_coverage.mat');
coverage = decoding_calc_coverage_over_exp_table(T,params_opt,nbins);
for ii_bat = 1:nBats
    bat_exp_IX = coverage.T.bat_num == T2.bat_num(ii_bat);
    T2.nEvents(ii_bat) = sum(coverage.n_seqs_all(bat_exp_IX,:,:),'all');
end
T2 = sortrows(T2,{'nEvents'},'descend')
for ii_bat = 1:nBats
    axes(panels{2}(ii_bat))
    cla
    hold on
    bat_exp_IX = coverage.T.bat_num == T2.bat_num(ii_bat);
    c = coverage.coverage_all(bat_exp_IX,:,:,:);
    c = normalize(c,4,'norm',1);
    nBins = size(c,4);
    x = linspace(0,1,nBins);
    lw = .2;
    plot(x, squeeze(c(:,1,1,:)),'-','LineWidth',lw,'Color',directions_clrs{1});
    plot(x, squeeze(c(:,2,1,:)),'--','LineWidth',lw,'Color',directions_clrs{1});
    plot(x, squeeze(c(:,1,2,:)),'-','LineWidth',lw,'Color',directions_clrs{2});
    plot(x, squeeze(c(:,2,2,:)),'--','LineWidth',lw,'Color',directions_clrs{2});
    hax=gca;
    hax.XTick = [0 1];
    hax.XLim = [0 1];
    hax.XRuler.TickLabelGapOffset = -1.5;
    hax.YRuler.TickLabelGapOffset = 1;
%     if ii_bat == 1
%         hl=legend({'sleep dir 1','rest dir 1','sleep dir 2','rest dir 2'},'NumColumns',2);
%         hl.Units = 'normalized';
%         hl.Position([1 2]) = [0.22 0.715];
%         hl.Position([3 4]) = [0.1 0.025];
%         hl.Box = 'off';
%     end
    if ii_bat == nBats
        xlabel('Position (norm.)', 'Units','normalized', 'Position',[0.5 -0.11]);
%         ylabel({'Replay coverage';'(prob.)'}, 'Units','normalized', 'Position',[-0.1 0.5]);
    end
    text(-1.9, .5, { ...
        "bat "+T2.bat_num(ii_bat); ...
        T2.recordingArena{ii_bat}; ...
        "n_{sessions} = " + T2.GroupCount(ii_bat); ...
        "n_{events} = " + T2.nEvents(ii_bat); ...
        },'FontSize',7,'HorizontalAlignment','center','Units','normalized','FontWeight','bold');
end

%% panel E legend
if exist('panels_coverage_legend','var')
    delete(panels_coverage_legend);
end
panels_coverage_legend(1) = axes('position', [panels{2}(1).Position([1 2])+[1 1.25] 0.5 0.5]);
panels_coverage_legend(2) = axes('position', [panels{2}(1).Position([1 2])+[4 1.25] 0.5 0.5]);
arrow_str = {'\rightarrow','\leftarrow'};
for ii_dir = 1:2
    axes(panels_coverage_legend(ii_dir))
    cla
    hold on
    t = linspace(0,1,100);
    x = pulstran(t,linspace(0,1,3),'rectpuls',1/6);
    x(x>0) = nan;x(~isnan(x)) = 0;
    clear h
    plot(  t  ,   x,     'color', directions_clrs{ii_dir}, 'LineWidth',lw,'Clipping','off','LineWidth',1.5);
    plot([0 1], [.5 .5], 'color', directions_clrs{ii_dir}, 'LineWidth',lw,'Clipping','off','LineWidth',1.5);
    text(1.3, 0.5, "Sleep dir "+ii_dir,'FontSize',7,'HorizontalAlignment','left');
    text(1.3, 0, "Awake dir "+ii_dir,'FontSize',7,'HorizontalAlignment','left');
    text(0.2, -0.3, arrow_str{ii_dir},'Color',directions_clrs{ii_dir},'FontWeight','bold', 'FontSize',10);
    xlim([0 1])
    ylim([0 1])
    axis off
end




%% ========================================================================
% 2 bats

%% 2 bats data - load data
exp_list_2bats = {
%     'b2299_d191202', % <70% accuracy, excluded
%     'b2299_d191203', % <70% accuracy, excluded
    'b2299_d191204',
    'b2299_d191205',
    'b2299_d191208',
    'b2299_d191209',
%     'b2299_d191210', % <70% accuracy, excluded
    'b2299_d191213',
    };

seqs_all_per_session_2bats = {};
for ii_epoch_type = 1:length(epoch_types)
    for ii_exp = 1:length(exp_list_2bats)
        exp_ID = exp_list{ii_exp};
        exp = exp_load_data(exp_ID,'details');
        epoch_type = epoch_types{ii_epoch_type};
        params_opt = 11;
        event_type = 'posterior';
        [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
        seqs = decoding_apply_seq_inclusion_criteria([events.seq_model]);
        seqs_all_per_session_2bats{ii_epoch_type,ii_exp} = seqs;
    end
end

%% panels F-I: features hist (2bats data)
axes(panels{1}(end))
hold on
clrs = {[.6 .1 .8],[.1 .8 .1]};
bat_seqs = [seqs_all_per_session_2bats{:}];
for ii_fn = 1:nFeatures
    fn = features_names{ii_fn};
    axes(panels{1}(end,ii_fn));
    cla
    hold on
    lw = 1.1;
    X = [bat_seqs.(fn)];
    histogram(X,'DisplayStyle','stairs','Normalization','pdf','LineStyle','-', 'EdgeColor','k','LineWidth',lw);
    hax=gca;
    hax.TickLength(1) = [0.035];
    hax.XRuler.TickLabelGapOffset = -1;
    hax.YRuler.TickLabelGapOffset = 1;
    hax.XLim = panels_xlim(ii_fn,:);
    hax.XTick = panels_xticks{ii_fn};
    hax.XTickLabelRotation = 0;
    xlabel_pos = [0.45 0.5 0.55 0.5];
    xlabel(xlable_strs{ii_fn}, 'Units','normalized', 'Position',[xlabel_pos(ii_fn) -0.3]);
    if ii_fn == 1
        ylabel({'Probability';'density'}, 'Units','normalized', 'Position',[-0.2 .5]);
    end
    linkaxes(panels{1}(:,ii_fn),'x');
end
axes(panels{1}(end,1));
text(-1.1, .5, { ...
    "bat 2299"; ...
    "135m, 2 bats"; ...
    "n_{sessions} = " + length(exp_list_2bats); ...
    "n_{events} = " + length(bat_seqs); ...
    },'FontSize',7,'HorizontalAlignment','center','Units','normalized','FontWeight','bold');


%% panel J - coverage per session 2 bats)
T_2bats = table();
T_2bats.exp_ID = exp_list_2bats;
coverage = decoding_calc_coverage_over_exp_table(T,params_opt,nbins);
axes(panels{2}(end))
cla
hold on
c = coverage.coverage_all;
c = normalize(c,4,'norm',1);
nBins = size(c,4);
x = linspace(0,1,nBins);
lw = .2;
plot(x, squeeze(c(:,1,1,:)),'-','LineWidth',lw,'Color',directions_clrs{1});
plot(x, squeeze(c(:,2,1,:)),'--','LineWidth',lw,'Color',directions_clrs{1});
plot(x, squeeze(c(:,1,2,:)),'-','LineWidth',lw,'Color',directions_clrs{2});
plot(x, squeeze(c(:,2,2,:)),'--','LineWidth',lw,'Color',directions_clrs{2});
hax=gca;
hax.XTick = [0 1];
hax.XLim = [0 1];
hax.XRuler.TickLabelGapOffset = -1.5;
hax.YRuler.TickLabelGapOffset = 1;
xlabel('Position (norm.)', 'Units','normalized', 'Position',[0.5 -0.11]);

%% add panel letters
font_size = 11;
axes(panels{1}(1,1))
text(-0.3,1.3, 'A', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(1,2))
text(-0.3,1.3, 'B', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(1,3))
text(-0.3,1.3, 'C', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(1,4))
text(-0.3,1.3, 'D', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{2}(1))
text(-0.1,1.3, 'E', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{1}(end,1))
text(-0.3,1.3, 'F', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(end,2))
text(-0.3,1.3, 'G', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(end,3))
text(-0.3,1.3, 'H', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(end,4))
text(-0.3,1.3, 'I', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{2}(end))
text(-0.1,1.3, 'J', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = sprintf('%s_decoding_opt_%d',fig_name_str, params_opt);
file_out = fullfile(res_dir, fig_name);
print(fig, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%
