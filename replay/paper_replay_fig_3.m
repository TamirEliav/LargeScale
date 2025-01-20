%% Replay - Fig 5 - 2 bats crossovers + replay directionality (behavioral relevance)
%%
clear 
clc
close all

%% plotting options
params_opt = 11; % decoding opt 
novelty_session_num_thr = 5;

behavior_ex_opt = 2 % number 2 was chosen

replay_ex_opt = 1;
replay_examples_list = {
    {'rest',11,'b2299_d191213',20}
    };
replay_examples_list = cellfun(@(c)cell2struct(c,{'epoch_type','params_opt','exp_ID','event_num'},2), replay_examples_list)

%% graphics params
epoch_types = {'sleep','rest'};
epoch_type_clrs = {[.6 .1 .8],[.1 .8 .1]};
directions_clrs = {[0    0.3843    0.7451];[ 0.5216    0.2471         0]};

% timediff_max = inf;
timediff_max = 100;
exp_types = {'1 bat','2 bats'};
panels_xlim = [-0.0950    1.9950; -1.9500   40.9500; -0.4500   56.5500; -0.0210    0.4410];
panels_xticks = {[0 0.5 1 1.5],[0 20 40],[0:10:50],[0 0.2 0.4]};
smooth_method = 'movmean';
smooth_n_points = 19;

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Figure_5';
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
panels{1}(1,1) = axes('position', [2 15 3.5 8]);
panels{1}(1,2) = axes('position', [6.5 15 3.5 8]);
panels{2}(1) = axes('position', [12 19.4 4 3.5]);
% panels{2}(2) = axes('position', [12 15 4 3]);
panels{3}(1) = axes('position', [11.5 15 4 3]);
panels{3}(2) = axes('position', [11.5 17.6 .2 .3]);
panels{4}(1) = axes('position', [2 9.5 5 3]);
panels{4}(2) = axes('position', [2 12.85 5 0.4]);
panels{5}(1) = axes('position', [8.8 9.5 3 3]);
panels{5}(2) = axes('position', [8.8 12.5 3 .5]);
panels{6}(1,1) = axes('position', [2 4.5 3 3]);
panels{6}(1,2) = axes('position', [2+2 4.5+3 1 0.7]);
panels{6}(2,1) = axes('position', [6.5 4.5 3 3]);
panels{6}(2,2) = axes('position', [6.5+2 4.5+3 1 0.7]);

%% ========================================================================
%% arrange sessions to load (novelty bats)
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
% novelty_exposure_bats = [184 2382];
bats = unique(T.bat_num);
bats_colors = [1 0 0; 0 1 0; 0 0 1; 0 0 0; 1 0 1; 0 1 1; 0.8 0.8 0];
bats_clr_map = containers.Map(num2cell(bats),num2cell(bats_colors,2));
% novelty_exposure_bats_first_session = {'b0184_d191127','b2382_d190623'};
% bats_list = [34 148 194 2289 9861 184 2382];
% num_training_sessions = [6 11 20 13 11 0 0];
% first_rec_sessions = {'b0034_d180227', 'b0148_d170605', 'b0194_d180429', 'b2289_d180514', 'b9861_d180519', 'b0184_d191127','b2382_d190623'};
% % bat_novelty_session_mapping = containers.Map(novelty_exposure_bats,novelty_exposure_bats_first_session);
% bat_num_training_sessions_mapping = containers.Map(bats_list, num_training_sessions);
% bat_first_rec_exp_ID_mapping = containers.Map(bats_list, first_rec_sessions);
% TF_novelty_bats = ismember(T.bat_num, novelty_exposure_bats);

%% get session number (from exposure) - old
% exp_summary_t = DS_get_exp_summary();
% exp_summary_t(~ismember(exp_summary_t.batNum, bats_list),:)=[];
% for ii_exp = 1:height(exp_summary_t)
%     exp = exp_summary_t(ii_exp,:);
%     num_training_sessions = bat_num_training_sessions_mapping(exp.batNum);
%     first_rec_exp_ID = bat_first_rec_exp_ID_mapping(exp.batNum);
%     first_rec_date = exp_summary_t.date(first_rec_exp_ID);
%     TF = exp_summary_t.batNum ==exp.batNum &...
%         exp_summary_t.date >= first_rec_date &...
%         exp_summary_t.date <= exp.date;
%     session_num = sum(TF) + num_training_sessions;
%     exp_summary_t{ii_exp,'session_num_from_exposure'} = session_num;
% end
% T = innerjoin(T,exp_summary_t,'Keys','exp_ID');

%% get session number (from exposure) - corrected
exp_summary_t = DS_get_exp_summary();
exp_summary_t(~ismember(exp_summary_t.batNum, bats),:)=[];
training_filename = "L:\Analysis\Code\inclusion_lists\training_history_per_bat.xlsx";
for ii_exp = 1:height(exp_summary_t)
    exp = exp_summary_t(ii_exp,:);
    T2 = readtable(training_filename,'Sheet',num2str(exp.batNum));
    surgery_date = T2.date(ismember(T2.comments,'surgery'));
    exp_summary_t{ii_exp,'session_num_from_exposure'} = sum(T2.long_flight_session == 1 & T2.date <= exp.date);
    exp_summary_t{ii_exp,'session_num_from_surgery'} = sum(T2.long_flight_session == 1 & T2.date <= exp.date & T2.date > surgery_date);
end
T = innerjoin(T,exp_summary_t,'Keys','exp_ID');

%% arrange to early and late sessions
% early_days_IX = find(T.session_num_from_exposure <= novelty_session_num_thr);
% late_days_IX = find(T.session_num_from_exposure > novelty_session_num_thr);
% early_late_IX = {early_days_IX,late_days_IX};
% exp_list = T.exp_ID;

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
if ~exist('directionality_contrast_index','var')
minSeqsThr = 20;
resampling_n = 20;
nrep = 1000;
directionality_contrast_index = [];
directionality_fraction = [];
directionality_binom_pval = [];
directionality_binom_surprise = [];
directionality_binom_surprise2 = nan(length(epoch_types)+1, length(exp_list));
directionality_binom_surprise2_sem = nan(length(epoch_types)+1, length(exp_list));
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
        n1 = sum([seqs.direction]==-1);
        n2 = sum([seqs.direction]==1);
        directionality_contrast_index(ii_epoch_type,ii_exp) = (n1-n2)/(n1+n2);
        directionality_fraction(ii_epoch_type,ii_exp) = (n1)/(n1+n2) .* sign(n1-n2);
        binom_pval = myBinomTest(n1,n1+n2,0.5);
        binom_pval(binom_pval==0) = realmin;
        directionality_binom_pval(ii_epoch_type,ii_exp) = binom_pval;
        directionality_binom_surprise(ii_epoch_type,ii_exp) = -log10(binom_pval) .* sign(n1-n2);
        nSeqs(ii_epoch_type,ii_exp) = length(seqs);
        if length(seqs)>=resampling_n
            sdf=zeros(1,nrep);
            for ii=1:nrep
                IX = randperm(length(seqs),resampling_n);
                n1 = sum([seqs(IX).direction]==-1);
                n2 = sum([seqs(IX).direction]==1);
                sdf(ii) = sign(n1-n2) * -log10(myBinomTest(n1,n1+n2,0.5));
            end
            directionality_binom_surprise2(ii_epoch_type,ii_exp) = nanmedian(sdf);
        end
    end
end
directionality_vals = cat(3,directionality_contrast_index,directionality_fraction,directionality_binom_pval,directionality_binom_surprise);
directionality_labels = {'contrast index','fraction','binom pval','binom surprise'};
% bat_sym_map = containers.Map(num2cell([184;2382]),{'+','x'});
end

%% save replay directionality results
filename = fullfile(res_dir, 'replay_directionality');
save(filename, ...
    "directionality_contrast_index", ...
    "directionality_fraction", ...
    "directionality_binom_pval", ...
    "directionality_binom_surprise", ...
    "directionality_binom_surprise2", ...
    "directionality_binom_surprise2_sem", ...
    "nSeqs","epoch_types", ...
    "minSeqsThr","resampling_n","nrep", ...
    "bats","bats_clr_map","T");

%% replay directionality bias - example session
% replay_directionality_bias_ex_exp_ID_list = {'b0184_d191129','b0184_d191130'};
replay_directionality_bias_ex_exp_ID_list = {'b0184_d191130'};
for ii_ex = 1:length(replay_directionality_bias_ex_exp_ID_list)
    replay_directionality_bias_ex_exp_ID = replay_directionality_bias_ex_exp_ID_list{ii_ex};
    ii_exp = find(strcmp(T.exp_ID,replay_directionality_bias_ex_exp_ID));
    example_session_num_form_exposure = T.session_num_from_exposure(replay_directionality_bias_ex_exp_ID);
    for ii_epoch_type = 1:length(epoch_types)
        axes(panels{1}(ii_ex,ii_epoch_type))
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
        text(0,-0.09,'Door','Units','normalized','FontSize',7,'HorizontalAlignment','center')
        switch epoch_types{ii_epoch_type}
            case 'sleep'
                title_str = 'Sleep replays';
            case 'rest'
                title_str = 'Awake replays';
        end
        title(title_str, 'Units','normalized', 'Position',[0.5 0.99],'FontWeight','normal');
    end
    title_str = {
        sprintf('Example session #%d:', example_session_num_form_exposure); ...
        'over-representation of replay directionality'};
    text(-.2, 1.12, title_str,'FontSize',9,'HorizontalAlignment','center','Units','normalized')
    for ii_dir = 1:2
        arrow_str = {'\rightarrow','\leftarrow'};
        arrow_offset_y = [0.02 0];
        text(-0.4, 0.99+arrow_offset_y(ii_dir), arrow_str{ii_dir},'Color',directions_clrs{ii_dir},'FontWeight','bold', 'FontSize',10,'Units','normalized');
    end
    
end

%% replay directionality bias - plot trend over exposure to enviroenment
clrs = [epoch_type_clrs,'k'];
surprise_clipping = 50;
for ii_epoch_type = 3%1:length(epoch_types)
    axes(panels{2}(1));
    cla reset
    hold on
    yline(0,'-','Color',[1 1 1]*0.8,'LineWidth',0.2);
%     c = clrs{ii_epoch_type};
    x = T.session_num_from_exposure;
    y = directionality_binom_surprise(ii_epoch_type,:);
%     y(~ismember(T.batNum,novelty_exposure_bats)) = nan;
    y(nSeqs(ii_epoch_type,:)<minSeqsThr) = nan;
    y(y>surprise_clipping) = surprise_clipping;
   y(y<-surprise_clipping) = -surprise_clipping;
   ylimits = [-1 1]*surprise_clipping;
    for ii_bat = 1:length(bats)
        bat_num = bats(ii_bat);
        c = bats_clr_map(bat_num);
        IX = T.bat_num==bat_num;
        xx = x(IX);
        yy = y(IX);
        invalid = isnan(yy);
        xx(invalid) = [];
        yy(invalid) = [];
        plot(xx,yy,'o-','color',c,'DisplayName',"bat "+bat_num,'MarkerFaceColor',c,'MarkerSize',2);
    end

%     nPointsSmooth = 3;
%     k = (nPointsSmooth-1)/2;
%     xi = [-1 1].*k + [(1-k):(max(x)+k)]';
%     xx=[];
%     yy=[];
%     for ii_xi = 1:size(xi,1)
%         TF = x>=xi(ii_xi,1) & x<=xi(ii_xi,2);
%         xx(ii_xi) = nanmedian(x(TF));
%         yy(ii_xi) = nanmedian(y(TF));
%     end
%     plot(xx,yy,'-','Color',c);
    ylim(ylimits)
    ylabel('Replay directionality index')
%     text(0.5,0.8,epoch_types{ii_epoch_type},'Units','normalized')
end
hax=gca;
text(.5, 2.8, {'Replay directionality';'all sessions'},'FontSize',9,'HorizontalAlignment','center','Units','normalized')
% legend
% x = [0.3 2]+32;
% y = linspace(0.35,0,3)+0.6;
% plot(x,y(1)*[1 1],'Color',clrs{1},'LineWidth',1.5,'Clipping','off')
% plot(x,y(2)*[1 1],'Color',clrs{2},'LineWidth',1.5,'Clipping','off')
% plot(x,y(3)*[1 1],'Color',clrs{3},'LineWidth',1.5,'Clipping','off')
% x = x(end)+1;
% text(x,y(1), "Sleep", 'FontSize',7)
% text(x,y(2), "Awake", 'FontSize',7)
% text(x,y(3), "Combined", 'FontSize',7)
xlabel('Session no.','Units','normalized','Position',[0.5 -0.05]);
ylim([-1 1]*55)
yticks([-50 0 50])
xticks([1 40])
hax=gca;
hax.XRuler.TickLength(1) = 0.02;
hax.YRuler.TickLength(1) = 0.035;
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = 0;
% text(5.5,-0.75,"\leftarrow"+"session #"+example_session_num_form_exposure,'FontSize',8)
h=annotation('textarrow');
h.Parent=hax;
h.X = [6 5];
h.Y = [50 50];
% h.String = {'  SAS=4'};
h.FontSize = 7; h.HeadLength = 4; h.HeadWidth = 3; h.HeadStyle = 'cback2';
h.HorizontalAlignment = 'left';
h.VerticalAlignment = 'middle';
h=annotation('textarrow');
h.Parent=hax;
h.X = 16 + [0 0.1];
h.Y = -16 + [0 0.1];
% h.String = {'SAS=8  '};
h.FontSize = 7; h.HeadLength = 4; h.HeadWidth = 3; h.HeadStyle = 'cback2';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'top';
h=annotation('textarrow');
h.Parent=hax;
h.X = 25 + [0 1];
h.Y = -50 + [0 0];
% h.String = {'SAS=27  '};
h.FontSize = 7; h.HeadLength = 4; h.HeadWidth = 3; h.HeadStyle = 'cback2';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'middle';
h=annotation('textarrow');
h.Parent=hax;
h.X = 30 + [0 -0.1];
h.Y = -50 + [0 0];
% h.String = {'   SAS=19'};
h.FontSize = 7; h.HeadLength = 4; h.HeadWidth = 3; h.HeadStyle = 'cback2';
h.HorizontalAlignment = 'left';
h.VerticalAlignment = 'middle';

%% cell revision 
bats = unique(T.bat_num);
bats_colors = [1 0 0; 0 1 0; 0 0 1; 0 0 0; 1 0 1; 0 1 1; 0.8 0.8 0];
bats_clr_map = containers.Map(num2cell(bats),num2cell(bats_colors,2));
surprise_clipping = 50;
fig2=figure(WindowState="maximized");
tiledlayout(4,3);
for ii=1:4
    for ii_epoch_type = 1:length(epoch_types)+1
        nexttile
        hold on
        x = T.session_num_from_exposure;
       switch ii
           case {1,2}
               y = directionality_contrast_index(ii_epoch_type,:);
               ylabel_str = {'Replay directionality';'contrast index'};
               ylimits = [-1 1];
           case 3
               y = directionality_binom_surprise(ii_epoch_type,:);
               ylabel_str = {'Replay directionality';'surprise'};
               y(y>surprise_clipping) = surprise_clipping;
               y(y<-surprise_clipping) = -surprise_clipping;
               ylimits = [-1 1]*surprise_clipping;
           case 4
               y = directionality_binom_surprise2(ii_epoch_type,:);
               ylabel_str = {'Replay directionality';'surprise (resampling)'};
               ylimits = [-10 10];
       end
        y(nSeqs(ii_epoch_type,:)<minSeqsThr) = nan;
        yline(0)
        for ii_bat = 1:length(bats)
            bat_num = bats(ii_bat);
            c = bats_clr_map(bat_num);
            IX = find(T.bat_num==bat_num);
            xx = x(IX);
            yy = y(IX);
            invalid = isnan(yy);
            xx(invalid) = [];
            yy(invalid) = [];
            IX(invalid) = [];
            if ii==2
                sym='.-';
            else
                sym='o-';
            end
            hp = plot(xx,yy,sym,'color',c,'DisplayName',"bat "+bat_num,'MarkerFaceColor',c);
            if ~isempty(hp)
                hp.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('bat num',repelem(bat_num,length(xx)));
                hp.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('exp ID',T.exp_ID(IX));
                hp.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('#replays',nSeqs(ii_epoch_type,IX));
                hp.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('#fromExposure',T.session_num_from_exposure(IX));
                hp.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('#fromSurgery',T.session_num_from_surgery(IX));
            end
        end
        if ii==2
            sz = nSeqs(ii_epoch_type,:);
            sz(isnan(y)) = nan;
%             sz = interp1([0 50],[0 1],sz,'linear','extrap');
            c = arrayfun(@(bat)bats_clr_map(bat),T.bat_num,'UniformOutput',0);
            c = cat(1,c{:});
            hbb=bubblechart(x,y,sz,c);
            bubblelim([min(sz) max(sz)])
            bubblesize([3 20])
            bubblelegend('#replays','location','southeast')
        end
%         legend('NumColumns',length(bats),'Location',[0.1 0.9 0.8 0.1])
        legend('NumColumns',2)
%         nPointsSmooth = 3;
%         k = (nPointsSmooth-1)/2;
%         xi = [-1 1].*k + [(1-k):(max(x)+k)]';
%         xx=[];
%         yy=[];
%         for ii_xi = 1:size(xi,1)
%             TF = x>=xi(ii_xi,1) & x<=xi(ii_xi,2);
%             xx(ii_xi) = nanmedian(x(TF));
%             yy(ii_xi) = nanmedian(y(TF));
%         end
%         plot(xx,yy,'-','Color','k','LineWidth',2);
        titles_strs = [epoch_types, 'pooled'];
        title(titles_strs{ii_epoch_type});
        ylabel(ylabel_str)
        xlabel('Session no.')
        ylim(ylimits)
        xlim([0 max(x)+1])
    end
end
fig2_filename = 'E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\cell_revision\replay_directionality_population_options';
saveas(fig2,fig2_filename,'fig');
exportgraphics(fig2,[fig2_filename '.pdf'],'BackgroundColor','white');
% saveas(fig2,fig2_filename,'pdf');


%% ========================================================================




%% panels C - experimental setup
axes(panels{3}(1))
hold on
tunnel_2bats_image = imread('E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\resources\tunnel_2bats.jpg');
imshow(tunnel_2bats_image);

axes(panels{3}(2))
cla reset
hold on
axis off
plot([0 1],[1 1],'LineWidth',2)
plot([0 1],[0 0],'LineWidth',1)
text(1.5,1,'Recorded bat','Units','normalized','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
text(1.5,0,'Other bat','Units','normalized','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
xlim([0 1])
ylim([0 1])

%% panels D - Behavior example
axes(panels{4}(1))
cla
hold on
switch behavior_ex_opt
    case 1
        exp_ID = 'b2299_d191205';
        ti_seconds_in_session = [5200 5600];
    case 2
        exp_ID = 'b2299_d191213';
        ti_seconds_in_session = 5850 +[0 4*60]+[0.4 -0.1]+[0.35 -0.13].*60; %[5850 6200];
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
axes(panels{4}(2))
cla
hold on
plot([0 0.1],[0.8 0.8],'-','LineWidth',lw);
plot([0 0.1],[0.0 0.0],'-','LineWidth',1);
plot(0.6,0.8,'xk','MarkerSize',8);
plot([0.57 0.63],[0 0],'-k','LineWidth',1.);
plot(0.57,0,'.k','MarkerSize',7);
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

% TT = exp.ripples.stats.best_TT;
% [LFP.signal, LFP.ts, LFP.fs, LFP.params] = LFP_load(exp_ID,TT,'band','ripple','limits_ts',ti);
% LFP.avg_signal = nanmean(LFP.signal,[2 3]);

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
axes(panels{5}(2));
cla
hold on
IX = get_data_in_ti(decode.time, ti);
prob_t = decode.time(IX);
prob_state = squeeze(decode.posterior_state(event.state_num,IX));
plot(prob_t, prob_state, 'k','LineWidth',2);
hax = gca;
hax.XLim = prob_t([1 end]);
hax.YLim = [0 1];
hax.XTickLabel = [];
box on
hax.Colormap = cmap;
hax.XLim = seq_ti+[-1 1].*0.2*range(seq_ti);
hax.TickDir = 'out';
hax.TickLength = [0.02 0.02];
hax.XRuler.TickLabelGapOffset = -4;
ylabel('Prob.','Units','normalized','Position',[-0.22 0.5]);
rescale_plot_data('x',[1e-6 seq_ti(1)]);
    
%% plot posterior (position)
axes(panels{5}(1));
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
hax.Colormap = cmap;
hax.XLim = seq_ti+[-1 1].*0.2*range(seq_ti);
hax.TickDir = 'out';
hax.TickLength = [0.02 0.02];
hax.XRuler.TickLabelGapOffset = -2;
xlabel('Time (s)','Units','normalized','Position',[0.5 -0.105]);
ylabel('Position (m)','Units','normalized','Position',[-0.22 0.5]);
rescale_plot_data('x',[1e-6 seq_ti(1)]);
plot(hax.XLim([1 1]),hax.YLim,'k-')
plot(hax.XLim([2 2]),hax.YLim,'k-')

%% link x axes
linkaxes(panels{5},'x');

%% add colorbar
hcb = colorbar('east');
hcb.Units = 'centimeters';
cb_offset = 1.12;
hcb.Position(1) = panels{5}(1).Position(1) + panels{5}(1).Position(3)*cb_offset;
cb_offset_middle = (hcb.Position(1)+hcb.Position(3)/2-panels{5}(1).Position(1))/panels{5}(1).Position(3);
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
axes(panels{6}(1,1))
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
% text(0.8,0.3,"n = "+ data.stats.n,'Units','normalized','FontSize',9);
% text(0.6,0.2,"r = "+ sprintf('%.2g',data.stats.Pearson.r), 'Units','normalized','FontSize',7);
% text(0.6,0.1,"P = "+ sprintf('%.2g',data.stats.Pearson.p), 'Units','normalized','FontSize',7);
text(.05,1.1,"{\rho} = "+ sprintf('%.2f',data.stats.Spearman.r), 'Units','normalized','FontSize',7);
text(.05,1.0,"P = "+ sprintf('%.2g',data.stats.Spearman.p), 'Units','normalized','FontSize',7);
% text(.3,1.3,"{\itn} = "+ sprintf('%d',sum(data.TF)), 'Units','normalized','FontSize',7);
% text(0,-.4,data.msg_str, 'Units','normalized','FontSize',10);
h=refline(1,0);
h.Color = .8.*[1 1 1];
hax=gca;
hax.XRuler.TickLength(1) = 0.035;
hax.YRuler.TickLength(1) = 0.024;
hax.XRuler.TickLabelGapOffset = -.5;
hax.YRuler.TickLabelGapOffset = 1;

%% shuffling test (revision)
nreps = 1000;
rng(0);
r_shuffles = zeros(1,nreps);
for ii = 1:nreps
    IX = randperm(length(X));
    r_shuffles(ii) = corr(X(IX)',Y','type','Spearman');
end
r = corr(X',Y','type','Spearman');
z = (r-mean(r_shuffles))./std(r_shuffles);
pval_nonparam = mean(r<r_shuffles);
[~,pval_ttest] = ttest(r_shuffles,r,'Tail','left');
pval_z = (1-normcdf(z));
[~,pval_z2] = ztest(r_shuffles,r,std(r_shuffles),'Tail','left');
% fig3 = figure(Units="centimeters",Position=[5 5 20 20]);
% hold on
%%
axes(panels{6}(1,2))
cla reset
hold on
histogram(r_shuffles,'FaceColor',[1 1 1]*0.5);
xline(r,'-r');
xlabel('{\rho}','Units','normalized','Position',[0.5 .1])
% title({'2-bats crossover replay correlation, shuffle analysis';sprintf('pval non-parametric (vs shuffles) = %.2g (n=%d shuffles)',pval_nonparam,nreps)});
% fig3_filename = 'E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\cell_revision\2bats_shuffle_corr';
% exportgraphics(fig3,[fig3_filename '.pdf'],'BackgroundColor','white');
% saveas(fig3,fig3_filename,'pdf');
ylim([0 120])
xticks([])
yticks([])
text(0.5,1.25,"P = "+ sprintf('%.2g',pval_nonparam), 'Units','normalized','FontSize',7,'HorizontalAlignment','Center');

%% 2-back previous CO (revision)
% fig4 = figure(Units="centimeters",Position=[5 5 6 6]);
% hold on
% data = load(data_filename);
% % X = [data.seqs_all.prev_co_pos];
% X = [data.seqs_all.prev_2_co_pos];
% Y = data.y;
% X = X(data.TF)';
% Y = Y(data.TF)';
% [r,pval] = corr(X,Y,'type','Spearman',rows='pairwise',tail='right');
% scatter(X,Y,5,'k','filled');
% text(.3,1.05,"{\rho} = "+ sprintf('%.2f',r), 'Units','normalized','FontSize',7);
% text(.3,0.95,"P = "+ sprintf('%.2g',pval), 'Units','normalized','FontSize',7);
% h=refline(1,0);
% h.Color = .8.*[1 1 1];
% axis equal
% xlim([0 135])
% ylim([0 135])
% xticks(linspace(0,135,4))
% yticks(linspace(0,135,4))
% xlabel('2-back previous cross-over position (m)');
% ylabel('Replay position (m)');
% fig4_filename = 'E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\cell_revision\2bats_2back_prev_CO';
% saveas(fig4,fig4_filename,'pdf');

%% panels D - scatter plot (control - next crossover)
axes(panels{6}(2,1))
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
% text(0.8,0.3,"n = "+ data.stats.n,'Units','normalized','FontSize',9);
% text(0.8,0.2,"r = "+ sprintf('%.2g',stats.Pearson.r), 'Units','normalized','FontSize',7);
% text(0.8,0.1,"P = "+ sprintf('%.2g',stats.Pearson.p), 'Units','normalized','FontSize',7);
text(0.05,1.1,"{\rho} = "+ sprintf('%.2f',stats.Spearman.r), 'Units','normalized','FontSize',7);
text(0.05,1.0,"P = "+ sprintf('%.2f',stats.Spearman.p), 'Units','normalized','FontSize',7);
h=refline(1,0);
h.Color = .8.*[1 1 1];
hax=gca;
hax.XRuler.TickLength(1) = 0.035;
hax.YRuler.TickLength(1) = 0.024;
hax.XRuler.TickLabelGapOffset = -.5;
hax.YRuler.TickLabelGapOffset = 1;

%% shuffling test (revision)
nreps = 1000;
rng(0);
r_shuffles = zeros(1,nreps);
for ii = 1:nreps
    IX = randperm(length(X));
    r_shuffles(ii) = corr(X(IX),Y,'type','Spearman');
end
r = corr(X,Y,'type','Spearman');
z = (r-mean(r_shuffles))./std(r_shuffles);
pval_nonparam = mean(r<r_shuffles);
[~,pval_ttest] = ttest(r_shuffles,r,'Tail','left');
pval_z = (1-normcdf(z));
[~,pval_z2] = ztest(r_shuffles,r,std(r_shuffles),'Tail','left');
% fig3 = figure(Units="centimeters",Position=[5 5 20 20]);
% hold on
%%
axes(panels{6}(2,2))
cla reset
hold on
histogram(r_shuffles,'FaceColor',[1 1 1]*0.5);
xline(r,'-r');
xlabel('{\rho}','Units','normalized','Position',[0.5 .1])
% title({'2-bats crossover replay correlation, shuffle analysis';sprintf('pval non-parametric (vs shuffles) = %.2g (n=%d shuffles)',pval_nonparam,nreps)});
% fig3_filename = 'E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\cell_revision\2bats_shuffle_corr';
% exportgraphics(fig3,[fig3_filename '.pdf'],'BackgroundColor','white');
% saveas(fig3,fig3_filename,'pdf');
ylim([0 120])
xticks([])
yticks([])
text(0.5,1.25,"P = "+ sprintf('%.2g',pval_nonparam), 'Units','normalized','FontSize',7,'HorizontalAlignment','Center');

%% add panel letters
font_size = 11;
axes(panels{1}(1))
text(-0.3,1.1, 'A', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{2}(1))
text(-0.3,1.25, 'B', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{3}(1))
text(-0.17,1.2, 'C', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{3}(2))
axes(panels{4}(1))
text(-0.2,1.3, 'D', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{5}(1))
text(-0.38,1.3, 'E', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{6}(1,1))
text(-0.35,1.15, 'F', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{6}(2,1))
text(-0.35,1.15, 'G', 'Units','normalized','FontWeight','bold','FontSize',font_size);

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







%%
