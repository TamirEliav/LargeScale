%% Replay - Fig supp 6 - Novelty (days 1-5 vs 6+)
clear 
clc
close all

%% data options 
params_opt = 11; % decoding opt 
novelty_session_num_thr = 5;
minSeqsThr = 30;

%% plotting options
panels_xlim = [-0.0950    1.9950; -1.9500   40.9500; -0.4500   56.5500; -0.0210    0.4410];
panels_xticks = {[0 0.5 1 1.5],[0 20 40],[0:10:50],[0 0.2 0.4]};
panels_bin_size = [0.05 1 1 0.01];
boxplot_prctiles = [10 25 50 75 90];

%% graphics params
epoch_types = {'sleep','rest'};
days_types = {'Days 1-5','Days 6+'};
epoch_type_clrs = {[.6 .1 .8],[.1 .8 .1]};
directions_clrs = {[0    0.3843    0.7451];[ 0.5216    0.2471         0]};
% directions_clrs = [0    0.3843    0.7451; 0.5216    0.2471         0];
line_styles = {'-','--'};
epoch_types_str = {'Sleep','Awake'};

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Figure S6';
fig_caption_str = 'Novelty effects';
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
w = 3.5;
h = 2;
x = [3 7.5 12 16.5];
y = 21.5;
panels{1}(1,1) = axes('position', [x(1) y w h]);
panels{1}(1,2) = axes('position', [x(2) y w h]);
panels{1}(1,3) = axes('position', [x(3) y w h]);
panels{1}(1,4) = axes('position', [x(4) y w h]);
y= 18.5;
panels{1}(2,1) = axes('position', [x(1) y w h]);
panels{1}(2,2) = axes('position', [x(2) y w h]);
panels{1}(2,3) = axes('position', [x(3) y w h]);
panels{1}(2,4) = axes('position', [x(4) y w h]);

% w = 3;
% h = 2.5;
% y = 14;
% panels{2}(1) = axes('position', [x(1) y w h]);
% panels{2}(2) = axes('position', [x(2) y w h]);
% panels{2}(3) = axes('position', [x(3) y w h]);
% panels{2}(4) = axes('position', [x(4) y w h]);

h = 2.5;
w = 2.9;
y = 13.7;
panels{3}(1,1) = axes('position', [x(1) y w h]);
panels{3}(1,2) = axes('position', [x(2) y w h]);
panels{3}(1,3) = axes('position', [x(3) y w h]);
panels{3}(1,4) = axes('position', [x(4) y w h]);
y = 10;
panels{3}(2,1) = axes('position', [x(1) y w h]);
panels{3}(2,2) = axes('position', [x(2) y w h]);
panels{3}(2,3) = axes('position', [x(3) y w h]);
panels{3}(2,4) = axes('position', [x(4) y w h]);

panels{3}(1,2).Position(1) = panels{3}(1,2).Position(1)+0.2;
panels{3}(2,2).Position(1) = panels{3}(2,2).Position(1)+0.2;
panels{3}(1,4).Position(1) = panels{3}(1,4).Position(1)+0.2;
panels{3}(2,4).Position(1) = panels{3}(2,4).Position(1)+0.2;

total_offset = [0 0];
for ii = 1:length(panels)
    subpanels = panels{ii};
    subpanels = subpanels(:);
    for jj = 1:length(subpanels)
        subpanels(jj).Position([1 2]) = subpanels(jj).Position([1 2]) + total_offset;
    end
end

%% properties to plot
features_names = {'duration';'compression';'distance';'distance_norm';};
fn_label_strs = {
    'Replay duration (s)';
    {'Compression ratio';'(replay speed / flight speed)'};
    'Replay distance (m)';
    {'Replay distance','(norm. to environment size)'};
    };

%% arrange sessions to load
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

%% plot histograms
for ii_EL = 1:length(early_late_IX)
    ii_EL
    for ii_fn = 1:length(features_names)
        fn = features_names{ii_fn};
        axes(panels{1}(ii_EL,ii_fn));
        cla
        hold on
        lw = 1.1;
        clear hh;
        for ii_epoch_type = 1:length(epoch_types)
            events = [events_all_per_session{ii_epoch_type,early_late_IX{ii_EL}}];
            seqs = [events.seq_model];
            X = [seqs.(fn)];
            hh(ii_epoch_type) = histogram(X,'DisplayStyle','stairs','Normalization','pdf','EdgeColor',epoch_type_clrs{ii_epoch_type},'LineWidth',lw,'BinWidth',panels_bin_size(ii_fn));
            text(0.5,0.9-ii_epoch_type*0.14, "n_{" + epoch_types_str{ii_epoch_type} + "} = "+length(X),'units','normalized','FontSize',7)
        end
        linkprop(hh,'BinEdges');
        xlim(panels_xlim(ii_fn,:))
        xticks(panels_xticks{ii_fn})
        if ii_EL == 2
            xlabel(fn_label_strs{ii_fn});
        end
        if ii_fn == 1
            text(-0.5,0.5,days_types{ii_EL},'Units','normalized','HorizontalAlignment','center','FontSize',8);
            ylabel('Probability density','Units','normalized','Position',[-0.13 0.5]);
        end
        hax=gca;
    %     hax.TickLength(1) = [0.025];
        hax.XRuler.TickLength(1) = 0.03;
        hax.XRuler.TickLabelGapOffset = -1.2;
        hax.YRuler.TickLabelGapOffset = 1;
        fprintf('%s: median %.2g, mean = %.2g\n',fn,median(X),mean(X))
    end
end

%% add sleep/awake legend
if exist('panels_hist_legend','var')
    delete(panels_hist_legend);
end
panels_hist_legend = axes('position', [3.8 23.5 0.3 0.25]);
cla
hold on
plot([0 1], [1 1], 'color', epoch_type_clrs{1}, 'LineWidth',lw,'Clipping','off');
plot([0 1], [0 0], 'color', epoch_type_clrs{2}, 'LineWidth',lw,'Clipping','off');
text(1.3, 1, 'Sleep','FontSize',7,'HorizontalAlignment','left');
text(1.3, 0, 'Awake','FontSize',7,'HorizontalAlignment','left');
xlim([0 1])
ylim([0 1])
axis off

%% plot boxplots
boxplot_panels_ylimits = [0 0.8; 0 25; 3 23; .02 .18];
% pvals = [];
% for ii_fn = 1:length(features_names)
%     axes(panels{2}(ii_fn));
%     cla
%     hold on
%     fn = features_names{ii_fn};
%     X={};
%     G={};
%     for ii_epoch_type = 1:length(epoch_types)
%         for ii_EL = 1:length(early_late_IX)
%             IX = early_late_IX{ii_EL};
%             events = [events_all_per_session{ii_epoch_type,IX}];
%             seqs = [events.seq_model];
%             x = [seqs.(fn)];
%             g = (ii_EL+2*(ii_epoch_type-1));
%             X{ii_epoch_type,ii_EL} = x;
%             G{ii_epoch_type,ii_EL} = ones(size(seqs)).*g;
%             m1 = prctile(x,50);
%             m2 = prctile(x,[25 75]);
%             m3 = prctile(x,[5 95]);
%             w = 0.2;
%             lw = 1.3;
%             plot([g-w g+w],[m1 m1],'Color',epoch_type_clrs{ii_epoch_type},'LineWidth',lw);
%             plot([g g],m3,'Color',epoch_type_clrs{ii_epoch_type},'LineWidth',lw);
%             rectangle('Position',[g-w m2(1) 2*w diff(m2)],'EdgeColor',epoch_type_clrs{ii_epoch_type},'FaceColor','none','LineWidth',lw);
%         end
%         pval = ranksum(X{ii_epoch_type,1},X{ii_epoch_type,2});
%         str = genSignifStrAstricks(pval);
%         xx = [g g-1];
%         yy = boxplot_panels_ylimits(ii_fn,[2 2]);
%         plot(xx,yy,'k-');
%         font_size = 10;
%         if strcmp(str,'n.s.')
%             font_size = 8;
%             yy = yy.*1.02;
%         end
%         text(mean(xx),mean(yy),str,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',font_size);
%         pvals(ii_fn,ii_epoch_type) = pval;
%     end
%     xlim([0.2 4.8])
%     ylim(boxplot_panels_ylimits(ii_fn,:))
%     xlabel('');
%     ylabel_pos_x = [-.23 -.21 -.21 -.25];
%     ylabel(fn_label_strs{ii_fn},'units','normalized','position',[ylabel_pos_x(ii_fn) 0.5]);
%     xticks(1:4);
%     xticklabels(["Sleep, days 1-5", "Sleep, days 6+","Awake, days 1-5", "Awake, days 6+"]);
%     xtickangle(50);    
% end


%% boxplots per day
maxDays2plot = 10;
for ii_fn = 1:length(features_names)
    for ii_epoch_type = 1:length(epoch_types)
        axes(panels{3}(ii_epoch_type,ii_fn));
        cla reset
        hold on
        fn = features_names{ii_fn};
        events = [events_all_per_session{ii_epoch_type,:}];
        session_num = [events.session_num_from_exposure];
        g = session_num;
        g(g>maxDays2plot) = maxDays2plot;
        [G,ID] = findgroups(g);
        seqs = [events.seq_model];
        x = [seqs.(fn)];
        m = 0;
        for ii = 1:length(ID)
            id = ID(ii);
            IX = G==id;
            xg = x(IX);
            xg_prctl = prctile(xg,boxplot_prctiles);
            c = epoch_type_clrs{ii_epoch_type};
            w = 0.2;
            xx = ii+w*[-1 1];
            plot(xx, xg_prctl([3 3]),'Color',c);
            plot(xx, xg_prctl([2 2]),'Color',c);
            plot(xx, xg_prctl([4 4]),'Color',c);
            plot(xx([1 1]), xg_prctl([2 4]),'Color',c);
            plot(xx([2 2]), xg_prctl([2 4]),'Color',c);
            plot([ii ii], xg_prctl([1 5]),'Color',c);
            m = max(m,xg_prctl(5));
        end
        ylimits = [0 m*1.05];
        ylim(ylimits)
%         boxplot(x,g,'PlotStyle','traditional','BoxStyle','outline','Colors',epoch_type_clrs{ii_epoch_type},'Symbol','','Notch','off');
        box off
%         ylim(boxplot_panels_ylimits(ii_fn,:))
        hax=gca;
        hax.XTickLabelRotation = 0;
        hax.XTick = 1:length(ID);
        hax.XTickLabel{end} = ['   ' char(8805) num2str(maxDays2plot)];
%         hax.XTickLabel{end} = [num2str(maxDays2plot) '+'];
        hax.XRuler.FontSize = 7;
        hax.XRuler.TickLength(1) = 0.03;
        hax.XRuler.TickLabelGapOffset = -1.;
%         hax.YRuler.TickLabelGapOffset = 1;
        if ii_epoch_type == 1
            ylabel_xpos = [-0.21 -0.21 -0.21 -0.25];
            ylabel(fn_label_strs{ii_fn},'units','normalized','position',[ylabel_xpos(ii_fn) -0.3])
        end
        if ii_epoch_type == 2
            xlabel('Session no.','units','normalized','position',[0.5 -0.2])
            if ii_fn==1
                text(.5,-0.4,'(Day no.)','Units','normalized','FontSize',7.7,'HorizontalAlignment','center');
            end
        end
        if ii_fn==1
            legend_strs = {'Sleep';'Awake'};
            str = legend_strs{ii_epoch_type};
            x = hax.XLim(1)+[0.6 0.75].*range(hax.XLim);
            y = hax.YLim(1)+[1.05 1.05].*range(hax.YLim);
            plot(x,y,'Color',epoch_type_clrs{ii_epoch_type},'LineWidth',1.5,'Clipping','off');
%             plot(x(2),y(2),'.','Color',epoch_type_clrs{ii_epoch_type},'LineWidth',1.5)
            text(x(2)+0.2*diff(x),y(1),str,'FontSize',7)
        end
    end
end



%% add panel letters
font_size = 11;
axes(panels{1}(1,1))
text(-0.3,1.2, 'A', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{3}(1))
text(-0.4,1.1, 'B', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
% fig_name = sprintf('%s_decoding_opt_%d',fig_name_str, params_opt);
fig_name = sprintf('%s_boxplot_prctiles_%d_%d_%d_%d_%d',fig_name_str, boxplot_prctiles);

file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')






%% playground...
exp_ID = 'b0184_d191202';
exp = exp_load_data(exp_ID,'pos','rest');
figure
tiledlayout("flow")
nexttile
plot(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos)
nexttile
hold on
rest_durations = [exp.rest.events.duration];
rest_ball_num = [exp.rest.events.ball_num];
N = splitapply(@sum,rest_durations,rest_ball_num)
ecdf(rest_durations(rest_ball_num==1))
ecdf(rest_durations(rest_ball_num==2))
nexttile
hold on
histogram(rest_durations(rest_ball_num==1),'Normalization','count','DisplayStyle','stairs')
histogram(rest_durations(rest_ball_num==2),'Normalization','count','DisplayStyle','stairs')
[~,ks_pval] = kstest2(rest_durations(rest_ball_num==1),rest_durations(rest_ball_num==2))

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