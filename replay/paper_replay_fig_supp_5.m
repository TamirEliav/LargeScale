%% Replay - Fig supp 5 - Novelty (days 1-5 vs 6+)
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
fig_name_str = 'Extended_Data_Fig_5';
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

w = 3;
h = 2.5;
y = 14;
panels{2}(1) = axes('position', [x(1) y w h]);
panels{2}(2) = axes('position', [x(2) y w h]);
panels{2}(3) = axes('position', [x(3) y w h]);
panels{2}(4) = axes('position', [x(4) y w h]);

h = 2.5;
w = 3;
y = 10;
panels{3}(1,1) = axes('position', [x(1) y w h]);
panels{3}(1,2) = axes('position', [x(2) y w h]);
panels{3}(1,3) = axes('position', [x(3) y w h]);
panels{3}(1,4) = axes('position', [x(4) y w h]);
y = 7;
panels{3}(2,1) = axes('position', [x(1) y w h]);
panels{3}(2,2) = axes('position', [x(2) y w h]);
panels{3}(2,3) = axes('position', [x(3) y w h]);
panels{3}(2,4) = axes('position', [x(4) y w h]);

w = 3;
w2 = 2;
h = 2.5;
y = 3.5;
panels{4}(1) = axes('position', [x(1) y w h]);
panels{4}(2) = axes('position', [x(1)+w+0.5 y w2 h]);
% panels{4}(1) = axes('position', [x(1) y w h]);
% panels{4}(2) = axes('position', [x(2) y w h]);
% panels{4}(3) = axes('position', [x(3) y w h]);
% panels{4}(4) = axes('position', [x(4) y w h]);

w = 8;
h = 3
y = -1;
panels{5}(1) = axes('position', [x(1) y w h]);
panels{5}(2) = axes('position', [x(3) y w h]);

total_offset = [0 3];
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
            text(0.5,0.9-ii_epoch_type*0.1, "{\itn}_{" + epoch_types_str{ii_epoch_type} + "} = "+length(X),'units','normalized','FontSize',7)
        end
        linkprop(hh,'BinEdges');
        xlim(panels_xlim(ii_fn,:))
        xticks(panels_xticks{ii_fn})
        if ii_EL == 2
            xlabel(fn_label_strs{ii_fn});
        end
        if ii_fn == 1
            text(-0.5,0.5,days_types{ii_EL},'Units','normalized','HorizontalAlignment','center');
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

%% plot boxplots
boxplot_panels_ylimits = [0 0.8; 0 25; 3 23; .02 .18];
pvals = [];
for ii_fn = 1:length(features_names)
    axes(panels{2}(ii_fn));
    cla
    hold on
    fn = features_names{ii_fn};
    X={};
    G={};
    for ii_epoch_type = 1:length(epoch_types)
        for ii_EL = 1:length(early_late_IX)
            IX = early_late_IX{ii_EL};
            events = [events_all_per_session{ii_epoch_type,IX}];
            seqs = [events.seq_model];
            x = [seqs.(fn)];
            g = (ii_EL+2*(ii_epoch_type-1));
            X{ii_epoch_type,ii_EL} = x;
            G{ii_epoch_type,ii_EL} = ones(size(seqs)).*g;
            m1 = prctile(x,50);
            m2 = prctile(x,[25 75]);
            m3 = prctile(x,[5 95]);
            w = 0.2;
            lw = 1.3;
            plot([g-w g+w],[m1 m1],'Color',epoch_type_clrs{ii_epoch_type},'LineWidth',lw);
            plot([g g],m3,'Color',epoch_type_clrs{ii_epoch_type},'LineWidth',lw);
            rectangle('Position',[g-w m2(1) 2*w diff(m2)],'EdgeColor',epoch_type_clrs{ii_epoch_type},'FaceColor','none','LineWidth',lw);
        end
        pval = ranksum(X{ii_epoch_type,1},X{ii_epoch_type,2});
        str = genSignifStrAstricks(pval);
        xx = [g g-1];
        yy = boxplot_panels_ylimits(ii_fn,[2 2]);
        plot(xx,yy,'k-');
        font_size = 10;
        if strcmp(str,'n.s.')
            font_size = 8;
            yy = yy.*1.02;
        end
        text(mean(xx),mean(yy),str,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',font_size);
        pvals(ii_fn,ii_epoch_type) = pval;
    end
    xlim([0.2 4.8])
    ylim(boxplot_panels_ylimits(ii_fn,:))
    xlabel('');
    ylabel_pos_x = [-.23 -.21 -.21 -.25];
    ylabel(fn_label_strs{ii_fn},'units','normalized','position',[ylabel_pos_x(ii_fn) 0.5]);
    xticks(1:4);
    xticklabels(["Sleep, days 1-5", "Sleep, days 6+","Awake, days 1-5", "Awake, days 6+"]);
    xtickangle(50);    
end

%% add sleep/awake legend
if exist('panels_hist_legend','var')
    delete(panels_hist_legend);
end
panels_hist_legend = axes('position', [3.8 26.5 0.3 0.25]);
cla
hold on
plot([0 1], [1 1], 'color', epoch_type_clrs{1}, 'LineWidth',lw,'Clipping','off');
plot([0 1], [0 0], 'color', epoch_type_clrs{2}, 'LineWidth',lw,'Clipping','off');
text(1.3, 1, 'Sleep','FontSize',7,'HorizontalAlignment','left');
text(1.3, 0, 'Awake','FontSize',7,'HorizontalAlignment','left');
xlim([0 1])
ylim([0 1])
axis off

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
        seqs = [events.seq_model];
        x = [seqs.(fn)];
        boxplot(x,g,'PlotStyle','traditional','BoxStyle','outline','Colors',epoch_type_clrs{ii_epoch_type},'Symbol','','Notch','off');
        box off
        ylim(boxplot_panels_ylimits(ii_fn,:))
        hax=gca;
        hax.XTickLabelRotation = 0;
        hax.XTickLabel{end} = [char(8805) num2str(maxDays2plot)];
        ylabel(fn_label_strs{ii_fn})
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

%% replay directionality bias - plot trend over exposure to enviroenment
clrs = [epoch_type_clrs,'k'];
axes(panels{4}(1));
cla reset
hold on
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

%% replay directionality bias - plot non-novelty hist
axes(panels{4}(2));
cla reset
hold on
for ii_epoch_type = 1:length(epoch_types)+1
    c = clrs{ii_epoch_type};
    y = directionality_contrast_index(ii_epoch_type,:);
    TF = true(size(y));
    TF = TF & ~(T.session_num_from_exposure <= novelty_session_num_thr)';
%     TF = TF & ~ismember(T.bat_num,novelty_exposure_bats)';
    TF = TF & nSeqs(ii_epoch_type,:)>=minSeqsThr;
    y = y(TF);
    histogram(y,'Normalization','count','DisplayStyle','stairs','EdgeColor',c,'BinWidth',0.1);
end
xlim([-1 1])
view(90,-90)

%% replay directionality bias - example session
replay_directionality_bias_ex_exp_ID = 'b0184_d191130';
ii_exp = find(strcmp(T.exp_ID,replay_directionality_bias_ex_exp_ID))
for ii_epoch_type = 1:length(epoch_types)
    axes(panels{5}(ii_epoch_type))
    cla reset
    hold on
    events = events_all_per_session{ii_epoch_type,ii_exp};
    seqs = [events.seq_model];
    epoch_sep = find(diff([events.epoch_num])~=0)+0.5;
    seqs_edges = [seqs.start_pos_norm; seqs.end_pos_norm];
    seqs_IX = 1:length(seqs);
    if strcmp(epoch_types(ii_epoch_type),'sleep')
        plot(repmat(epoch_sep,2,1),[0 1],':','LineWidth',1.5,'Color',0.5*[1 1 1]);
    end
    h=plot([seqs_IX;seqs_IX],seqs_edges,'-','LineWidth',.55);
    dir_1_IX = [seqs.state_direction]==1;
    dir_2_IX = [seqs.state_direction]==-1;
    [h(dir_1_IX).Color] = disperse(repelem(directions_clrs(1),length(dir_1_IX)));
    [h(dir_2_IX).Color] = disperse(repelem(directions_clrs(2),length(dir_2_IX)));
    scatter(seqs_IX,seqs_edges(1,:), 3, [seqs.state_direction]==-1, "filled");
    hax=gca;
    hax.Colormap = cell2mat(directions_clrs);
    hax.YTick = [0 1];
    hax.XTick = [1 10*ceil(length(seqs)/10)];
    hax.YLim = [0 1];
    hax.XRuler.TickLabelGapOffset = -1;
    hax.YRuler.TickLabelGapOffset = 1;
    ylabel('Position (norm.)', 'Units','normalized', 'Position',[-0.017 0.5]);
    xlabel('Replay event no.', 'Units','normalized', 'Position',[0.5 -0.06]);
%     title('Sleep replays', 'Units','normalized', 'Position',[0.47 0.99],'FontWeight','normal');
%     text(1.1, 1.05, 'Example session: all individual replays','FontSize',9,'HorizontalAlignment','center','Units','normalized')
end

%% add panel letters
font_size = 11;
axes(panels{1}(1,1))
text(-0.3,1.1, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{2}(1))
text(-0.4,1.1, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = sprintf('%s_decoding_opt_%d',fig_name_str, params_opt);
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