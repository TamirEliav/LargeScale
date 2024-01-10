%% Replay - Fig supp 5 - Novelty (days 1-5 vs 6+)
clear 
clc
close all

%% data options 
params_opt = 11; % decoding opt 

%% plotting options
panels_xlim = [-0.0950    1.9950; -1.9500   40.9500; -0.4500   56.5500; -0.0210    0.4410];
panels_xticks = {[0 0.5 1 1.5],[0 20 40],[0:10:50],[0 0.2 0.4]};
panels_bin_size = [0.05 1 1 0.01];

%% graphics params
epoch_types = {'sleep','rest'};
days_types = {'Days 1-5','Days 6+'};
epoch_type_clrs = {[.6 .1 .8],[.1 .8 .1]};
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
panels{1}(1,1) = axes('position', [3   20 3.5 3]);
panels{1}(1,2) = axes('position', [7.5 20 3.5 3]);
panels{1}(1,3) = axes('position', [12 20 3.5 3]);
panels{1}(1,4) = axes('position', [16.5 20 3.5 3]);
panels{1}(2,1) = axes('position', [3  16 3.5 3]);
panels{1}(2,2) = axes('position', [7.5  16 3.5 3]);
panels{1}(2,3) = axes('position', [12 16 3.5 3]);
panels{1}(2,4) = axes('position', [16.5 16 3.5 3]);
panels{2}(1) = axes('position', [3 9 3 4]);
panels{2}(2) = axes('position', [7.5 9 3 4]);
panels{2}(3) = axes('position', [12 9 3 4]);
panels{2}(4) = axes('position', [16.5 9 3 4]);

%% properties to plot
features_names = {'duration';'compression';'distance';'distance_norm';};
xlable_strs = {
    'Replay duration (s)';
    {'Compression ratio';'(replay speed / flight speed)'};
    'Replay distance (m)';
    {'Replay distance','(norm. to environment size)'};
    };

%% arrange sessions to load
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
novelty_exposure_bats = [184 2382];
T(~ismember(T.bat_num, novelty_exposure_bats),:)=[];
exp_list = T.exp_ID;
early_days_exp_list = {
    'b0184_d191129',
    'b0184_d191130',
    'b0184_d191201',
    'b2382_d190623',
    'b2382_d190624',
    'b2382_d190627',
    };
early_days_IX = find(ismember(T.exp_ID, early_days_exp_list));
late_days_IX = find(~ismember(T.exp_ID, early_days_exp_list));
early_late_IX = {early_days_IX,late_days_IX};

%% load data
events_all_per_session = {};
flight_speed_all = zeros(size(exp_list));
for ii_epoch_type = 1:length(epoch_types)
    for ii_exp = 1:length(exp_list)
        exp_ID = exp_list{ii_exp};
        exp = exp_load_data(exp_ID,'details','flight');
        flight_speed_all(ii_exp) = mean(abs([exp.flight.speed_traj.vel_median_median]));
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
            events = [events_all_per_session{ii_epoch_type}{early_late_IX{ii_EL}}];
            seqs = [events.seq_model];
            X = [seqs.(fn)];
            hh(ii_epoch_type) = histogram(X,'DisplayStyle','stairs','Normalization','pdf','EdgeColor',epoch_type_clrs{ii_epoch_type},'LineWidth',lw,'BinWidth',panels_bin_size(ii_fn));
            text(0.5,0.9-ii_epoch_type*0.1, "{\itn}_{" + epoch_types_str{ii_epoch_type} + "} = "+length(X),'units','normalized','FontSize',7)
        end
        linkprop(hh,'BinEdges');
        xlim(panels_xlim(ii_fn,:))
        xticks(panels_xticks{ii_fn})
        if ii_EL == 2
            xlabel(xlable_strs{ii_fn});
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
            events = [events_all_per_session{ii_epoch_type}{IX}];
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
    ylabel(xlable_strs{ii_fn},'units','normalized','position',[ylabel_pos_x(ii_fn) 0.5]);
    xticks(1:4);
    xticklabels(["Sleep, days 1-5", "Sleep, days 6+","Awake, days 1-5", "Awake, days 6+"]);
    xtickangle(50);    
end

%% add sleep/awake legend
if exist('panels_hist_legend','var')
    delete(panels_hist_legend);
end
panels_hist_legend = axes('position', [3.8 23 0.3 0.25]);
cla
hold on
plot([0 1], [1 1], 'color', epoch_type_clrs{1}, 'LineWidth',lw,'Clipping','off');
plot([0 1], [0 0], 'color', epoch_type_clrs{2}, 'LineWidth',lw,'Clipping','off');
text(1.3, 1, 'Sleep','FontSize',7,'HorizontalAlignment','left');
text(1.3, 0, 'Awake','FontSize',7,'HorizontalAlignment','left');
xlim([0 1])
ylim([0 1])
axis off

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