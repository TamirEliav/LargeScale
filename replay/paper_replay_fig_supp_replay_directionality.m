%% Replay - Fig supp 9 - replay directionality
clear 
clc
close all

%% data options 
params_opt = 11; % decoding opt 
novelty_session_num_thr = 5;
minSeqsThr = 30;

%% plotting options

%% graphics params
epoch_types = {'sleep','rest'};
epoch_type_clrs = {[.6 .1 .8],[.1 .8 .1]};
directions_clrs = {[0    0.3843    0.7451];[ 0.5216    0.2471         0]};

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Extended_Data_Fig_9';
fig_caption_str = 'Replay directinoality';
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
% panels{1}(2,1) = axes('position', [2 4 3.5 8]);
% panels{1}(2,2) = axes('position', [6.5 4 3.5 8]);
panels{2}(1) = axes('position', [12 19.4 4 3.5]);

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
        n1 = sum([seqs.direction]==-1);
        n2 = sum([seqs.direction]==1);
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
        'Over-representation of replay directionality'};
    text(-.2, 1.12, title_str,'FontSize',9,'HorizontalAlignment','center','Units','normalized')
end

%% replay directionality bias - plot trend over exposure to enviroenment
axes(panels{2}(1));
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
% legend
x = [0.3 2]+32;
y = linspace(0.35,0,3)+0.6;
plot(x,y(1)*[1 1],'Color',clrs{1},'LineWidth',1.5,'Clipping','off')
plot(x,y(2)*[1 1],'Color',clrs{2},'LineWidth',1.5,'Clipping','off')
plot(x,y(3)*[1 1],'Color',clrs{3},'LineWidth',1.5,'Clipping','off')
x = x(end)+1;
text(x,y(1), "Sleep", 'FontSize',7)
text(x,y(2), "Awake", 'FontSize',7)
text(x,y(3), "Pooled", 'FontSize',7)
xlabel('Session no.','Units','normalized','Position',[0.5 -0.05]);
ylabel('Replay directionality index')
ylim([-1 1])
xticks([1 40])
hax=gca;
hax.XRuler.TickLength(1) = 0.02;
hax.YRuler.TickLength(1) = 0.035;
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = 0;
% text(5.5,-0.75,"\leftarrow"+"session #"+example_session_num_form_exposure,'FontSize',8)
h=annotation('textarrow');
h.Parent=hax;
h.X = 3.5*[1 1];
h.Y = [1 0.85];
h.String = {'                  Sessions #3 and #4'};
h.FontSize = 7; h.HeadLength = 4; h.HeadWidth = 4;

%% add panel letters
font_size = 11;
axes(panels{1}(1))
text(-0.3,1.1, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{2}(1))
text(-0.3,1.1, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = sprintf('%s_decoding_opt_%d',fig_name_str, params_opt);
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%
