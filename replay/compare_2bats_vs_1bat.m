%%
clear
clc

%% data options 
params_opt = 11; % decoding opt 

%% plotting options
panels_xlim = [-0.0950    1.9950; -1.9500   40.9500; -0.4500   56.5500; -0.0210    0.4410];
panels_xticks = {[0 0.5 1 1.5],[0 20 40],[0:10:50],[0 0.2 0.4]};

%% graphics params
epoch_types = {'sleep','rest'};
exp_types = {'1 bat','2 bats'};
epoch_type_clrs = {[.6 .1 .8],[.1 .8 .1]};
line_styles = {'-','--'};

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'comparing_2bats_vs_1bat_data';
fig_caption_str = 'Two bats vs 1 bats experiemtn - replay properties comparison';
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

ii = 1;
y = 23;
x = [3 7.5 12 16.5];
w = 3; h = 3;
panels{ii}(1) = axes('position', [x(1) y w h]);
panels{ii}(2) = axes('position', [x(2) y w h]);
panels{ii}(3) = axes('position', [x(3) y w h]);
panels{ii}(4) = axes('position', [x(4) y w h]);

ii = 2;
y = 18;
x = [3 7.5 12 16.5];
w = 3; h = 3;
panels{ii}(1) = axes('position', [x(1) y w h]);
panels{ii}(2) = axes('position', [x(2) y w h]);
panels{ii}(3) = axes('position', [x(3) y w h]);
panels{ii}(4) = axes('position', [x(4) y w h]);

ii = 3;
y = 13;
x = [3 7.5 12 16.5];
w = 3; h = 3;
panels{ii}(1) = axes('position', [x(1) y w h]);
panels{ii}(2) = axes('position', [x(2) y w h]);
panels{ii}(3) = axes('position', [x(3) y w h]);
panels{ii}(4) = axes('position', [x(4) y w h]);

ii = 4;
y = 8;
x = [3 7.5 12 16.5];
w = 3; h = 3;
panels{ii}(1) = axes('position', [x(1) y w h]);
panels{ii}(2) = axes('position', [x(2) y w h]);
panels{ii}(3) = axes('position', [x(3) y w h]);
panels{ii}(4) = axes('position', [x(4) y w h]);

ii = 5;
y = 3;
x = [3 7.5 12 16.5];
w = 3; h = 3;
panels{ii}(1) = axes('position', [x(1) y w h]);
panels{ii}(2) = axes('position', [x(2) y w h]);
panels{ii}(3) = axes('position', [x(3) y w h]);
panels{ii}(4) = axes('position', [x(4) y w h]);

%% load data - 1 bat / 2 bats experiments
epoch_types = {'sleep','rest'};
exp_list_1bat = decoding_get_inclusion_list();
exp_list_2bats = {
    'b2299_d191202',
    'b2299_d191203',
    'b2299_d191204',
    'b2299_d191205',
    'b2299_d191208',
    'b2299_d191209',
    'b2299_d191210',
    'b2299_d191213',
    };
exp_lists = {exp_list_1bat,exp_list_2bats};

events_all_per_session = {};
for ii_exp_type = 1:length(exp_types)
    exp_list = exp_lists{ii_exp_type};
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
            events_all_per_session{ii_exp_type,ii_epoch_type}{ii_exp} = events;
        end
    end
end

%%
features_names = {'duration';'compression';'distance';'distance_norm';};
xlable_strs = {
    'Replay duration (s)';
    {'Compression ratio';'(replay speed / flight speed)'};
    'Replay distance (m)';
    {'Replay distance','(norm. to environment size)'};
    };

%% 2bats replay inclusion for awake replays
data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115.mat';
data_2bats_rest_inclusion = load(data_filename);

%% plot boxplots - apply 2 bats inclusino criteria (only for awake)
boxplot_panels_ylimits = [0 0.8; 0 25; 3 23; .02 .18];
pvals = [];
for ii_fn = 1:length(features_names)
    axes(panels{1}(ii_fn));
    cla
    hold on
    fn = features_names{ii_fn};
    X={};
    G={};
    for ii_epoch_type = 1:length(epoch_types)
        for ii_exp_type = 1:length(exp_types)
            epoch_type = epoch_types{ii_epoch_type};
            exp_type = exp_types{ii_exp_type};
            events = [events_all_per_session{ii_exp_type,ii_epoch_type}{:}];
            if strcmp(epoch_type,'rest') && strcmp(exp_type,'2 bats')
                events(~data_2bats_rest_inclusion.TF) = [];
            end
            seqs = [events.seq_model];
            x = [seqs.(fn)];
            g = (ii_exp_type+2*(ii_epoch_type-1));
            X{ii_epoch_type,ii_exp_type} = x;
            G{ii_epoch_type,ii_exp_type} = ones(size(seqs)).*g;
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
%         hax= gca;
%         yy = hax.YLim([2 2]);
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
    ylabel(xlable_strs{ii_fn},'units','normalized','position',[-0.23 0.5]);
    xticks(1:4);
    sdf = string(epoch_types) + ", " + string(exp_types)';
    sdf = strrep(sdf(:),'rest','awake');
    xticklabels(sdf);
    xtickangle(50);    
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
        for ii_exp_type = 1:length(exp_types)
            events = [events_all_per_session{ii_exp_type,ii_epoch_type}{:}];
            seqs = [events.seq_model];
            x = [seqs.(fn)];
            g = (ii_exp_type+2*(ii_epoch_type-1));
            X{ii_epoch_type,ii_exp_type} = x;
            G{ii_epoch_type,ii_exp_type} = ones(size(seqs)).*g;
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
%         hax= gca;
%         yy = hax.YLim([2 2]);
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
    ylabel(xlable_strs{ii_fn},'units','normalized','position',[-0.23 0.5]);
    xticks(1:4);
    sdf = string(epoch_types) + ", " + string(exp_types)';
    sdf = strrep(sdf(:),'rest','awake');
    xticklabels(sdf);
    xtickangle(50);    
end

%% plot boxplots - pool sleep/rest
boxplot_panels_ylimits = [0 0.8; 0 25; 3 23; .02 .18];
pvals = [];
for ii_fn = 1:length(features_names)
    axes(panels{3}(ii_fn));
    cla
    hold on
    fn = features_names{ii_fn};
    X={};
    G={};
    for ii_exp_type = 1:length(exp_types)
        events1 = [events_all_per_session{ii_exp_type,1}{:}];
        events2 = [events_all_per_session{ii_exp_type,2}{:}];
        seqs1 = [events1.seq_model];
        seqs2 = [events2.seq_model];
        seqs = [seqs1 seqs2];
        x = [seqs.(fn)];
        g = ii_exp_type;
        X{ii_exp_type} = x;
        G{ii_exp_type} = ones(size(seqs)).*g;
        m1 = prctile(x,50);
        m2 = prctile(x,[25 75]);
        m3 = prctile(x,[5 95]);
        w = 0.2;
        lw = 1.3;
        clr = 'k';
        plot([g-w g+w],[m1 m1],'Color',clr,'LineWidth',lw);
        plot([g g],m3,'Color',clr,'LineWidth',lw);
        rectangle('Position',[g-w m2(1) 2*w diff(m2)],'EdgeColor',clr,'FaceColor','none','LineWidth',lw);
    end
    pval = ranksum(X{1},X{2});
    str = genSignifStrAstricks(pval);
    xx = [g g-1];
    yy = boxplot_panels_ylimits(ii_fn,[2 2]);
%         hax= gca;
%         yy = hax.YLim([2 2]);
    plot(xx,yy,'k-');
    font_size = 10;
    if strcmp(str,'n.s.')
        font_size = 8;
        yy = yy.*1.02;
    end
    text(mean(xx),mean(yy),str,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',font_size);
    pvals(ii_fn,ii_epoch_type) = pval;
    xlim([0.2 2.8])
    ylim(boxplot_panels_ylimits(ii_fn,:))
    xlabel('');
    ylabel(xlable_strs{ii_fn},'units','normalized','position',[-0.23 0.5]);
    xticks(1:length(exp_types));
    xticklabels(exp_types);
    xtickangle(50);    
end

%% plot boxplots - pool sleep/rest - compare only 130m (remove 200m data)
boxplot_panels_ylimits = [0 0.8; 0 25; 3 23; .02 .18];
pvals = [];
for ii_fn = 1:length(features_names)
    axes(panels{4}(ii_fn));
    cla
    hold on
    fn = features_names{ii_fn};
    X={};
    G={};
    for ii_exp_type = 1:length(exp_types)
        events1 = [events_all_per_session{ii_exp_type,1}{:}];
        events2 = [events_all_per_session{ii_exp_type,2}{:}];
        seqs1 = [events1.seq_model]; [seqs1.recordingArena] = disperse({events1.recordingArena});
        seqs2 = [events2.seq_model]; [seqs2.recordingArena] = disperse({events2.recordingArena});
        seqs = [seqs1 seqs2];
        TF_200m = {seqs.recordingArena}=="200m";
        seqs(TF_200m)=[];
        x = [seqs.(fn)];
        g = ii_exp_type;
        X{ii_exp_type} = x;
        G{ii_exp_type} = ones(size(seqs)).*g;
        m1 = prctile(x,50);
        m2 = prctile(x,[25 75]);
        m3 = prctile(x,[5 95]);
        w = 0.2;
        lw = 1.3;
        clr = 'k';
        plot([g-w g+w],[m1 m1],'Color',clr,'LineWidth',lw);
        plot([g g],m3,'Color',clr,'LineWidth',lw);
        rectangle('Position',[g-w m2(1) 2*w diff(m2)],'EdgeColor',clr,'FaceColor','none','LineWidth',lw);
    end
    pval = ranksum(X{1},X{2});
    str = genSignifStrAstricks(pval);
    xx = [g g-1];
    yy = boxplot_panels_ylimits(ii_fn,[2 2]);
%         hax= gca;
%         yy = hax.YLim([2 2]);
    plot(xx,yy,'k-');
    font_size = 10;
    if strcmp(str,'n.s.')
        font_size = 8;
        yy = yy.*1.02;
    end
    text(mean(xx),mean(yy),str,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',font_size);
    pvals(ii_fn,ii_epoch_type) = pval;
    xlim([0.2 2.8])
    ylim(boxplot_panels_ylimits(ii_fn,:))
    xlabel('');
    ylabel(xlable_strs{ii_fn},'units','normalized','position',[-0.23 0.5]);
    xticks(1:length(exp_types));
    xticklabels(exp_types);
    xtickangle(50);    
end

%% histograms
for ii_fn = 1:length(features_names)
    fn = features_names{ii_fn};
    axes(panels{5}(ii_fn));
    cla
    hold on
    lw = 1.1;
    clear hh;
    for ii_exp_type = 2%1:length(exp_types)
        for ii_epoch_type = 1:length(epoch_types)
            events = [events_all_per_session{ii_exp_type,ii_epoch_type}{:}];
            seqs = [events.seq_model];
            X = [seqs.(fn)];
            hh(ii_exp_type,ii_epoch_type) = histogram(X,'DisplayStyle','stairs','Normalization','pdf','EdgeColor',epoch_type_clrs{ii_epoch_type},'LineWidth',lw);
            text(0.5,0.9-ii_epoch_type*0.1, "\itn_{" + epoch_types{ii_epoch_type} + "} = "+length(X),'units','normalized','FontSize',7)
        end
    end
    linkprop(hh(:),'BinEdges');
    xlim(panels_xlim(ii_fn,:))
    xticks(panels_xticks{ii_fn})
    xlabel(xlable_strs{ii_fn});
    if ii_fn == 1
        ylabel('Probability density','Units','normalized','Position',[-0.13 0.5]);
    end
    hax=gca;
%     hax.TickLength(1) = [0.025];
    hax.XRuler.TickLength(1) = 0.03;
    hax.XRuler.TickLabelGapOffset = -1.2;
    hax.YRuler.TickLabelGapOffset = 1;
    fprintf('%s: median %.2g, mean = %.2g\n',fn,median(X),mean(X))
end




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










%%
