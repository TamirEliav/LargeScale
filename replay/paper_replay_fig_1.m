%% Replay - Fig 1 - Behavior and replay examples
%%
% clear 
clc
close all

%% plotting options
err_dist_normalization = "cdf";
% err_dist_normalization = "probability";

% replay examples options
replay_examples_list = {
    {'sleep',11,'b0184_d191203',34 }
    {'sleep',11,'b0184_d191130',139}
    {'sleep',11,'b0184_d191202',60 }
    {'sleep',11,'b0184_d191203',31 }
    {'sleep',11,'b0184_d191204',45 }
    
    {'sleep',11,'b0184_d191203',34 }
    {'sleep',11,'b0184_d191130',139}
    {'sleep',11,'b0184_d191202',60 }
    {'sleep',11,'b0184_d191203',31 }
    {'sleep',11,'b0184_d191204',45 }
    };
replay_examples_list = cellfun(@(c)cell2struct(c,{'epoch_type','params_opt','exp_ID','event_num'},2), replay_examples_list)

%% graphics params
sleep_clr = [.6 .1 .8];
rest_clr = [.1 .8 .1];
flight_clr = 'r';

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Fig_1';
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
panels{1}(1) = axes('position', [2 23 17 2]);
panels{2}(1) = axes('position', [2 19.3 2.2 2.5]);
panels{2}(2) = axes('position', [4.8 19.3 2.2 2.5]);
panels{3}(1) = axes('position', [8.2 19.3 2.5 2.5]);
panels{4}(1) = axes('position', [12.7 19.3 3.3 2.5]);
panels{4}(2) = axes('position', [17 19.3 2 2.5]);
offsets_x = linspace(0,14,5)+2;
offsets_y = linspace(0,6,2)+7;
offsets_y = flip(offsets_y);
for ii=1:2
    for jj=1:5
        offset_x = offsets_x(jj);
        offset_y = offsets_y(ii);
        offset = [offset_x offset_y];
        panels{5}(jj,ii,1) = axes('position', [offset+[0 0] 2.3 3.5]);
        panels{5}(jj,ii,2) = axes('position', [offset+[0 3.5] 2.3 .5]);
        panels{5}(jj,ii,3) = axes('position', [offset+[0 4] 2.3 .5]);
    end
end
panels{6}(1,1) = axes('position', [2 1.5 4 4]);
panels{6}(2,1) = axes('position', [8 1.5 4 4]);
panels{6}(1,2) = axes('position', [4.7 4 1.5 1.5]);
panels{6}(2,2) = axes('position', [10.7 4 1.5 1.5]);
panels{6}(1,3) = axes('position', [2.35 4.8 .3 .4]);
panels{6}(2,3) = axes('position', [8.35 4.8 .3 .4]);

total_offset = [.5 1];
for ii = 1:length(panels)
    subpanels = panels{ii};
    subpanels = subpanels(:);
    for jj = 1:length(subpanels)
        subpanels(jj).Position([1 2]) = subpanels(jj).Position([1 2]) + total_offset;
    end
end

%% panel A
axes(panels{1})
cla reset
hold on
exp_ID = 'b0184_d191208';
exp = exp_load_data(exp_ID,'details','pos','flight','rest');
t = exp.pos.proc_1D.ts;
pos = exp.pos.proc_1D.pos_csaps;
nanpos = isnan(pos);
pos = interp1(t(~nanpos),pos(~nanpos),t,'linear','extrap');
win_s = 3;
fs = exp.pos.proc_1D.fs;
win_samples = round(win_s*fs);
% pos = smoothdata(pos,2,"movmedian",win_samples);
plot(t,pos,'-k','LineWidth',1)
m1 = min(pos);
m2 = max(pos);
ylimits = [m1 m2];
ylimits = ylimits + [-1 1].*range(ylimits)*0.01;
ylim(ylimits)

% add shaded area for different epochs
sleep_ti = exp_get_sessions_ti(exp_ID,'Sleep1','Sleep2');
rest_ti = exp.rest.ti;
flight_ti = [exp.flight.FE.start_ts; exp.flight.FE.end_ts]';
rest_ti(end,end) = 61550899128; % TODO: temp, need to fix the ts
sleep_ti(end,1) = 61604440361; % TODO: temp, need to fix the ts
epochs_ti = [sleep_ti; rest_ti];
epochs_clr = {repelem({sleep_clr},size(sleep_ti,1))
              repelem({rest_clr},size(rest_ti,1))};
epochs_clr = [epochs_clr{:}];
ylimits = get(gca,'YLim');
for ii_epoch = 1:size(epochs_ti,1)
    epoch_ti = epochs_ti(ii_epoch,:);
    clr = epochs_clr{ii_epoch};
    area(epoch_ti,ylimits([2 2]),'FaceColor',clr,'FaceAlpha',0.1,'EdgeColor','none','ShowBaseLine','off');
end
xticks([])
xlim(sleep_ti([1 end]));%+60e6.*[-1 1])
rescale_plot_data('x',[1e-6/60 0]); % change units to minutes
axis off
box off
% add scale bars
scale_bar_time_min = 2;
scale_bar_size_m = 50;
hax=gca;
xlimits = hax.XLim;
ylimits = hax.YLim;
x = xlimits(1) + 0.025*range(xlimits);
y = ylimits(1) + 0.35*range(ylimits);
plot(x+[0 scale_bar_time_min], y+[0 0],'k-','LineWidth',2,'Clipping','off');
plot(x+[0 0], y+[0 scale_bar_size_m],'k-','LineWidth',2,'Clipping','off');
text(x+scale_bar_time_min/2, y, sprintf('%dmin',scale_bar_time_min) ,'FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
text(x, y+scale_bar_size_m/2, sprintf('%dm',scale_bar_size_m) ,'FontSize',7,'HorizontalAlignment','center','VerticalAlignment','bottom','Rotation',90);
hax.Clipping='off';
ah = annotation('textarrow','headStyle','cback1','HeadLength',4,'HeadWidth',4);
ah.Parent = gca;
ah.Position = [xlimits(1) ylimits(1) 0.05*range(xlimits) 0];
ah = annotation('textarrow','headStyle','cback1','HeadLength',4,'HeadWidth',4);
ah.Parent = gca;
ah.Position = [xlimits(1) ylimits(1) 0 0.4*range(ylimits)];
text(xlimits(1)+0.00*range(xlimits), ylimits(1)-0.15*range(ylimits),'time in session','FontSize',8,'Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle');
text(xlimits(1)-0.02*range(xlimits), ylimits(1)+0.*range(ylimits),'Position','FontSize',8,'Rotation',90,'HorizontalAlignment','left','VerticalAlignment','middle');

%% load data for panels B and C
exp_ID = 'b0184_d191208';
exp = exp_load_data(exp_ID,'details','path','flight');
dec_param_opt = 4;
decode = decoding_load_data(exp_ID, 'flight', dec_param_opt);
CM = load(fullfile('F:\sequences\decoded_figs\flight\conf_mat\',sprintf('%s_flight_decoding_opt_%d.mat',exp_ID,dec_param_opt)));

%% panel B - single flight examples
FEs = exp.flight.FE;
FEs([FEs.distance]<100)=[];
examples_flight_IX = [27 28];
for ii=1:2
    axes(panels{2}(ii));
    cla
    hold on
    FE = FEs(examples_flight_IX(ii));
    ti = [FE.start_ts FE.end_ts];
    IX = get_data_in_ti(decode.time, ti);
    decode.MAP_pos(IX);
    plot([FE.ts], [FE.pos],'k','LineWidth',3);
    plot(decode.time(IX), decode.MAP_pos(IX),'r','LineWidth',.5);
    rescale_plot_data('x',[1e-6 ti(1)])
    xlabel('Time (s)','Units','normalized','Position',[0.5 -0.15]);
    if ii==1
        ylabel('Position (m)','Units','normalized','Position',[-0.25 0.5]);
    end
    xlim([0 20])
    ylim([0 150])
    xticks([0:10:20])
    yticks([0:50:150])
    hax=gca;
    hax.XRuler.TickLength(1) = 0.035;
    hax.YRuler.TickLength(1) = 0.035;
    hax.XRuler.TickLabelGapOffset = -1;
    hax.YRuler.TickLabelGapOffset = 0;
    title("Flight "+ii,'FontWeight','normal')
end

%% panel C - single sessions confusion matrix example
x = CM.res_raw.pos_real;
y = CM.res_raw.pos_predict;
bin_centers = decode.pos;
[bin_edges,bin_size] = centers2edges(bin_centers);
N = histcounts2(x,y,bin_edges,bin_edges)'; % transpose so yaxis(rows)=predict and xaxis(cols)=real
N_norm_by_real = N ./ sum(N,1);
% N_norm_by_predict = N ./ sum(N,2);
axes(panels{3});
cla
hold on
imagesc(bin_centers,bin_centers,N_norm_by_real);
colormap(gca,flip(colormap('bone')));
hcb = colorbar('Location','eastoutside');
hcb.Units = 'centimeters';
hcb.Position = [panels{3}.Position([1 2]) + [panels{3}.Position([3]) 0] + [0.1 0] ...
                0.25 panels{3}.Position(4)];
hcb.Limits = [0 1];
hcb.Label.String = 'Probability';
hcb.Label.Rotation = -90;
hcb.Label.Position = [2.5 0.5 0];
hcb.Label.FontSize = 8;
hcb.Ticks = [0 1];
hcb.TickLength = 0.015;
hcb.TickDirection='out';
axis xy
axis square
axis equal
xlim(bin_edges([1 end]));
ylim(bin_edges([1 end]));
pos_ticks = 0:20:200;
xticks(pos_ticks)
yticks(pos_ticks)
xlabel('Real position (m)','Units','normalized','Position',[0.5 -0.25]);
ylabel('Decoded position (m)','Units','normalized','Position',[-0.25 0.5]);
title("Example session",'FontWeight','normal')
hax=gca;
hax.TickDir = 'out';
hax.XRuler.TickLength(1) = 0.035;
hax.YRuler.TickLength(1) = 0.035;
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = 0;
hax.XRuler.TickLabelRotation = 45;

%% load data - error CDF plots (all sessions)
if ~exist('pos_err_median_all','var')
    [exp_list,T] = decoding_get_inclusion_list();
    T = T(exp_list,:);
    err_CDF_edges = 0:0.1:100;
    err_prob_edges = [0:0.2:5 inf];
    err_CDF = nan(length(err_CDF_edges)-1, length(exp_list));
    err_prob = nan(length(err_prob_edges)-1, length(exp_list));
    pos_err_median_all = [];
    for ii_exp = 1:length(exp_list)
        exp_ID = exp_list{ii_exp};
        filename = fullfile('F:\sequences\decoded_figs\flight\conf_mat\', ...
                            sprintf('%s_flight_decoding_opt_4.mat', exp_ID));
        CM = load(filename);
        err_CDF(:,ii_exp) = histcounts(CM.res_raw.pos_err,'Normalization','cdf','BinEdges',err_CDF_edges);
        err_prob(:,ii_exp) = histcounts(CM.res_raw.pos_err,'Normalization','probability','BinEdges',err_prob_edges);
        pos_err_median_all(ii_exp) = CM.res.pos_err_median;
    end
    err_CDF = [zeros(1,size(err_CDF,2));err_CDF]; % add zeros (cdf should start from 0)
end

%% panel D - error curves (CDF)
if err_dist_normalization == "cdf"
axes(panels{4})
hax=gca;
cla
hold on
plot(err_CDF_edges, err_CDF,'LineWidth',0.3);
shadedErrorBar(err_CDF_edges, err_CDF', {@mean,@nansem},'lineprops',{'k','linewidth',3});
xlabel('Positional decoding error (m)','Units','normalized', Position=[0.5 -0.22]);
title("{\itn} = "+size(err_CDF,2)+" sessions",'FontWeight','normal');
grid on
ylabel('Cumulative fraction')
hax.XScale = 'log';
% hax.XScale = 'linear';
hax.XTick = [.1 1 10 100];
hax.TickDir = 'out';
hax.XRuler.TickLength(1) = 0.035;
hax.YRuler.TickLength(1) = 0.035;
hax.XRuler.TickLabelGapOffset = -1.5;
hax.YRuler.TickLabelGapOffset = 0;
end

%% panel D - error curves (probability)
if err_dist_normalization == "probability"
axes(panels{4})
hax=gca;
cla
hold on
err_prob_centers = edges2centers(err_prob_edges);
err_prob_centers(end) = err_prob_centers(end-1) + median(diff(err_prob_centers));
xxx = err_prob_edges(1:end-1);
yyy = err_prob;
plot(xxx,yyy,'LineWidth',0.3);
shadedErrorBar(xxx, yyy', {@nanmean,@nansem},'lineprops',{'k','linewidth',3});
xlabel('Positional decoding error (m)','Units','normalized', Position=[0.5 -0.22]);
title("{\itn} = "+size(err_CDF,2)+" sessions",'FontWeight','normal');
grid on
ylabel('Probability')
hax.XScale = 'linear';
% hax.YScale = 'log';
hax.YScale = 'linear';
hax.XTick = [0:1:5];
hax.XLim = [0 7];
hax.TickDir = 'out';
hax.XRuler.TickLength(1) = 0.035;
hax.YRuler.TickLength(1) = 0.035;
hax.XRuler.TickLabelGapOffset = -1.5;
hax.YRuler.TickLabelGapOffset = 0;
end

%% median error hist
axes(panels{4}(2))
hax=gca;
cla
hold on
h = histogram(pos_err_median_all);
h.BinLimits = [0 3];
h.BinWidth = 0.2;
h.FaceColor = [1 1 1].*0.5;
hax=gca;
hax.XTick = [0:1:5];
hax.YTick = hax.YTick([1 end]);
hax.XLim = [0 3];
hax.TickDir = 'out';
hax.XRuler.TickLength(1) = 0.035;
hax.YRuler.TickLength(1) = 0.035;
hax.XRuler.TickLabelGapOffset = -1.5;
hax.YRuler.TickLabelGapOffset = 0;
xlabel({'Median decoding';'error per session (m)'},'Units','normalized', Position=[0.5 -0.15]);
ylabel('Counts','Units','normalized', Position=[-0.15 0.5]);

%% some stats to report
fprintf('median error: %.2g+-%.2g (mean=-std)\n',mean(pos_err_median_all),std(pos_err_median_all));

for CDF_thr = [0.5 0.6 0.7 0.75 0.8 0.85 0.9 0.95]
    mean_CDF_error = err_CDF_edges(find(mean(err_CDF,2)>CDF_thr,1,'first'));
    fprintf('mean cumulative error at %d%%: %gm\n',CDF_thr*100,mean_CDF_error);
end



%% replay examples (panels E+F)
panels_ex = panels{5}
panels_ex = reshape(panels_ex, size(panels_ex,1)*size(panels_ex,2), size(panels_ex,3))
cmap = bone;
cmap = flipud(cmap);
for ii_ex = 1:size(panels_ex,1)
    
    %% load data
    ex = replay_examples_list(ii_ex);
    addFieldsToWorkspace(ex);
    decode = decoding_load_data(exp_ID, epoch_type, params_opt);
    exp = exp_load_data(exp_ID,'details','path','MUA','ripples');
    events = decoding_load_events_quantification(exp_ID,epoch_type,params_opt,"posterior");
    event = events(event_num);
    seq = event.seq_model;
    seq_ti = [event.start_ts event.end_ts];
    t0 = mean(seq_ti);
    ti = t0 + [-1 1].*win_s*1e6;
    
    TT = exp.ripples.stats.best_TT;
    [LFP.signal, LFP.ts, LFP.fs, LFP.params] = LFP_load(exp_ID,TT,'band','ripple','limits_ts',ti);
    LFP.avg_signal = nanmean(LFP.signal,[2 3]);

    %% plot LFP
    axes(panels_ex(ii_ex,3));
    cla
    hold on
    plot(LFP.ts, LFP.avg_signal,'k');
    xlim(seq_ti+[-1 1].*0.2*range(seq_ti))
    xticks([])
    yticks([])
    rescale_plot_data('x',[1e-6 seq_ti(1)]);
    axis off
    title(sprintf('%s_%s_%d',epoch_type,exp_ID,event_num),'Interpreter','none','FontWeight','normal','FontSize',6);
    
    %% plot posterior (state)
    axes(panels_ex(ii_ex,2));
    cla
    hold on
    IX = get_data_in_ti(decode.time, ti);
    prob_t = decode.time(IX);
    prob_state = squeeze(decode.posterior_state(event.state_num,IX));
    plot(prob_t, prob_state, 'k','LineWidth',2);
    hax = gca;
    hax.XLim = prob_t([1 end]);
    hax.YLim = [0 1];
    box on
    colormap(cmap);
    hax.XLim = seq_ti+[-1 1].*0.2*range(seq_ti);
    hax.TickDir = 'out';
    hax.TickLength = [0.02 0.02];
    hax.XRuler.TickLabelGapOffset = -4;
    rescale_plot_data('x',[1e-6 seq_ti(1)]);
        
    %% plot posterior (position)
    axes(panels_ex(ii_ex,1));
    cla
    hold on
    IX = get_data_in_ti(decode.time, ti);
    prob_t = decode.time(IX);
    prob_pos = squeeze(decode.posterior(:,event.state_num,IX));
    imagesc(prob_t, decode.pos, prob_pos);
    plot([seq.start_ts seq.end_ts],[seq.start_pos seq.end_pos],'-r','LineWidth',0.01);
    hax = gca;
    hax.CLim = quantile(prob_pos(:),[0.01 0.99]);
    hax.XLim = prob_t([1 end]);
    hax.YLim = [min(decode.pos) max(decode.pos)] + [-1 1].*median(diff(decode.pos))*1;
    box on
    colormap(cmap);
    hax.XLim = seq_ti+[-1 1].*0.2*range(seq_ti);
    hax.TickDir = 'out';
    hax.TickLength = [0.02 0.02];
    hax.XRuler.TickLabelGapOffset = -1;
    rescale_plot_data('x',[1e-6 seq_ti(1)]);
    xlabel('Time (s)')
    
    %% link x axes
    linkaxes(panels_ex(ii_ex,:),'x');

end
%%
axes(panels{5}(1,1,1))
ylabel('Position (m)','Units','normalized','Position',[-0.3 0.6]);
axes(panels{5}(1,2,1))
ylabel('Position (m)','Units','normalized','Position',[-0.3 0.6]);

%% replay-triggered MUA/ripple-power
replay_trig_data_file = "L:\paper_replay\figures\Fig_replay_trig_MUA_ripples_10000ms_bats_all.mat";
data = load(replay_trig_data_file);
data = data.data_to_save;
fns = {'ripple','mua'};
labels = {'Ripple power (z)','Multiunit firing rate (Hz)'};

for ii = 1:2

    fn = fns{ii};
    axes(panels{6}(ii,1))
    cla
    hold on
    shadedErrorBar(data.(fn).t, data.(fn).sleep.mean, data.(fn).sleep.sem,'lineprops',{'Color',sleep_clr});
    shadedErrorBar(data.(fn).t, data.(fn).rest.mean, data.(fn).rest.sem,'lineprops',{'Color',rest_clr});
    xlim([-10 10])
    xticks([-10:5:10])
    hax=gca;
    hax.XRuler.TickLength(1) = 0.02;
    hax.YRuler.TickLength(1) = 0.02;
    hax.XRuler.TickLabelGapOffset = -1;
    hax.YRuler.TickLabelGapOffset = 0;
    xlabel('Time from replay (s)','Units','normalized','Position',[0.5 -0.1])
    ylabel(labels{ii},'Units','normalized','Position',[-0.15 0.5])
    
    axes(panels{6}(ii,2))
    cla
    hold on
    shadedErrorBar(data.(fn).t, data.(fn).sleep.mean, data.(fn).sleep.sem,'lineprops',{'Color',sleep_clr});
    shadedErrorBar(data.(fn).t, data.(fn).rest.mean, data.(fn).rest.sem,'lineprops',{'Color',rest_clr});
    xlim([-2 2])
    xticks([-2:1:2])
    xtickangle(0)
    hax=gca;
    hax.XRuler.TickLength(1) = 0.04;
    hax.YRuler.TickLength(1) = 0.04;
    hax.XRuler.TickLabelGapOffset = -1;
    hax.YRuler.TickLabelGapOffset = -1;
    
    axes(panels{6}(ii,3))
    cla
    hold on
    axis off
    plot([0 1],[1 1],'Color',sleep_clr,'LineWidth',2);
    plot([0 1],[0 0],'Color',rest_clr,'LineWidth',2);
    text(1.2,1,'Sleep','FontSize',9)
    text(1.2,0,'Rest','FontSize',9)
    xlim([0 1])
    ylim([0 1])
    
end



%% add panel letters
font_size = 11;
axes(panels{1}(1));
text(-0.06,1.12, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{2}(1));
text(-0.5,1.17, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{3}(1))
text(-0.42,1.17, 'c', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{4}(1))
text(-0.35,1.17, 'd', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{5}(1,1,end));
text(-0.35,1.9, 'e', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{5}(1,2,end));
text(-0.35,1.9, 'f', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{6}(1,1))
text(-0.28,1.12, 'g', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{6}(2,1))
text(-0.28,1.12, 'h', 'Units','normalized','FontWeight','bold','FontSize',font_size);


%% print/save the figure
fig_name_out = fullfile(res_dir, sprintf('%s_%s',fig_name_str,err_dist_normalization));

print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% exportgraphics(gcf,fullfile(res_dir,'testexport.pdf'),'BackgroundColor','none','ContentType','vector');

disp('figure was successfully saved to pdf/tiff/fig formats');

%%

