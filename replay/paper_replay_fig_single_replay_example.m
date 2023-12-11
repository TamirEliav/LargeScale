%% Replay - Fig 1 - single replay example
function paper_replay_fig_single_replay_example(exp_ID,epoch_type,params_opt,event_num,win_s,opts)
arguments
    %%
    exp_ID = 'b0184_d191205';
    epoch_type = 'sleep';
    params_opt = 11;
    event_num = 98;
    win_s = 0.5;
    opts.res_dir =  'L:\paper_replay\figures\Fig_replay_examples'
    opts.filename_prefix = ''
    opts.title_str_prefix = ''
end

%% define output files
mkdir(opts.res_dir)
fig_name_str = sprintf('Fig_replay_example_%s_%d_%s_event_%d',epoch_type,params_opt,exp_ID,event_num);
fig_name_str = [opts.filename_prefix fig_name_str];
title_str = sprintf('%s_%s_%d',epoch_type,exp_ID,event_num);
title_str = [opts.title_str_prefix title_str];

%% params
cmap = bone;
cmap = flipud(cmap);

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
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
% annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

% create panels
clear panels
panels{1}(1) = axes('position', [5 20 3 1]);
panels{1}(2) = axes('position', [5 16 3 4]);

panels{2}(1) = axes('position', [5 10.5 3 1]);
panels{2}(2) = axes('position', [5 10 3 .5]);
panels{2}(3) = axes('position', [5  6 3 4]);

panels{3}(1) = axes('position', [9 21.5 3 .5]);
panels{3}(2) = axes('position', [9 20.5 3 1]);
panels{3}(3) = axes('position', [9 20 3 .5]);
panels{3}(4) = axes('position', [9 16 3 4]);

panels{4}(1) = axes('position', [8.7 11.5 3 .5]);
panels{4}(2) = axes('position', [8.7 10.5 3 1]);
panels{4}(3) = axes('position', [8.7 10 3 .5]);
panels{4}(4) = axes('position', [8.7  6 3 4]);

%% load data
decode = decoding_load_data(exp_ID, epoch_type, params_opt );
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

%% ========================= version 1 ====================================
%% plot LFP
axes(panels{1}(1));
cla
hold on
plot(LFP.ts, LFP.avg_signal,'k');
xlim(ti)
xticks([])
yticks([])
rescale_plot_data('x',[1e-6 t0]);
axis off
title(title_str,'Interpreter','none');

%% plot posterior (position)
axes(panels{1}(2));
cla
hold on
IX = get_data_in_ti(decode.time, ti);
prob_t = decode.time(IX);
prob_pos = decode.posterior_pos(:,IX);
imagesc(prob_t, decode.pos, prob_pos);
plot([seq.start_ts seq.end_ts],[seq.start_pos seq.end_pos],'-r','LineWidth',0.01);
hax = gca;
hax.CLim = quantile(prob_pos(:),[0.01 0.99]);
hax.XLim = prob_t([1 end]);
hax.YLim = [min(decode.pos) max(decode.pos)] + [-1 1].*median(diff(decode.pos))*1;
box on
colormap(cmap);
alpha = 0.05;
area([hax.XLim(1) event.start_ts],hax.YLim([2 2]), 'FaceColor','r','FaceAlpha',alpha);
area([event.end_ts hax.XLim(2)],hax.YLim([2 2]),   'FaceColor','r','FaceAlpha',alpha);
area([event.start_ts event.end_ts],hax.YLim([2 2]),   'FaceColor','g','FaceAlpha',alpha);
rescale_plot_data('x',[1e-6 t0]);
hax.XLim = [-1 1].*win_s;
hax.XTick = [-win_s 0 win_s];
hax.TickDir = 'out';
hax.TickLength = [0.02 0.02];
hax.XRuler.TickLabelGapOffset = -4;

%% link x axes
linkaxes(panels{1}(:),'x');

%% ========================= version 2 ====================================
%% plot LFP
axes(panels{2}(1));
cla
hold on
plot(LFP.ts, LFP.avg_signal,'k');
xlim(seq_ti+[-1 1].*0.2*range(seq_ti))
xticks([])
yticks([])
rescale_plot_data('x',[1e-6 seq_ti(1)]);
axis off
title(title_str,'Interpreter','none');

%% plot posterior (state)
axes(panels{2}(2));
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
axes(panels{2}(3));
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
hax.XRuler.TickLabelGapOffset = -4;
rescale_plot_data('x',[1e-6 seq_ti(1)]);

%% link x axes
linkaxes(panels{2}(:),'x');


%% ========================= version 3 ====================================
%% plot MUA
axes(panels{3}(1));
cla
hold on
IX = get_data_in_ti(exp.MUA.t,ti);
x = exp.MUA.t(IX);
y = exp.MUA.FR(IX);
area(x,y,'FaceColor','k');
xlim(ti)
xticks([])
yticks([])
rescale_plot_data('x',[1e-6 t0]);
axis off
title(title_str,'Interpreter','none');

%% plot LFP
axes(panels{3}(2));
cla
hold on
plot(LFP.ts, LFP.avg_signal,'k');
xlim(ti)
xticks([])
yticks([])
rescale_plot_data('x',[1e-6 t0]);
axis off

%% plot posterior (state)
axes(panels{3}(3));
cla
hold on
IX = get_data_in_ti(decode.time, ti);
prob_t = decode.time(IX);
h=plot(decode.time(IX), decode.posterior_state(:,IX),'LineWidth',1);
h(event.state_num).LineWidth=2; % highlight the event state
hax = gca;
hax.XLim = prob_t([1 end]);
hax.YLim = [0 1];
box on
colormap(cmap);
hax.XLim = seq_ti+[-1 1].*0.2*range(seq_ti);
hax.TickDir = 'out';
hax.TickLength = [0.02 0.02];
hax.XRuler.TickLabelGapOffset = -4;
rescale_plot_data('x',[1e-6 t0]);



%% plot posterior (position)
axes(panels{3}(4));
cla
hold on
IX = get_data_in_ti(decode.time, ti);
prob_t = decode.time(IX);
prob_pos = decode.posterior_pos(:,IX);
imagesc(prob_t, decode.pos, prob_pos);
plot([seq.start_ts seq.end_ts],[seq.start_pos seq.end_pos],'-r','LineWidth',0.01);
hax = gca;
hax.CLim = quantile(prob_pos(:),[0.01 0.99]);
hax.XLim = prob_t([1 end]);
hax.YLim = [min(decode.pos) max(decode.pos)] + [-1 1].*median(diff(decode.pos))*1;
box on
colormap(cmap);
alpha = 0.05;
area([hax.XLim(1) event.start_ts],hax.YLim([2 2]), 'FaceColor','r','FaceAlpha',alpha);
area([event.end_ts hax.XLim(2)],hax.YLim([2 2]),   'FaceColor','r','FaceAlpha',alpha);
area([event.start_ts event.end_ts],hax.YLim([2 2]),   'FaceColor','g','FaceAlpha',alpha);
rescale_plot_data('x',[1e-6 t0]);
hax.XLim = [-1 1].*win_s;
hax.XTick = [-win_s 0 win_s];
hax.TickDir = 'out';
hax.TickLength = [0.02 0.02];
hax.XRuler.TickLabelGapOffset = -4;

%% link x axes
linkaxes(panels{3}(:),'x');


%% ========================= version 4 ====================================
%% plot MUA
axes(panels{4}(1));
cla
hold on
IX = get_data_in_ti(exp.MUA.t,ti);
x = exp.MUA.t(IX);
y = exp.MUA.FR(IX);
area(x,y,'FaceColor','k');
xticks([])
yticks([])
rescale_plot_data('x',[1e-6 t0]);
axis off
title(title_str,'Interpreter','none');

%% plot LFP
axes(panels{4}(2));
cla
hold on
plot(LFP.ts, LFP.avg_signal,'k');
xlim(seq_ti+[-1 1].*0.2*range(seq_ti))
xticks([])
yticks([])
rescale_plot_data('x',[1e-6 t0]);
axis off

%% plot posterior (state)
axes(panels{4}(3));
cla
hold on
IX = get_data_in_ti(decode.time, ti);
prob_t = decode.time(IX);
h=plot(decode.time(IX), decode.posterior_state(:,IX),'LineWidth',1);
h(event.state_num).LineWidth=2; % highlight the event state
hax = gca;
hax.XLim = prob_t([1 end]);
hax.YLim = [0 1];
box on
colormap(cmap);
hax.XLim = seq_ti+[-1 1].*0.2*range(seq_ti);
hax.TickDir = 'out';
hax.TickLength = [0.02 0.02];
hax.XRuler.TickLabelGapOffset = -4;
rescale_plot_data('x',[1e-6 t0]);


%% plot posterior (position)
axes(panels{4}(4));
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
hax.XRuler.TickLabelGapOffset = -4;
rescale_plot_data('x',[1e-6 t0]);

%% link x axes
linkaxes(panels{4}(:),'x');
xlim([-1 1].*win_s)



%% ============ print/save the figure==============
fig_name_out = fullfile(opts.res_dir, fig_name_str);

saveas(fig,fig_name_out,'pdf');
saveas(fig,fig_name_out,'jpeg');
% saveas(fig,fig_name_out,'meta');
disp('figure was successfully saved to pdf/tiff/fig formats');
close(fig)


