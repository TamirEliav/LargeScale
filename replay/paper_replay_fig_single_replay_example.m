%% Replay - Fig 1 - single replay example
function paper_replay_fig_single_replay_example(exp_ID,epoch_type,params_opt,event_num,win_s)
arguments
    %%
    exp_ID = 'b0184_d191205';
    epoch_type = 'sleep';
    params_opt = 11;
    event_num = 98;
    win_s = 0.5;
end

%% define output files
res_dir =  'L:\paper_replay\figures\Fig_replay_examples';
mkdir(res_dir)
fig_name_str = sprintf('Fig_replay_example_%s_%d_%s_event_%d',epoch_type,params_opt,exp_ID,event_num);

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
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

% create panels
clear panels
panels(1) = axes('position', [8 18 3 1]);
panels(2) = axes('position', [8 15 3 3]);

%% load data
decode = decoding_load_data(exp_ID, epoch_type, params_opt );
exp = exp_load_data(exp_ID,'details','path','MUA','ripples');
TT = exp.ripples.stats.best_TT;
[LFP.signal, LFP.ts, LFP.fs, LFP.params] = LFP_load(exp_ID,TT,'band','ripple');
LFP.avg_signal = nanmean(LFP.signal,[2 3]);
events = decoding_load_events_quantification(exp_ID,epoch_type,params_opt,"posterior");
event = events(event_num);
seq = event.seq_model;
ti = [event.start_ts event.end_ts];
t0 = mean(ti);
ti = mean(ti) + [-1 1].*win_s*1e6;

%% plot LFP
axes(panels(1));
cla
hold on
plot(LFP.ts, LFP.avg_signal,'k');
xlim(ti)
xticks([])
yticks([])
rescale_plot_data('x',[1e-6 t0]);
axis off
title(sprintf('%s_%s_%d',epoch_type,exp_ID,event_num),'Interpreter','none');

%% plot posterior (position)
axes(panels(2));
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
cmap = bone;
cmap = flipud(cmap);
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
linkaxes(panels(:),'x');

%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);

saveas(fig,fig_name_out,'pdf');
saveas(fig,fig_name_out,'jpeg');
disp('figure was successfully saved to pdf/tiff/fig formats');


