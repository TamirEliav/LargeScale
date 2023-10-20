%% Replay - Fig supp 1 - replay example deconstructed
%%
clear 
clc
close all

%% plotting options
example_opt_IX = 10 % I chose option 5
example_options =[
struct(exp_ID = 'b0184_d191130', epoch_type = 'sleep', params_opt = 11, event_num = 139)
struct(exp_ID = 'b0184_d191201', epoch_type = 'sleep', params_opt = 11, event_num = 100)
struct(exp_ID = 'b0184_d191202', epoch_type = 'sleep', params_opt = 11, event_num = 204)
struct(exp_ID = 'b0184_d191202', epoch_type = 'sleep', params_opt = 11, event_num = 60)
struct(exp_ID = 'b0184_d191203', epoch_type = 'sleep', params_opt = 11, event_num = 31)
struct(exp_ID = 'b0184_d191203', epoch_type = 'sleep', params_opt = 11, event_num = 34)
struct(exp_ID = 'b0184_d191204', epoch_type = 'sleep', params_opt = 11, event_num = 45)
struct(exp_ID = 'b0184_d191205', epoch_type = 'sleep', params_opt = 11, event_num = 227)
struct(exp_ID = 'b0184_d191212', epoch_type = 'sleep', params_opt = 11, event_num = 56)
struct(exp_ID = 'b2382_d190712', epoch_type = 'sleep', params_opt = 11, event_num = 24)
];
example = example_options(example_opt_IX)
addFieldsToWorkspace(example);

%% graphics params
win_s = 1;
cmap = bone;
cmap = flipud(cmap);

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Fig_supp_1';
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
panels{1}(1) = axes('position', [4 23.8 8 1.5]);
panels{1}(2) = axes('position', [4 22.1 8 1.5]);
panels{1}(3) = axes('position', [4 20.5 8 1.5]);
panels{1}(4) = axes('position', [4 12.0 8 8]);
panels{1}(5) = axes('position', [4  3.5 8 8]);

%% load data
if ~exist('decode','var')
decode = decoding_load_data(exp_ID, epoch_type, params_opt );
exp = exp_load_data(exp_ID,'details','path','MUA','ripples');
events = decoding_load_events_quantification(exp_ID,epoch_type,params_opt,"posterior");
event = events([events.num] ==event_num);
seq = event.seq_model;
seq_ti = [event.start_ts event.end_ts];
t0 = mean(seq_ti);
ti = t0 + [-1 1].*win_s*1e6;

TT = exp.ripples.stats.best_TT;
[LFP.signal, LFP.ts, LFP.fs, LFP.params] = LFP_load(exp_ID,TT,'band','ripple','limits_ts',ti);
LFP.avg_signal = nanmean(LFP.signal,[2 3]);
end

%% rename decoder states
states = decode.state';
states(states=="Inbound-empirical_movement") = "Movement state dir 2";
states(states=="Inbound-identity") = "Stationary state dir 2";
states(states=="Inbound-uniform") = "Fragmented state dir 2";
states(states=="Outbound-empirical_movement") = "Movement state dir 1";
states(states=="Outbound-identity") = "Stationary state dir 1";
states(states=="Outbound-uniform") = "Fragmented state dir 1";

%% legend
if exist('panels_legend','var')
    delete(panels_legend);
end
panels_legend = axes('position', [13. 20.8 5 0.8]);
cla
hold on
lw=2;
t = [0 0.05];
x = [.62 .62 .62 0 0 0];
y = [1 .5 0 1 .5 0];
for ii=1:6
    plot(x(ii)+t,y([ii ii]),'-','LineWidth',lw,'Clipping','off');
    text(x(ii)+t(2)+0.02, y(ii), states(ii),'FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
end

xlim([0 1])
ylim([0 1])
axis off
% title('States legend:','Units','normalized','Position',[0.5 1.4])

%% plot MUA
axes(panels{1}(1));
cla
hold on
IX = get_data_in_ti(exp.MUA.t,ti);
x = exp.MUA.t(IX);
y = exp.MUA.FR(IX);
area(x,y,'FaceColor','k');
xlim(ti)
hax=gca;
m = roundToNearestValue(max(y), 20, @ceil)
hax.XTick = [];
hax.YTick = [0 m];
hax.YLim = [0 m]
rescale_plot_data('x',[1e-6 t0]);
% axis off
ylabel({'Multiunit';'firing-rate';'(Hz)'}, 'Units','normalized', 'Position',[-0.02 .5]);
title(sprintf('%s_%s_%d',epoch_type,exp_ID,event_num),'Interpreter','none','Units','normalized','Position',[0.5 1.4]);

%% plot LFP
axes(panels{1}(2));
cla
hold on
plot(LFP.ts, LFP.avg_signal,'k');
xlim(ti)
xticks([])
yticks([])
rescale_plot_data('x',[1e-6 t0]);
axis off

%% plot posterior (state)
axes(panels{1}(3));
cla
hold on
IX = get_data_in_ti(decode.time, ti);
prob_t = decode.time(IX);
h=plot(decode.time(IX), decode.posterior_state(:,IX),'LineWidth',1);
h(event.state_num).LineWidth=2.5; % highlight the event state
yline(0.8,'r--')
hax = gca;
hax.XLim = prob_t([1 end]);
hax.YLim = [0 1];
box off
colormap(cmap);
hax.XLim = seq_ti+[-1 1].*0.2*range(seq_ti);
hax.XTick = [-win_s 0 win_s];
hax.XTickLabel=[];
hax.YTick = [0 1];
hax.TickDir = 'out';
hax.TickLength = [0.02 0.02];
hax.XRuler.TickLabelGapOffset = -4;
rescale_plot_data('x',[1e-6 t0]);
ylabel({'State';'prob.'}, 'Units','normalized', 'Position',[-0.02 .5]);

%% plot posterior (position)
axes(panels{1}(4));
cla
hold on
IX = get_data_in_ti(decode.time, ti);
prob_t = decode.time(IX);
prob_pos = decode.posterior_pos(:,IX);
imagesc(prob_t, decode.pos, prob_pos);
plot([seq.start_ts seq.end_ts],[seq.start_pos seq.end_pos],'-r','LineWidth',1);
hax = gca;
hax.CLim = quantile(prob_pos(:),[0.01 0.99]);
hax.XLim = prob_t([1 end]);
hax.YLim = [min(decode.pos) max(decode.pos)] + [-1 1].*median(diff(decode.pos))*1;
box on
colormap(cmap);
alpha = 0.05;
% area([hax.XLim(1) event.start_ts],hax.YLim([2 2]), 'FaceColor','r','FaceAlpha',alpha);
% area([event.end_ts hax.XLim(2)],hax.YLim([2 2]),   'FaceColor','r','FaceAlpha',alpha);
% area([event.start_ts event.end_ts],hax.YLim([2 2]),   'FaceColor','g','FaceAlpha',alpha);
rescale_plot_data('x',[1e-6 t0]);
hax.XLim = [-1 1].*win_s;
hax.XTick = [-win_s 0 win_s];
hax.XTickLabel=[];
hax.YTick=[];
hax.TickDir = 'out';
hax.TickLength = [0.02 0.02];
hax.XRuler.TickLabelGapOffset = -4;
ylabel('Replay position (m)', 'Units','normalized', 'Position',[-0.07 .5]);

%% plot posterior (position)
axes(panels{1}(5));
cla
hold on
IX = get_data_in_ti(decode.time, ti);
prob_t = decode.time(IX);
prob_pos = squeeze(decode.posterior(:,event.state_num,IX));
imagesc(prob_t, decode.pos, prob_pos);
plot([seq.start_ts seq.end_ts],[seq.start_pos seq.end_pos],'-r','LineWidth',1);
hax = gca;
hax.CLim = quantile(prob_pos(:),[0.01 0.99]);
hax.XLim = prob_t([1 end]);
hax.YLim = [min(decode.pos) max(decode.pos)] + [-1 1].*median(diff(decode.pos))*1;
box on
colormap(cmap);
alpha = 0.05;
% area([hax.XLim(1) event.start_ts],hax.YLim([2 2]), 'FaceColor','r','FaceAlpha',alpha);
% area([event.end_ts hax.XLim(2)],hax.YLim([2 2]),   'FaceColor','r','FaceAlpha',alpha);
% area([event.start_ts event.end_ts],hax.YLim([2 2]),   'FaceColor','g','FaceAlpha',alpha);
rescale_plot_data('x',[1e-6 t0]);
hax.XLim = [-1 1].*win_s;
hax.XTick = [-win_s 0 win_s];
hax.YTick = [];
hax.TickDir = 'out';
hax.TickLength = [0.02 0.02];
hax.XRuler.TickLabelGapOffset = 0;
xlabel('Time (s)', 'Units','normalized', 'Position',[0.5 -0.08]);
ylabel('Replay position (m)', 'Units','normalized', 'Position',[-0.07 .5]);

%% link x axes
linkaxes(panels{1}(:),'x');

%% add panel letters
font_size = 11;
axes(panels{1}(1))
text(-0.13,1.1, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(2))
text(-0.13,0.9, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(3))
text(-0.13,1.1, 'c', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(4))
text(-0.13,1.0, 'd', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(5))
text(-0.13,1.0, 'e', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = sprintf('%s_opt_%d_%s_%s_%d_%d',fig_name_str,example_opt_IX,exp_ID,epoch_type,params_opt,event_num);
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%
