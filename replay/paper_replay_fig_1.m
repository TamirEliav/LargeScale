%% Replay - Fig 1 - Behavior and replay examples
%%
clear 
clc

%% plotting options
err_dist_normalization = 'cdf';
% err_dist_normalization = 'pdf';

%% graphics params
sleep_clr = [.6 .1 .8];
rest_clr = [.1 .8 .1];
flight_clr = 'r';

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Fig_behavior_decoding';
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
panel_A = axes('position', [1 22 19 1.5]);
panel_B = [  axes('position', [5 17 3 2.5]) ...
             axes('position', [12 17 3 2.5])];
panel_C = axes('position', [5 10 4 4]);
panel_D = axes('position', [12 10 4 3]);
% panel_E = [];
% panel_E_pos = [3 7];
% for ii=1:2
%     for jj=1:5
%         offset_x = (jj-1)*3.5;
%         offset_y = (ii-1)*5.5;
%         offset = panel_E_pos + [offset_x offset_y];
%         panel_E(ii,jj,1) = axes('position', [offset+[0 0] 3 3]);
%         panel_E(ii,jj,2) = axes('position', [offset+[0 3] 3 1]);
%     end
% end
% panel_GH = [  axes('position', [3 2 3 3]) ...
%              axes('position', [8 2 3 3])];

%% panel A
axes(panel_A)
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

% add shaded area for different epochs
sleep_ti = exp_get_sessions_ti(exp_ID,'Sleep1','Sleep2');
rest_ti = exp.rest.ti;
flight_ti = [exp.flight.FE.start_ts; exp.flight.FE.end_ts]';
rest_ti(end,end) = 61550899128; % TODO: temp, need to fix the ts
sleep_ti(end,1) = 61604440361; % TODO: temp, need to fix the ts
epochs_ti = [sleep_ti; rest_ti; flight_ti];
epochs_clr = {repelem({sleep_clr},size(sleep_ti,1))
              repelem({rest_clr},size(rest_ti,1))
              repelem({flight_clr},size(flight_ti,1))};
epochs_clr = [epochs_clr{:}];
ylimits = get(gca,'YLim');
for ii_epoch = 1:size(epochs_ti,1)
    epoch_ti = epochs_ti(ii_epoch,:);
    clr = epochs_clr{ii_epoch};
    area(epoch_ti,ylimits([2 2]),'FaceColor',clr,'FaceAlpha',0.1,'EdgeColor','none','ShowBaseLine','off');
end
xticks([])
xlim(sleep_ti([1 end])+60e6.*[-1 1])
rescale_plot_data('x',[1e-6/60 0]); % change units to minutes
axis off
box off
% add scale bars
scale_bar_time_min = 2;
scale_bar_size_m = 50;
hax=gca;
xlimits = hax.XLim;
ylimits = hax.YLim;
x = xlimits(1) - 0.00001*range(xlimits);
y = ylimits(1) - 0.1*range(ylimits);
plot(x+[0 scale_bar_time_min], y+[0 0],'k-','LineWidth',2,'Clipping','off');
plot(x+[0 0], y+[0 scale_bar_size_m],'k-','LineWidth',2,'Clipping','off');
text(x+scale_bar_time_min/2, y, sprintf('%dmin',scale_bar_time_min) ,'FontSize',7,'HorizontalAlignment','center','VerticalAlignment','top');
text(x, y+scale_bar_size_m/2, sprintf('%dm',scale_bar_size_m) ,'FontSize',7,'HorizontalAlignment','center','VerticalAlignment','bottom','Rotation',90);
% text(-0.04,1.05, 'A', 'Units','normalized','FontWeight','bold');
hax.XLim = xlimits;
hax.YLim = ylimits;

%% load data for panels B and C
if ~exist('decode','var')
    exp_ID = 'b0184_d191208';
    exp = exp_load_data(exp_ID,'details','path','flight');
    dec_param_opt = 4;
    decode = decoding_load_data(exp_ID, 'flight', dec_param_opt);
    CM = load(fullfile('F:\sequences\decoded_figs\flight\conf_mat\',sprintf('%s_flight_decoding_opt_%d.mat',exp_ID,dec_param_opt)));
end

%% single flight examples
FEs = exp.flight.FE;
FEs([FEs.distance]<100)=[];
examples_flight_IX = [27 28];
for ii=1:2
    axes(panel_B(ii));
    cla
    hold on
    FE = FEs(examples_flight_IX(ii));
    ti = [FE.start_ts FE.end_ts];
    IX = get_data_in_ti(decode.time, ti);
    decode.MAP_pos(IX);
    plot([FE.ts], [FE.pos],'k','LineWidth',3);
    plot(decode.time(IX), decode.MAP_pos(IX),'r','LineWidth',.5);
    rescale_plot_data('x',[1e-6 ti(1)])
    xlabel('Time (s)')
    ylabel('Position (m)')
    xlim([0 20])
    ylim([0 150])
    xticks([0:10:20])
    yticks([0:50:150])
end
% axes(panel_B(1));
% text(-0.4,1.15, 'B', 'Units','normalized','FontWeight','bold');

%% single sessions confusion matrix example
x = CM.res_raw.pos_real;
y = CM.res_raw.pos_predict;
bin_centers = decode.pos;
[bin_edges,bin_size] = centers2edges(bin_centers);
N = histcounts2(x,y,bin_edges,bin_edges)'; % transpose so yaxis(rows)=predict and xaxis(cols)=real
N_norm_by_real = N ./ sum(N,1);
% N_norm_by_predict = N ./ sum(N,2);
axes(panel_C);
cla
hold on
imagesc(bin_centers,bin_centers,N_norm_by_real);
colormap(gca,flip(colormap('bone')));
hcb = colorbar('Location','eastoutside');
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
xlabel('Real position (m)')
ylabel('Predicted position (m)')
hax=gca;
hax.TickDir = 'out';
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = -1;
% text(-0.22,1.225, 'B', 'Units','normalized','FontWeight','bold');

%% error CDF plots (all sessions)
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
err_dist_edges = 0:0.1:100;
err_dist = nan(length(err_dist_edges)-1, length(exp_list));
pos_err_median_all = [];
for ii_exp = 1:length(exp_list)
    exp_ID = exp_list{ii_exp};
    filename = fullfile('F:\sequences\decoded_figs\flight\conf_mat\', ...
                        sprintf('%s_flight_decoding_opt_4.mat', exp_ID));
    CM = load(filename);
    err_dist(:,ii_exp) = histcounts(CM.res_raw.pos_err,'Normalization',err_dist_normalization,'BinEdges',err_dist_edges);
    pos_err_median_all(ii_exp) = CM.res.pos_err_median;
end
err_dist = [zeros(1,size(err_dist,2));err_dist]; % add zeros (cdf should start from 0)

%% plot CDFs
axes(panel_D);
cla
hold on
plot(err_dist_edges, err_dist,'LineWidth',0.3);
shadedErrorBar(err_dist_edges, err_dist', {@mean,@nansem},'lineprops',{'k','linewidth',3});
xlabel('Positional decoding error (m)')
switch err_dist_normalization
    case 'cdf'
        CDF_thr = 0.7;
        mean_CDF_error = err_dist_edges(find(mean(err_dist,2)>CDF_thr,1,'first'));
%         xline(mean_CDF_error,'-r');
%         yline(CDF_thr,'-r');
        grid on
        ylabel('Cumulative fraction')
    case 'pdf'
        ylabel('Probability density')
end
% text(-0.3,1.22, 'C', 'Units','normalized','FontWeight','bold');
hax=gca;
hax.TickLength(1) = [0.015];
hax.XScale = 'log';
% hax.XScale = 'linear';

%% some stats to report
fprintf('median error: %.2g+-%.2g (mean=-std)\n',mean(pos_err_median_all),std(pos_err_median_all));

for CDF_thr = [0.5 0.6 0.7 0.75 0.8 0.85 0.9 0.95]
    mean_CDF_error = err_dist_edges(find(mean(err_dist,2)>CDF_thr,1,'first'));
    fprintf('mean cumulative error at %d%%: %gm\n',CDF_thr*100,mean_CDF_error);
end


%% print/save the figure
fig_name_out = fullfile(res_dir, sprintf('%s_%s',fig_name_str,err_dist_normalization));

print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% exportgraphics(gcf,fullfile(res_dir,'testexport.pdf'),'BackgroundColor','none','ContentType','vector');

disp('figure was successfully saved to pdf/tiff/fig formats');

%%

