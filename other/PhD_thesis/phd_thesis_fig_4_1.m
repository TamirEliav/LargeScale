%% PhD thesis figure 4.1 - deocding accuracy (flight)

%%
close all
clear 
clc

%% plotting options

%% define output files
res_dir = 'E:\Tamir\PhD\Thesis\resources\ch_4_seq';
mkdir(res_dir)
fig_name_str = 'Fig_4_1';
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
figure ;
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
panel_A(1) = axes('units', 'centimeters', 'position', [3  21 2 2]);
panel_A(2) = axes('units', 'centimeters', 'position', [3  18 2 2]);
panel_B(1) = axes('units', 'centimeters', 'position', [6.5 17.5 6 6]);
panel_C(1) = axes('units', 'centimeters', 'position', [13.64 18.76 4 4]);
panel_C(2) = axes('units', 'centimeters', 'position', [15.64 19.5 2 2]);

%% load data for panels A and B
exp_ID = 'b0184_d191208';
exp = exp_load_data(exp_ID,'details','path','flight');
decode = decoding_load_data(exp_ID, 'flight', 4);
CM = load('F:\sequences\decoded_figs\flight\conf_mat\b0184_d191208_flight_decoding_opt_4.mat');

%% single flight examples
FEs = exp.flight.FE;
FEs([FEs.distance]<100)=[];
examples_flight_IX = [27 28];
for ii=1:2
    axes(panel_A(ii));
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
end
axes(panel_A(1));
text(-0.6,1.3, 'A', 'Units','normalized','FontWeight','bold');

%% single sessions confusion matrix example
x = CM.res_raw.pos_real;
y = CM.res_raw.pos_predict;
bin_centers = decode.pos;
[bin_edges,bin_size] = centers2edges(bin_centers);
N = histcounts2(x,y,bin_edges,bin_edges)'; % transpose so yaxis(rows)=predict and xaxis(cols)=real
N_norm_by_real = N ./ sum(N,1);
% N_norm_by_predict = N ./ sum(N,2);
axes(panel_B);
cla
hold on
imagesc(bin_centers,bin_centers,N_norm_by_real);
colormap(gca,flip(colormap('bone')));
hcb = colorbar('Location','eastoutside');
hcb.Limits = [0 1];
hcb.Label.String = 'Probability';
hcb.Label.Rotation = -90;
hcb.Label.Position = [2.1 0.5 0];
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
text(-0.22,1.225, 'B', 'Units','normalized','FontWeight','bold');

%% error CDF plots (all sessions)
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
CDFs = [];
CDF_edges = 0:0.1:100;
pos_err_median_all = [];
for ii_exp = 1:length(exp_list)
    exp_ID = exp_list{ii_exp};
    filename = fullfile('F:\sequences\decoded_figs\flight\conf_mat\', ...
                        sprintf('%s_flight_decoding_opt_4.mat', exp_ID));
    CM = load(filename);
    CDFs(:,ii_exp) = histcounts(CM.res_raw.pos_err,'Normalization','cdf','BinEdges',CDF_edges);
    pos_err_median_all(ii_exp) = CM.res.pos_err_median;
end
CDFs = [zeros(1,size(CDFs,2));CDFs]; % add zeros 

%% some stats to report
fprintf('median error: %.2g+-%.2g (mean=-std)\n',mean(pos_err_median_all),std(pos_err_median_all));

CDF_thr = 0.5;
mean_CDF_error = CDF_edges(find(mean(CDFs,2)>CDF_thr,1,'first'));
fprintf('mean cumulative error at %d%%: %gm\n',CDF_thr*100,mean_CDF_error);
CDF_thr = 0.7;
mean_CDF_error = CDF_edges(find(mean(CDFs,2)>CDF_thr,1,'first'));
fprintf('mean cumulative error at %d%%: %gm\n',CDF_thr*100,mean_CDF_error);
CDF_thr = 0.75;
mean_CDF_error = CDF_edges(find(mean(CDFs,2)>CDF_thr,1,'first'));
fprintf('mean cumulative error at %d%%: %gm\n',CDF_thr*100,mean_CDF_error);
CDF_thr = 0.8;
mean_CDF_error = CDF_edges(find(mean(CDFs,2)>CDF_thr,1,'first'));
fprintf('mean cumulative error at %d%%: %gm\n',CDF_thr*100,mean_CDF_error);

%% plot CDFs
axes(panel_C(1));
cla
hold on
plot(CDF_edges, CDFs,'LineWidth',0.3);
shadedErrorBar(CDF_edges, CDFs', {@mean,@nansem},'lineprops',{'k','linewidth',3});
xlabel('Position error (m)')
ylabel('Cumulative fraction')
text(-0.3,1.22, 'C', 'Units','normalized','FontWeight','bold');
hax=gca;
hax.TickLength(1) = [0.015];

%% plot CDFs (zoom-in inset)
axes(panel_C(2));
cla
hold on
plot(CDF_edges, CDFs,'LineWidth',0.3);
shadedErrorBar(CDF_edges, CDFs', {@mean,@nansem},'lineprops',{'k','linewidth',3});
CDF_thr = 0.5;
m = CDF_edges(find(mean(CDFs,2)>CDF_thr,1,'first'));
plot([0 m],CDF_thr*[1 1],'r-')
plot([m m],[0 CDF_thr],'r-')
xlim([0 5])
ylim([0 1])
xticks([0 m 5])
yticks([0 0.5 1])
hax=gca;
hax.TickLength(1) = [0.03];

%% save fig
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
disp('figure was successfully saved to pdf/tiff/fig formats');
