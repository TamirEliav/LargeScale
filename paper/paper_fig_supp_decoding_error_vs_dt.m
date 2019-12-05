%% Large Scale - Fig. S10 - model gamma fit + decoder dependency on dt (integration time)

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'Fig_S10';
fig_caption_str = 'Theoretical analysis - error vs. integration time (dt)';
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
disp([fig_name_str ': ' fig_caption_str]);
disp('======================================================');
disp('');

%% create figure
% =========================================================================
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
pause(0.2); % workaround to solve matlab automatically changing the axes positions...

% create panels
panel_A(1) = axes('position', [ 3 21 4 3]);
panel_A(2) = axes('position', [ 8 21 4 3]);
panel_B    = axes('position', [ 3 14 7 5]);
panel_legend = axes('position', [10.6 16.8 0.5 2]);

%% load data
load('L:\paper_figures\pop_dist_fields_size.mat');

%% fit gamma to field size distribution
gamma_phat = gamfit(fields_size);

%% panel A - field size distribution gamma fit
axes(panel_A(1));
cla
hold on
text(-0.3,1.1, 'A', 'Units','normalized','FontWeight','bold');

h = histogram(fields_size);
h.FaceColor = 0.5*[1 1 1];
h.BinWidth = 1;
h.Normalization = 'pdf';

X = linspace(1,50,100);
Y = gampdf(X,gamma_phat(1), gamma_phat(2));
plot(X,Y,'-r', 'LineWidth', 1.2);

ha= gca;
ha.XLim = [0 35];
% ha.YLim = [0.8 310];
% ha.XTick = [0:5:35];
% ha.YTick = [1 10 100];
% ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel('Field size (m)')
ylabel('PDF','Units','normalized','Position',[-0.15 0.5])
ha.XScale = 'linear';
ha.YScale = 'linear';

text(0.4,0.9, 'Gamma distribution',                   'Units','normalized','FontSize',7, 'HorizontalAlignment','left');
text(0.4,0.8, ['\alpha = ' num2str(gamma_phat(1),3)], 'Units','normalized','FontSize',7, 'HorizontalAlignment','left');
text(0.4,0.7, ['\beta= '   num2str(gamma_phat(2),3)], 'Units','normalized','FontSize',7, 'HorizontalAlignment','left');

%% panel A - field size distribution gamma fit (log scale)
axes(panel_A(2));
cla
hold on

h = histogram(fields_size);
h.FaceColor = 0.5*[1 1 1];
h.BinWidth = 1;
h.Normalization = 'pdf';

X = linspace(1,50,100);
Y = gampdf(X,gamma_phat(1), gamma_phat(2));
plot(X,Y,'-r', 'LineWidth', 1.2);

ha= gca;
ha.XLim = [0 35];
% ha.YLim = [0.8 310];
% ha.XTick = [0:5:35];
% ha.YTick = [1 10 100];
% ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel('Field size (m)')
ylabel('PDF','Units','normalized','Position',[-0.15 0.5])
ha.XScale = 'log';
ha.YScale = 'linear';

%% load decoding data
load('L:\Theory\Yonatan_code_data\Multiscale_PV_ML_decoder_5_Summary_1.mat');
nL = length(L);
ds = 20 ;
% decoder schemes colors
clr  = [1.0 0.0 0.5   ; ... % 1
        0.0 1.0 0.5   ; ... % 2
        1.0 0.5 0.0   ; ... % 3
        0.25 0.75 1.0 ; ... % 4
        0.5 0.0 1.0   ; ... % 5
        0.5 0.9 1.0   ; ... % 6
        0.85 0.4 1.0  ; ... % 7
        1.0 0.75 0.25 ; ... % 8
        0   0    0    ; ... % 9
        0.7 0.7 0.7 ] ;     % 10

%% panel B - decoding error vs. dt
axes(panel_B);
cla
hold on
text(-0.2,1.1, 'B', 'Units','normalized','FontWeight','bold');

jN = find(ismember(N, [50]));
jL = find(ismember(L, [20000]));
plot(dt.*1e3, squeeze(meMLA(jN,:,jL)), 'Color', clr(1,:), 'LineWidth', 2);
plot(dt.*1e3, squeeze(meMLB(jN,:,jL)), 'Color', clr(3,:), 'LineWidth', 2);
plot(dt.*1e3, squeeze(meMLI(jN,:,jL)), 'Color', clr(2,:), 'LineWidth', 2);
plot(dt.*1e3, squeeze(meMLF(jN,:,jL)), 'Color', clr(6,:), 'LineWidth', 2);
plot(dt.*1e3, squeeze(meMLG(jN,:,jL)), 'Color', clr(7,:), 'LineWidth', 2);

ha= gca;
ha.XLim = [20 500];
% ha.YLim = [0.8 310];
ha.XTick = [20 100:100:500];
% ha.YTick = [1 10 100];
% ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
ha.TickDir='out';
ha.TickLength = [0.02 0.02];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel('Integration window (ms)')
ylabel('Mean decoding error (m)','Units','normalized','Position',[-0.12 0.5]);
ha.XScale = 'linear';
ha.YScale = 'log';

%% legend panel
axes(panel_legend);
cla
hold on
plot([0 1], 5*[1 1], 'Color',clr(1,:),'LineWidth',2) ; hold on ; 
plot([0 1], 4*[1 1], 'Color',clr(3,:),'LineWidth',2) ; hold on ; 
plot([0 1], 3*[1 1], 'Color',clr(2,:),'LineWidth',2) ; hold on ; 
plot([0 1], 2*[1 1], 'Color',clr(6,:),'LineWidth',2) ; hold on ; 
plot([0 1], 1*[1 1], 'Color',clr(7,:),'LineWidth',2) ; hold on ; 
text(1.2,5,'1: Single small fields', 'FontSize', 10);
text(1.2,4,'2: Single large fields', 'FontSize', 10);
text(1.2,3,'3: multiple small fields (Rich et al. 2014)', 'FontSize', 10);
text(1.2,2,'4: multi-scale (population)', 'FontSize', 10);
text(1.2,1,'5: multi-scale (single-cell)', 'FontSize', 10);
axis off
xlim([0 1])
ylim([1 5])


%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

