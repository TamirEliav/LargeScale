%% Large Scale - Fig. S8 - model gamma fit

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'Fig_S8_field_size_hist_gamma_fit';
fig_caption_str = 'Theoretical analysis - fields size distribution gamma fit';
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
panel_A(2) = axes('position', [ 8.5 21 4 3]);

%% load data
load('L:\paper_figures\pop_dist_fields_size.mat');

%% fit gamma / log-normal to field size distribution
gamma_phat = gamfit(fields_size);
logn_phat = lognfit(fields_size);

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
Y = gampdf(X, gamma_phat(1), gamma_phat(2));
plot(X,Y,'-r', 'LineWidth', 1.2);
Y = lognpdf(X, logn_phat(1), logn_phat(2));
plot(X,Y,'-g', 'LineWidth', 1.2);

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

text(0.4,1.0, 'Gamma distribution',                   'Units','normalized','FontSize',7, 'HorizontalAlignment','left');
text(0.4,0.9, ['\alpha = ' num2str(gamma_phat(1),3)], 'Units','normalized','FontSize',7, 'HorizontalAlignment','left');
text(0.4,0.8, ['\beta= '   num2str(gamma_phat(2),3)], 'Units','normalized','FontSize',7, 'HorizontalAlignment','left');

text(0.4,0.6, 'Log-Normal distribution',              'Units','normalized','FontSize',7, 'HorizontalAlignment','left');
text(0.4,0.5, ['\mu= ' num2str(logn_phat(1),3)], 'Units','normalized','FontSize',7, 'HorizontalAlignment','left');
text(0.4,0.4, ['\sigma= '   num2str(logn_phat(2),3)], 'Units','normalized','FontSize',7, 'HorizontalAlignment','left');

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
Y = lognpdf(X, logn_phat(1), logn_phat(2));
plot(X,Y,'-g', 'LineWidth', 1.2);

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


%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

