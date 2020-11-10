%% Large Scale - Fig. S6 - model gamma/log-normal fits

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S7';
fig_caption_str = 'Theoretical analysis - fields size distribution gamma/log-normal fit';
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
panel_A(1) = axes('position', [ 3.0 21 4 3]);
panel_A(2) = axes('position', [ 9.5 21 4 3]);

%% load data
load('L:\paper_figures\pop_dist_fields_size.mat');

%% fit gamma / log-normal to field size distribution
gamma_phat = gamfit(fields_size);
logn_phat = lognfit(fields_size);

fprintf('Log-normal fit: mu=%.5f, sigma=%.5f\n', logn_phat(:) );
fprintf('Gamma fit: alpha=%.5f, beta=%.5f\n', gamma_phat(:) )

%% panel A - field size distribution gamma fit
axes(panel_A(1));
cla
hold on
% text(-0.3,1.1, 'A', 'Units','normalized','FontWeight','bold');

h = histogram(fields_size);
h.FaceColor = 0.7*[1 1 1];
h.BinWidth = 1;
h.Normalization = 'pdf';

X = linspace(1,50,100);
Y = gampdf(X, gamma_phat(1), gamma_phat(2));
plot(X,Y,'-g', 'LineWidth', 1.5);
Y = lognpdf(X, logn_phat(1), logn_phat(2));
plot(X,Y,'-r', 'LineWidth', 1.5);

ha= gca;
ha.XLim = [0 35];
% ha.YLim = [0.8 310];
% ha.XTick = [0:5:35];
% ha.YTick = [1 10 100];
% ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.2;
xlabel('Field size (m)','Units','normalized','Position',[0.5 -0.16])
ylabel('Probability density function','Units','normalized','Position',[-0.22 0.5])
ha.XScale = 'linear';
ha.YScale = 'linear';

plot([12 16],[0.20 0.20], 'r', 'LineWidth', 2)
plot([12 16],[0.12 0.12], 'g', 'LineWidth', 2)
text(0.5,1.0, 'Log-Normal distribution',              'Units','normalized','FontSize',7, 'HorizontalAlignment','left');
text(0.55,0.9, ['\mu = ' num2str(logn_phat(1),3)], 'Units','normalized','FontSize',7, 'HorizontalAlignment','left');
text(0.55,0.8, ['\sigma = '   num2str(logn_phat(2),3)], 'Units','normalized','FontSize',7, 'HorizontalAlignment','left');
text(0.5,0.6, 'Gamma distribution',                   'Units','normalized','FontSize',7, 'HorizontalAlignment','left');
text(0.55,0.5, ['\alpha = ' num2str(gamma_phat(1),3)], 'Units','normalized','FontSize',7, 'HorizontalAlignment','left');
text(0.55,0.4, ['\beta = '   num2str(gamma_phat(2),3)], 'Units','normalized','FontSize',7, 'HorizontalAlignment','left');


%% panel A - field size distribution gamma fit (log scale)
axes(panel_A(2));
cla
hold on

h = histogram(fields_size);
h.FaceColor = 0.7*[1 1 1];
h.BinWidth = 1;
h.Normalization = 'pdf';

X = linspace(1,50,100);
Y = gampdf(X, gamma_phat(1), gamma_phat(2));
plot(X,Y,'-g', 'LineWidth', 1.5);
Y = lognpdf(X, logn_phat(1), logn_phat(2));
plot(X,Y,'-r', 'LineWidth', 1.5);

ha= gca;
ha.XLim = [0 35];
% ha.YLim = [0.8 310];
% ha.XTick = [0:5:35];
ha.XTick = [1 10];
ha.XTickLabel = {'10^{ 0}';'10^{ 1}'};
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.2;
xlabel('Field size (m)','Units','normalized','Position',[0.5 -0.16])
ylabel('Probability density function','Units','normalized','Position',[-0.22 0.5])
ha.XScale = 'log';
ha.YScale = 'linear';

%% Goodnes-of-fit
pd_logn = makedist('Lognormal','mu',logn_phat(1),'sigma',logn_phat(2));
pd_gamma = makedist('Gamma','a',gamma_phat(1),'b',gamma_phat(2));
[~,pval_logn, ksstat_logn ] = kstest(fields_size,'CDF',pd_logn);
[~,pval_gamma,ksstat_gamma] = kstest(fields_size,'CDF',pd_gamma);
fprintf('logn:  P=%.2g,\t\tks2stat=%.2g\n',pval_logn,  ksstat_logn);
fprintf('gamma: P=%.2g,\tks2stat=%.2g\n',pval_gamma ,ksstat_gamma);
% [cdf_data,X] = ecdf(fields_size);
% figure
% hold on
% plot(X,cdf_data,'k');
% plot(X,pd_gamma.cdf(X),'b');
% plot(X,pd_logn.cdf(X),'r');
% legend('data','gamma','logn')

% we decided to report both pval and KS_stat in the fig legend

%% Goodness-of-fit - ignore!
% % % % [cdf_data,X] = ecdf(fields_size);
% % % % cdf_gamma = gamcdf(X, gamma_phat(1), gamma_phat(2));
% % % % cdf_logn = logncdf(X, logn_phat(1), logn_phat(2));
% % % % n = size(fields_size);
% % % % n = [1 1e4];
% % % % fields_size_gamma = gamrnd(gamma_phat(1), gamma_phat(2), n);
% % % % fields_size_logn = lognrnd(logn_phat(1), logn_phat(2), n);
% % % % [h1,p1,ks2stat1] = kstest2(fields_size, fields_size_logn);
% % % % [h2,p2,ks2stat2] = kstest2(fields_size, fields_size_gamma);
% % % % figure
% % % % hold on
% % % % plot(X,cdf_data,'k');
% % % % plot(X,cdf_gamma,'b');
% % % % plot(X,cdf_logn,'r');
% % % % fprintf('logn: P=%.2g, ks2stat=%.2g, manual ks2stat=%.2g, \n',p1,ks2stat1, max(abs(cdf_data-cdf_logn)));
% % % % fprintf('gamma: P=%.2g, ks2stat=%.2g, manual ks2stat=%.2g, \n',p2,ks2stat2, max(abs(cdf_data-cdf_gamma)));
% % % % legend('data','gamma','logn')



%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

