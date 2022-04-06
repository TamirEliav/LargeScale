%% PhD thesis figure 5.2 - discussion: dorso-ventral illustration

%%
close all
clear 
clc

%% define output files
res_dir = 'E:\Tamir\PhD\Thesis\resources\ch_5_discussion';
mkdir(res_dir)
fig_name_str = 'Fig_5_2_dorso_ventral';
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
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

% create panels
panels(1) = axes('position', [3 18 4 4]);
panels(2) = axes('position', [8 18 4 4]);

w = 0.4;

axes(panels(1))
hold on
patch([0.1 0.1 0.9 0.9],[0.1 0.1+w 0.9 0.9-w], 'k');
hax=gca;
hax.YTick = [];
hax.XTick = [0.2 0.8];
hax.XTickLabel = {'Dorsal','Ventral'};
hax.XAxis.FontSize = 10;
hax.TickLength = [0 0];
xlim([0 1]);
ylim([0 1]);
ylabel('Fields size','Units','normalized','Position',[-0.15 0.5],'FontSize',12);

axes(panels(2))
hold on
patch([0.1 0.1 0.9 0.9],[0.1 0.1+w 0.1+w 0.1+w], 'k');
hax=gca;
hax.YTick = [];
hax.XTick = [0.2 0.8];
hax.XTickLabel = {'Dorsal','Ventral'};
hax.XAxis.FontSize = 10;
hax.TickLength = [0 0];
xlim([0 1]);
ylim([0 1]);
% ylabel('Fields size')

%% save fig(s)
fig_name_out = fullfile(res_dir, [fig_name_str ]);
print(fig, fig_name_out, '-dpdf', '-cmyk', '-painters');

disp('figure was successfully saved to pdf/tiff/fig formats');
diary off

%%
