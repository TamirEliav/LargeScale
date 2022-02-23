%% PhD thesis figure 4.2 - replay examples (legend)

%%
close all
clear 
clc

%% define output files
res_dir = 'E:\Tamir\PhD\Thesis\resources\ch_4_seq';
mkdir(res_dir)
fig_name_str = 'Fig_4_2_replay_examples_legend';
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
panels_size = [1 1];
panels(1) = axes('position', [10 10 panels_size]);

%% create legend
axes(panels(1));
hold on
axis off
plot(nan(6,6),'-','LineWidth',2);
hl=legend( "Continuous (direction 1)","Stationary (direction 1)","Fragmented (direction 1)", ...
           "Continuous (direction 2)","Stationary (direction 2)","Fragmented (direction 2)", 'NumColumns',2);
hl.Box='off';

arrow_len = 0.03;

ar = annotation('arrow');
ar.Color = 'k';
ar.Position = [0.13 0.402 0 -arrow_len];
ar.HeadStyle = 'cback1';
ar.HeadWidth = 5;
ar.HeadLength = 5;

ar = annotation('arrow');
ar.Color = 'k';
ar.Position = [0.51 0.372 0 arrow_len];
ar.HeadStyle = 'cback1';
ar.HeadWidth = 5;
ar.HeadLength = 5;


%% save fig
fig_name_out = fullfile(res_dir, [fig_name_str]);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
disp('figure was successfully saved to pdf/tiff/fig formats');









%%
