function plot_seqs_single_session(exp_ID)

%% plotting options
epoch_types = {'sleep','rest'};
clrs = [1 0 0; 0 0 1];

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Fig_replay_pop_hist';
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
panels_size = [3 3];
clear panels_hist
panels_hist(1) = axes('position', [2 20 panels_size]);
panels_hist(2) = axes('position', [7 20 panels_size]);
panels_hist(3) = axes('position', [12 20 panels_size]);
panels_hist(4) = axes('position', [17 20 panels_size]);
panel_over_represent = axes('position', [5 14 10 2.5]);
panel_TakeLandPastFuture = axes('position', [5 10 3 2.5]);
panel_seq_dynamics(1) = axes('position', [5 4 3 3]);
panel_seq_dynamics(2) = axes('position', [11 4 3 3]);

