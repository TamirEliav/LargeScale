%% PhD thesis figure 4.3 - replay - population histograms

%%
close all
clear 
clc

%% options
% bats_to_include = [34 148 9861 2289];
% bats_to_include = [194 184 2382];
% bats_to_include = 34;
% bats_to_include = 148;
% bats_to_include = 9861;
% bats_to_include = 2289;
% bats_to_include = 194;
% bats_to_include = 184;
% bats_to_include = 2382;

%% define output files
res_dir = 'E:\Tamir\PhD\Thesis\resources\ch_4_seq';
mkdir(res_dir)
fig_name_str = 'Fig_4_3_replay_population';
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
panels_size = [4 4];
panels(1,1) = axes('position', [2 19 panels_size]);
panels(2,1) = axes('position', [7.5 19 panels_size]);
panels(3,1) = axes('position', [13 19 panels_size]);
% panels(1,2) = axes('position', [2 23 4 1]);
% panels(2,2) = axes('position', [7.5 23 4 1]);
% panels(3,2) = axes('position', [13 23 4 1]);
% panels(4) = axes('position', [2 13 panels_size]);
% panels(5) = axes('position', [7.5 13 panels_size]);
% panels(6) = axes('position', [13 13 panels_size]);
% panels(7) = axes('position', [2 7 panels_size]);
% panels(8) = axes('position', [7.5 7 panels_size]);
% panels(9) = axes('position', [13 7 panels_size]);

%% choose bats / sessions
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
clear exp_list
groupsummary(T,'bat_num')
if exist('bats_to_include','var')
    T = groupfilter(T,"bat_num",@(x)ismember(x,bats_to_include),'bat_num');
end
bats = unique(T.bat_num)

%% load data
events = {};
for ii_exp = 1:height(T)
    % load exp data
    exp_ID = T.exp_ID{ii_exp};
%     exp = exp_load_data(exp_ID,'details','path','ripples','MUA','PE');
    epoch_type = 'sleep';
    params_opt = 11;
    [events_session, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    events{ii_exp} = events_session;
end
T.nEvents = cellfun(@length,events)'; % note this is without filtering sequence by features!
sortrows( groupsummary(T,'bat_num',["median","mean","max","sum"],"nEvents"),"sum_nEvents", 'descend')

%% pool data
events = [events{:}];

%% apply inclusion criteria 
seqs = [events.seq_model];
% seqs([seqs.compression]<2)=[];
seqs([seqs.score]<0.5)=[];
seqs([seqs.distance]<3)=[];

%%
clr = 0.5*[1 1 1];

%% Duration
axes(panels(1,1));
cla
hold on
x = [seqs.duration];
histogram(x,'FaceColor',clr);
xline(median(x),'r',sprintf('%.2g',median(x)))
xlabel('Duration (s)')
ylabel('Counts')
text(-0.27,1.1, 'A', 'Units','normalized','FontWeight','bold');
fprintf('\nDuration (s) stats:\n');
report_stats(x);
% axes(panels(1,2));
% cla
% hold on
% h=boxplot(x,'Orientation','horizontal','BoxStyle','filled');
% axis off

%% Distance
axes(panels(2,1));
cla
hold on
x = [seqs.distance];
histogram(x,'FaceColor',clr);
xline(median(x),'r',sprintf('%.2g',median(x)))
xlabel('Distance (m)')
ylabel('Counts')
text(-0.27,1.1, 'B', 'Units','normalized','FontWeight','bold');
fprintf('\nDistance (m) stats:\n');
report_stats(x);
% axes(panels(2,2));
% cla
% hold on
% h=boxplot(x,'Orientation','horizontal','BoxStyle','filled');
% axis off

%% Compression
axes(panels(3,1));
cla
hold on
x = [seqs.compression];
histogram(x,'FaceColor',clr);
xline(median(x),'r',sprintf('%.2g',median(x)))
xlabel('Compression')
ylabel('Counts')
text(-0.27,1.1, 'C', 'Units','normalized','FontWeight','bold');
fprintf('\nCompression stats:\n');
report_stats(x);
% axes(panels(3,2));
% cla
% hold on
% h=boxplot(x,'Orientation','horizontal','BoxStyle','outline','PlotStyle','traditional','MedianStyle','line','Notch','on');
% axis off

%% link axes
% linkaxes(panels(1,:),'x');
% linkaxes(panels(2,:),'x');
% linkaxes(panels(3,:),'x');

%% save fig(s)
if exist('bats_to_include','var')
    bats_str = ['_bats_' char(strjoin(""+bats,'_'))];
else
    bats_str = '_bats_all';
end

fig_name_out = fullfile(res_dir, [fig_name_str bats_str]);
print(fig, fig_name_out, '-dpdf', '-cmyk', '-painters');

disp('figure was successfully saved to pdf/tiff/fig formats');



%%
function report_stats(x)
fprintf('mean: %.2g\n',mean(x));
fprintf('median: %.2g\n',median(x));
fprintf('std: %.2g\n',std(x));
fprintf('CV: %.2g\n',std(x)/mean(x));
fprintf('IQR: %.2g-%.2g\n',prctile(x,[25 75]));
fprintf('range: %.2g-%.2g\n',[min(x) max(x)]);
end






