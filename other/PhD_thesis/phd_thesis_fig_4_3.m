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
panels(2,1) = axes('position', [7.5 19 panels_size.*[2 1]]);
panels(3,1) = axes('position', [2 13 panels_size]);
panels(1,2) = axes('position', [4 20.5 2 2]);
% panels(2,2) = axes('position', [9.5 20.5 2 2]);
panels(3,2) = axes('position', [4 14.5 2 2]);
% panels(4,1) = axes('position', [2 13 panels_size]);
panels(4,1) = axes('position', [7.5 13 panels_size]);
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
FE_all = [];
FE_median_speed_all = [];
for ii_exp = 1:height(T)
    % load exp data
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details','flight');
    epoch_type = 'sleep';
%     epoch_type = 'rest';
    params_opt = 11;
    [events_session, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    events{ii_exp} = events_session;
    FE=exp.flight.FE;
    FE([FE.distance]<100)=[];
    FE_all = [FE_all FE];
    FE_median_speed_all(ii_exp) = median(abs([FE.vel]));
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
clr_replay = 0.*[1 1 1];
% clr_flight = [0.0863    0.7294    0.1294];
clr_flight = 0.5*[1 1 1];
lw_replay = 2;
lw_flight = 1;

%% Duration (replay)
axes(panels(1,1));
cla
hold on
x = [seqs.duration];
fprintf('\nReplay duration (s) stats:\n');
report_stats(x);
h=histogram(x);
% h.BinEdges = linspace(0,35,100);
h.DisplayStyle = 'stairs';
% h.FaceColor = clr_replay;
h.EdgeColor = clr_replay;
h.LineWidth = lw_replay;
m = median(x);
h=xline(m,'r',sprintf('%.2g',m));
h.LabelOrientation = 'horizontal';
h.FontSize = 8;
yticks([0 400])
xlabel('Replay duration (s)')
ylabel('Counts','Units','normalized','Position',[-0.1 0.5])
text(-0.27,1.1, 'A', 'Units','normalized','FontWeight','bold');

%% Duration (flight) - inset
axes(panels(1,2));
cla
hold on
x = [FE_all.duration];
fprintf('\nFlight duration (s) stats:\n');
report_stats(x);
h=histogram(x);
h.BinEdges = linspace(0,35,100);
h.DisplayStyle = 'stairs';
% h.FaceColor = clr_flight;
h.EdgeColor = clr_flight;
h.LineWidth = lw_flight;
m = median(x);
h=xline(m,'r',sprintf('%#.03g',m));
h.LabelOrientation = 'horizontal';
h.FontSize = 8;
yticks([0 400])
xlabel('Flight duration (s)','Units','normalized','Position',[0.5 -0.15])
ylabel('Counts','Units','normalized','Position',[-0.1 0.5])
hax=gca;
hax.XRuler.TickLabelGapOffset = -2;
hax.XAxis.TickLength(1) = 0.02;

%% Distance (replay)
axes(panels(2,1));
cla
hold on
x = [seqs.distance];
h=histogram(x);
% h.BinEdges = linspace(0,35,100);
h.DisplayStyle = 'stairs';
% h.FaceColor = clr_replay;
h.EdgeColor = clr_replay;
h.LineWidth = lw_replay;
h.Normalization='pdf';
m = median(x);
h=xline(m,'r',sprintf('%.2g',m));
h.LabelOrientation = 'horizontal';
h.FontSize = 8;
fprintf('\nReplay distance (m) stats:\n');
report_stats(x);

x = [FE_all.distance];
h=histogram(x);
h.BinWidth = 5;
% h.BinEdges = linspace(0,35,100);
% h.BinEdges = linspace(0,200,40);
h.DisplayStyle = 'stairs';
% h.FaceColor = clr_flight;
h.EdgeColor = clr_flight;
h.LineWidth = lw_flight;
h.Normalization='pdf';
m = median(x);
h=xline(m,'r',sprintf('%#.04g',m));
h.LabelOrientation = 'horizontal';
h.FontSize = 8;
h.LabelHorizontalAlignment = 'left';
fprintf('\nFlight distance (m) stats:\n');
report_stats(x);

yticks([0 0.1])
xlabel('Distance (m)')
ylabel('Probability density','Units','normalized','Position',[-0.03 0.5])

h1=plot(nan,nan,'Color',clr_replay,'LineWidth',lw_replay);
h2=plot(nan,nan,'Color',clr_flight,'LineWidth',lw_flight);
legend([h1 h2],'Replay','Flight','box','off')

text(-0.135,1.1, 'B', 'Units','normalized','FontWeight','bold');

%% Speed (replay)
axes(panels(3,1));
cla
hold on
x = [seqs.speed];
fprintf('\nReplay speed (m/s) stats:\n');
report_stats(x);
h=histogram(x);
% h.BinEdges = linspace(0,35,100);
h.DisplayStyle = 'stairs';
h.EdgeColor = clr_replay;
h.LineWidth = lw_replay;
m = median(x);
h=xline(m,'r',sprintf('%#.03g',m));
h.LabelOrientation = 'horizontal';
h.FontSize = 8;
yticks([0 500]);
xlabel('Replay speed (m/s)')
ylabel('Counts','Units','normalized','Position',[-0.1 0.5])
text(-0.27,1.1, 'C', 'Units','normalized','FontWeight','bold');

%% Speeed (flight) - inset
axes(panels(3,2));
cla
hold on
x = FE_median_speed_all; % per sessions
% x = arrayfun(@(x)(median([x.vel])), FE_all); % per individual flight epochs
fprintf('\nFlight speed (m/s) stats:\n');
report_stats(x);
h=histogram(x);
h.BinEdges = linspace(0,10,15);
h.DisplayStyle = 'stairs';
% h.FaceColor = clr_flight;
h.EdgeColor = clr_flight;
h.LineWidth = lw_flight;
m = median(x);
h=xline(m,'r',sprintf('%#.02g',m));
h.LabelOrientation = 'horizontal';
h.FontSize = 8;
yticks([0 40])
xlabel('Flight speed (m/s)','Units','normalized','Position',[0.5 -0.15])
ylabel('Counts','Units','normalized','Position',[-0.1 0.5])
hax=gca;
hax.XRuler.TickLabelGapOffset = -2;
hax.XAxis.TickLength(1) = 0.02;

%% Compression
axes(panels(4));
cla
hold on
x = [seqs.compression];
fprintf('\nReplay compression-ratio stats:\n');
report_stats(x);
h=histogram(x);
% h.BinEdges = linspace(0,35,100);
h.DisplayStyle = 'stairs';
h.EdgeColor = clr_replay;
h.LineWidth = lw_replay;
m = median(x);
h=xline(m,'r',sprintf('%#.03g',m));
h.LabelOrientation = 'horizontal';
h.FontSize = 8;
yticks([0 400]);
xlabel('Compression-ratio')
ylabel('Counts','Units','normalized','Position',[-0.1 0.5])
text(-0.27,1.1, 'D', 'Units','normalized','FontWeight','bold');

%% save fig(s)
if exist('bats_to_include','var')
    bats_str = ['bats_' char(strjoin(""+bats,'_'))];
else
    bats_str = 'bats_all';
end

fig_name_out = fullfile(res_dir, sprintf('%s_%s_%s',fig_name_str, epoch_type, bats_str));
print(fig, fig_name_out, '-dpdf', '-cmyk', '-painters');

disp('figure was successfully saved to pdf/tiff/fig formats');
diary off



%%
function report_stats(x)
fprintf('mean: %.2g\n',mean(x));
fprintf('median: %.2g\n',median(x));
fprintf('std: %.2g\n',std(x));
fprintf('CV: %.2g\n',std(x)/mean(x));
fprintf('IQR: %.2g-%.2g\n',prctile(x,[25 75]));
fprintf('range: %.2g-%.2g\n',[min(x) max(x)]);
end






