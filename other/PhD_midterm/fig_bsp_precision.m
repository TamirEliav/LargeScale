%%
clear
clc

%% load data
load('L:\BeSpoon\testing\test_20180530__YOM_KEF_200m_static+dynamic+discretization\30-05-2018__calib_test_dynamic+non-jitter_jitter_with_kalman\data\outside_perpendicular_error.mat')

%% plot
fig_size = [7 6];
figure('Units','centimeters','Position',[5 5 fig_size], 'PaperPosition', [0 0 fig_size])
err = perpendicular_error.*100;
h=histogram(err,'Normalization','probability');
hax = gca;
hax.XTick = [-.4: .2: .4].*100;
hax.XLim= [-.4 .4].*100;
hax.YTick = [0 .1 .2];
hax.YLim= [0 .2];
hax.FontSize = 10;
xlabel('X position error (cm)')
ylabel('Prob.')

% text(0.6,0.95,sprintf('mean = %.1f cm',mean(err)) , 'Units','normalized','FontSize',10)
% text(0.6,0.85,sprintf('median = %.1f cm',median(err)) , 'Units','normalized','FontSize',10)
text(0.5,0.9,sprintf('std = %.1f cm',std(err)) , 'Units','normalized','FontSize',10,'HorizontalAlignment','center')

box off

%% save figure
dir_out = 'L:\Analysis\Results\midterm';
mkdir(dir_out);
filename = fullfile(dir_out, 'fig_bsp_precision')
saveas(gcf,filename,'fig')
saveas(gcf,filename,'tif')
saveas(gcf,filename,'pdf')
