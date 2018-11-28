%%
clear
clc

%% load data
pos1 = load('L:\DATA\0148_Boson\calib\20170809__tunnel_midline\client\bsp_pos_tag_2365.mat');
pos2 = load('L:\DATA\9892_Rocky\calib\20160320__tunnel_center_line_calib\bsp_pos_tag_2839.mat');
pos3 = load('L:\DATA\9861_Somo\calib\20180606_calib_middle_line_top_14_anchors_active\bsp\client\bsp_pos_tag_708.mat')

%%
figure('units','normalized','outerposition',[0 0 1 1])
hold on
plot(pos1.bsp_pos.pos(:,1), pos1.bsp_pos.pos(:,2), '.b')
plot(pos2.bsp_pos.pos(:,1), pos2.bsp_pos.pos(:,2), '.r')
plot(pos3.bsp_pos.pos(:,1), pos3.bsp_pos.pos(:,2), '.m')
legend({'Boson - 09/08/2017';'Rocky - 20/03/2016';'Somo - 06/06/2018'})
title('comparing calibration biases over 1 year')
xlim([1360 1500])
ylim([2499.5 2501])
% fileout = 'L:\Analysis\Results\calib\compare_calib_over_year';
% saveas(gcf, fileout, 'tif')







%%
