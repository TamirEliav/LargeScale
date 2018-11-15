%%
clear all
close all
clc

%% params
smooth_p_turn = 0.01;
smooth_p_ybias = 0.01;
remove_bad_loc_opt = 2; % 1 - manual, 2 - from file
calib_day_opt = 1;
switch calib_day_opt
    case 1
        dir_out = 'L:\Analysis\Results\calib\20160320';
        load('L:\DATA\9892_Rocky\calib\20160320__tunnel_center_line_calib\bsp_pos_tag_2839.mat');
    case 2
    case 3
        dir_out = 'L:\Analysis\Results\calib\20170809';
        load('L:\DATA\0148_Boson\calib\20170809__tunnel_midline\client\bsp_pos_tag_2365.mat')
    case 4
    case 5
end
mkdir(dir_out)

%% remove bad localization points
figure
title('bad localization points removal')
hold on
hLine = plot(bsp_pos.pos(:,1), bsp_pos.pos(:,2), '.');
bad_loc_ix_filename = fullfile(dir_out, 'bad_loc_datapoints_ix');
switch remove_bad_loc_opt 
    case 1
        disp('chose data point to remove');
        while waitforbuttonpress==0; end
        brushedIdx = hLine.BrushData;
        disp([num2str(sum(brushedIdx)) ' bad data points to remove were chosen!']);
        
        save(bad_loc_ix_filename, 'brushedIdx')
    case 2
        load(bad_loc_ix_filename);
        hLine.BrushData = brushedIdx;
end

%% calibrate!
[curvexy, fit_res_all_seg, tunnel_length] = CALIB_fit_curve_200m(bsp_pos, smooth_p_turn, smooth_p_ybias);

% save calib results to mat file
clear calib_tunnel;
calib_tunnel.curvexy = curvexy;
calib_tunnel.tunnel_length = tunnel_length;
calib_tunnel.fitresults = fit_res_all_seg;
calib_filename = fullfile(dir_out, 'tunnel_calib');
save(calib_filename, 'calib_tunnel');

% save figures
figHandles = findobj('Type', 'figure');
figHandles = sort(figHandles,'descend');
for ii_fig = 1:length(figHandles)
    fig = figHandles(ii_fig);
    figure(fig);
    fig_filename = fullfile(dir_out, sprintf('fig_%d_%s', ii_fig, get(get(gca,'title'),'String'))); 
    saveas(gcf, fig_filename , 'fig')
    saveas(gcf, fig_filename , 'tif')
end

%%