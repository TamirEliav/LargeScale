%% besoon dynamic test

%%
clear all
clc

%%
% test_opt = 'inside';
test_opt = 'outside';

%%
main_dir = 'L:\BeSpoon\testing\test_20180530__YOM_KEF_200m_static+dynamic+discretization\30-05-2018__calib_test_dynamic+non-jitter_jitter_with_kalman';
nlg_dir = fullfile(main_dir, 'nlg');
nlx_dir = fullfile(main_dir, 'nlx');
bsp_dir = fullfile(main_dir, 'bsp');
sync_dir = fullfile(main_dir, 'sync');
results_dir = fullfile(main_dir, 'results');
mkdir(results_dir)
% Nlg2Nlx(main_dir);
% PRE_sync_bsp_to_nlg(bsp_dir, nlx_dir, sync_dir);

%% load bsp data
load('L:\BeSpoon\testing\test_20180530__YOM_KEF_200m_static+dynamic+discretization\30-05-2018__calib_test_dynamic+non-jitter_jitter_with_kalman\bsp\bsp_pos_tag_3570.mat');
load(['L:\BeSpoon\testing\test_20180530__YOM_KEF_200m_static+dynamic+discretization\30-05-2018__calib_test_dynamic+non-jitter_jitter_with_kalman\data\' test_opt '_tunnel_throws.mat'])

switch test_opt
    case 'inside'
        balls_meas_ts = [
        50729311 50767588;	% ball 1
        50628470 50699245	% ball 2
        ].*1e3;
    case 'outside'
        balls_meas_ts = [
        55127709 55162085;	% ball 1
        55185664 55218836% ball 2
        ].*1e3;
end
for ii_ball = 1:2
    IX = find(bsp_pos.ts_nlg_usec > balls_meas_ts(ii_ball,1) & bsp_pos.ts_nlg_usec < balls_meas_ts(ii_ball,2));
    balls_pos(ii_ball,:) = mean(bsp_pos.pos(IX,1:2));
end

%%
% figure
% hold on
% plot(bsp_pos.ts_nlg_usec, bsp_pos.pos(:,1:2),'.')
% plot(balls_meas_ts(:), 2500,'r*')
% plot(throws_comments_ts, 2470, 'db')
% IX = find(diff(throws_pos_ts(:,1))*1e-6 > 2);
% throws_num = zeros(1,size(throws_pos_ts,1));
% throws_num(1) = 1;
% throws_num(IX+1) = 1;
% throws_num = cumsum(throws_num);
% plot(throws_pos_ts(IX,1), throws_pos_ts(IX,2),'mo')
% plot(throws_pos_ts(:,1), throws_num +2450,'g.')

%%
% clear throws
% throws.ts = throws_pos_ts(:,1);
% throws.pos = interp1(bsp_pos.ts_nlg_usec, bsp_pos.pos, throws_pos_ts(:,1));
% throws.comments_ts = throws_comments_ts;
% throws.throws_num = throws_num;
% throws.throws_cat = interp1(1:length(throws.cat), throws.cat, throws.throws_num)
% save(['L:\BeSpoon\testing\test_20180530__YOM_KEF_200m_static+dynamic+discretization\30-05-2018__calib_test_dynamic+non-jitter_jitter_with_kalman\data\' test_opt  '_tunnel_throws'],'throws')

%%
[xy,distance,t_a] = distance2curve(balls_pos, throws.pos(:,1:2),'linear');
side = side_of_line(balls_pos(1,:),balls_pos(2,:), throws.pos(:,1:2));
perpendicular_error = side.*distance;
save(['L:\BeSpoon\testing\test_20180530__YOM_KEF_200m_static+dynamic+discretization\30-05-2018__calib_test_dynamic+non-jitter_jitter_with_kalman\data\' test_opt  '_perpendicular_error'],'perpendicular_error')

%%
figure
h=histogram(perpendicular_error)
h.Normalization = 'probability';
text(0.7,0.9, sprintf('mean=%.2gm',mean(perpendicular_error)), 'unit', 'norm')
text(0.7,0.8, sprintf('median=%.2gm',median(perpendicular_error)), 'unit', 'norm')
text(0.7,0.7, sprintf('std=%.2gm',std(perpendicular_error)), 'unit', 'norm')
xlabel('Error in perpendicular direction (m)')
ylabel('prob.')
switch test_opt
    case 'inside'
        title({'Dynamic test inside tunnel';'Linear movement in X direction';'(testing perpendicular errors in Y!)'});
    case 'outside'
        title({'Dynamic test outside tunnel';'Linear movement in Y direction';'(testing perpendicular errors in X!)'});
end
filename = sprintf('dynamic_test_%s', test_opt);
fileout = fullfile(results_dir, filename);
saveas(gcf, fileout,'tif')
saveas(gcf, fileout,'fig')

%% per category histograms
% % % % % % % % figure
% % % % % % % % p = panel();
% % % % % % % % p.pack('h',3);
% % % % % % % % p.margin = 20;
% % % % % % % % % accumarray(throws.throws_cat', perpendicular_error', [], @mean)
% % % % % % % % % accumarray(throws.throws_cat', perpendicular_error', [], @median)
% % % % % % % % % accumarray(throws.throws_cat', perpendicular_error', [], @std)
% % % % % % % % cat_strs = {'Excellent';'Good';'Bad'};
% % % % % % % % for ii_cat = 1:3
% % % % % % % %     p(ii_cat).select()
% % % % % % % %     err = perpendicular_error(throws.throws_cat==ii_cat);
% % % % % % % % %     subplot(1,3,ii_cat)
% % % % % % % %     h=histogram(err,25);
% % % % % % % %     h.Normalization = 'probability';
% % % % % % % %     text(0,0.95, sprintf('mean=%.2gcm',100*mean(err)), 'unit', 'norm')
% % % % % % % %     text(0,0.9, sprintf('median=%.2gcm',100*median(err)), 'unit', 'norm')
% % % % % % % %     text(0,0.85, sprintf('std=%.2gcm',100*std(err)), 'unit', 'norm')
% % % % % % % % end
% % % % % % % % p.xlabel('Error in perpendicular direction (m)')
% % % % % % % % p.ylabel('prob.')
% % % % % % % % 
% % % % % % % % switch test_opt
% % % % % % % %     case 'inside'
% % % % % % % %         suptitle({'Dynamic test inside tunnel';'Linear movement in X direction (testing perpendicular errors in Y!)'});
% % % % % % % %     case 'outside'
% % % % % % % %         suptitle({'Dynamic test outside tunnel';'Linear movement in Y direction (testing perpendicular errors in X!)'});
% % % % % % % % end

%% plot the raw data
figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
axis equal
hold on
plot(bsp_pos.pos(:,1), bsp_pos.pos(:,2), '.', 'color', 0.9.*[1 1 1])
% plot(balls_pos(1,1), balls_pos(1,2), 'ob', 'MarkerFaceColor', 'b')
% plot(balls_pos(2,1), balls_pos(2,2), 'ob', 'MarkerFaceColor', 'b')
plot(balls_pos(:,1), balls_pos(:,2),'o-b', 'LineWidth', 2, 'MarkerFaceColor', 'b')
plot(throws.pos(:,1), throws.pos(:,2), '.k')
switch test_opt
    case 'inside'
        xlim([1460 1470])
        ylim([2499 2501])
    case 'outside'
        xlim([1369 1372])
        ylim([2484 2495])
end
grid on
grid minor
xlabel('X (m)')
ylabel('Y (m)')
title(['Dynamic test raw data - ' test_opt]);
legend({...
    'All measurements';...
    'Balls';...
    'Selected measurements';...
    }, 'Location', 'BestOutside');
filename = sprintf('dynamic_test_raw_%s', test_opt);
fileout = fullfile(results_dir, filename);
saveas(gcf, fileout,'tif')
saveas(gcf, fileout,'fig')








%%







%%
