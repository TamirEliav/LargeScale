function [curvexy, fit_res_all_seg, tunnel_length] = CALIB_fit_curve_200m(pos, smooth_p_turn, smooth_p_ybias)
% pos   calibration localization points
% p     smoothing paramfor the y biases
% calib 

%%
% smooth_p_turn = 0.01;
% smooth_p_ybias = 0.01;

%% 1 - divide the data to 3 segments - long arm / turn / short arm
tunnel_parts_sep = [1365 1356];
turn_margins = [5 -10];
% turn_margins = [0 0];
% we add margins to avoid edge effects in the turn segment. 
% Note we will only use the points in the margins for the initial curve fitting!
IX1 = find( pos.pos(:,1) > tunnel_parts_sep(1));
IX2 = find( pos.pos(:,1) < tunnel_parts_sep(1) + turn_margins(1) &...
            pos.pos(:,1) > tunnel_parts_sep(2) + turn_margins(2) );
IX3 = find( pos.pos(:,1) < tunnel_parts_sep(2));

segments_data_IXs = {IX1, IX2, IX3};
num_seg = length(segments_data_IXs);

%% 2 - fit curves (for linearized position estimation)
smooth_param_linearized = [0 smooth_p_turn 0]; % zero is linear fit
fit_res_all_seg = [];
segments_colors = 'mgb';
figure
hold on
for ii_seg = [2 1 3]
    X = pos.pos(segments_data_IXs{ii_seg},1)';
    Y = pos.pos(segments_data_IXs{ii_seg},2)';
    
    % Set up fittype and options.
    [xData, yData] = prepareCurveData( X, Y );
    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline' );
    opts.SmoothingParam = smooth_param_linearized(ii_seg);
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    fit_res_all_seg{ii_seg} = fitresult;
    
    % plot fitted segments
    plot(X,Y, '.', 'color', segments_colors(ii_seg))
    h=plot( fitresult );
    h.LineWidth = 2;
    h.Color = 'r';
end
xlim([1320 1500])
ylim([2465 2502])
grid on
% axis equal
xlabel('X (m)')
ylabel('Y (m)')
legend({'segments 1 - data';'segments 1 - fit';...
        'segments 2 - data';'segments 2 - fit';...
        'segments 3 - data';'segments 3 - fit';}, 'Location','best')
title('Fitting calibration data', 'FontSize', 12)
    
%% re-sample the tunnel mid-line fit with many points
X_resample = 1500:-0.01:1320;
Y_resample = nan(size(X_resample));
IX1 = find( X_resample > tunnel_parts_sep(1));
IX2 = find( X_resample <= tunnel_parts_sep(1) &...
            X_resample >= tunnel_parts_sep(2) );
IX3 = find( X_resample < tunnel_parts_sep(2));
segments_fit_IXs = {IX1, IX2, IX3};

for ii_seg = 1:num_seg
    seg_IX = segments_fit_IXs{ii_seg};
    Y_resample(seg_IX) = feval(fit_res_all_seg{ii_seg}, X_resample(seg_IX) );
end
curvexy = [X_resample ; Y_resample]';
tunnel_length = sum( sqrt(sum(diff(curvexy).^2,2)) ); % calc tunnel total length

%%
figure
plot(X_resample , Y_resample,'.-' )
xlim([1320 1500])
ylim([2465 2502])
grid on
% axis equal
xlabel('X (m)')
ylabel('Y (m)')
title('Final calibrated curve', 'FontSize', 12)

end






