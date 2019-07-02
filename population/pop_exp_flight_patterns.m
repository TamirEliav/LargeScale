%%
clear
clc

%% load exp summary and choose exps
prm = PARAMS_GetAll();
exp_t = DS_get_exp_summary();
exp_t(~contains(exp_t.recordingArena, '200m'),:) = [];
exp_t(exp_t.position_data_exist==0,:) = [];
exp_t(exp_t.neural_data_exist==0,:) = [];
exp_t(~ismember(exp_t.batNum, [79,148,34,9861,2289] ),:) = [];
exp_t

%%
exps = cellfun(@(c)(exp_load_data(c,'details','flight')), exp_t.exp_ID, 'UniformOutput',0);
exps  = [exps{:}];
nExp = length(exps);

%%
bins_centers = exps(1).flight.speed_traj(1).bins_centers;
nBins = length(bins_centers);
vel_mean_all = nan(nExp, 2, nBins);         % nExp X direction X bins
vel_median_all = nan(nExp, 2, nBins);       % nExp X direction X bins
vel_median_median_all = nan(nExp,2);        % nExp X direction
vel_low_speed_edge_prc_all = nan(nExp,2,2); % nExp X direction X edges
for ii_exp = 1:nExp
    flight = exps(ii_exp).flight;
    for ii_dir = 1:2
        vel_mean_all(ii_exp, ii_dir, :) = exps(ii_exp).flight.speed_traj(ii_dir).vel_mean;
        vel_median_all(ii_exp, ii_dir, :) = exps(ii_exp).flight.speed_traj(ii_dir).vel_median;
        vel_low_speed_edge_prc_all(ii_exp, ii_dir, :) = exps(ii_exp).flight.speed_traj(ii_dir).vel_low_speed_edge_prc;
    end
end

%% figure for histograms of edge vel prc
figure('Units','normalized','Position',[0.2 0.2 0.5 0.5]);
subplot(2,1,1)
hold on
for ii_dir = 1:2
    for ii_edge = 1:2
        histogram(vel_low_speed_edge_prc_all(:,ii_dir,ii_edge),30);
    end
end
grps_str = {'ball1 takeoff';'ball1 landing';'ball2 takeoff';'ball2 landing'};
legend(grps_str);
xlabel('low velocity edge ratio from median velocity')
ylabel('counts');
title('all bats')
subplot(2,1,2)
sdf=reshape(vel_low_speed_edge_prc_all,nExp,4);
[~,lh] = violin(sdf,'xlabel',grps_str');
lh.FontSize = 10;
lh.Position = [0.85 0.4 0.1 0.1];
ylabel({'low velocity edge ratio';'from median velocity'})

figname = 'flight_patterns_low_vel_edge_ratio';
dir_out = 'L:\Analysis\Results\population\flight_pattern';
figname = fullfile(dir_out, figname);
saveas(gcf, figname, 'tif');
saveas(gcf, figname, 'fig');

%% create figure for speed trajectories
figure('Units','normalized','Position',[0 0 1 1]);
pnl = panel();
pnl.pack('v',2);
pnl.margin=30;
pnl.de.margin=10;
h=pnl.title('Velocity profiles');
h.FontSize = 18;
h.Position = [0.5 1.05];
exps_bat = arrayfun(@(x)(x.batNum), [exps.details]);
bats = unique(exps_bat);
colors = brewermap(length(bats),'Set1');
% add low velocity (invalid) area
for ii_pnl = 1:2
    pnl(ii_pnl).select(); hold on
    area([0 prm.fields.valid_speed_pos(1)],   [10 10], 'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 0.5, 'BaseValue',-10, 'LineStyle','none')
    area([prm.fields.valid_speed_pos(2) 200], [10 10], 'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 0.5, 'BaseValue',-10, 'LineStyle','none')
end
exp = exp_load_data(exps(1).details.exp_ID,'details','LM');
h_lgnd1=[];
h_lgnd2=[];
for ii_bat = 1:length(bats)
    bat = bats(ii_bat);
    bat_exp_IX = ismember(exps_bat,bat);
    x = bins_centers;
    c = colors(ii_bat,:);
    
    pnl(1).select(); hold on
    y1 = squeeze(vel_mean_all(bat_exp_IX,1,:));
    y2 = squeeze(vel_mean_all(bat_exp_IX,2,:));
    h = plot([x flip(x)], [y1 fliplr(y2)],'color',c);
    h_lgnd1(ii_bat) = h(1);
    ylabel('Velocity (m/s)')
    
    pnl(2).select(); hold on
    y1 = median(y1);
    y2 = median(y2);
    h = plot([x flip(x)], [y1 flip(y2)],'color',c, 'LineWidth',2);
    h_lgnd2(ii_bat) = h(1);
end
pnl(1).select(); hold on
plot_LM(exp.LM);
pnl(2).select(); hold on
plot_LM(exp.LM);
legend(h_lgnd1,"bat "+bats, 'Location','best');
legend(h_lgnd2,"bat "+bats, 'Location','best');
xlabel('Position (m)')
ylabel('Velocity (m/s)')
linkaxes(pnl.de.axis, 'xy')
% plotbrowser
figname = 'flight_patterns_vel_trajectories';
dir_out = 'L:\Analysis\Results\population\flight_pattern';
figname = fullfile(dir_out, figname);
saveas(gcf, figname, 'tif');
saveas(gcf, figname, 'fig');


%%

%% local function to plot lines for landmarks
function plot_LM(LM)

ylimits = get(gca,'ylim');
for ii_LM=1:length(LM)
    x = LM(ii_LM).pos_proj;
    plot(repelem(x,2) , ylimits, '-', 'color', 0.6.*[1 1 1], 'LineWidth',0.5)
    text(x, ylimits(2)+0.02*diff(ylimits), LM(ii_LM).name, 'Rotation', 45, 'FontSize',6)
end

end

