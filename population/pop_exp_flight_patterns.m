%%
clear
clc

%% load exp summary and choose exps
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
prm = PARAMS_GetAll();
bins_edges = prm.FR_map.bin_limits(1) : prm.FR_map.bin_size : prm.FR_map.bin_limits(2);
bins_centers = (bins_edges(1:end-1) + bins_edges(2:end))/2;
nBins = length(bins_centers);
vel_mean_all = nan(nExp, 2, nBins); % nExp X direction X bins
vel_median_all = nan(nExp, 2, nBins); % nExp X direction X bins
for ii_exp = 1:nExp
    flight = exps(ii_exp).flight;
    directions = [1 -1];
    for ii_dir = 1:2
        IX = [flight.FE.direction] == directions(ii_dir);
        flights = flight.FE(IX);
        IX = discretize([flights.pos], bins_edges);
        vel_mean = accumarray(IX', [flights.vel]',[nBins 1], @mean,nan);
        vel_median = accumarray(IX', [flights.vel]',[nBins 1], @median,nan);
        vel_mean_all(ii_exp, ii_dir, :) = vel_mean;
        vel_median_all(ii_exp, ii_dir, :) = vel_median;
    end
end

%% create figure
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
% add LM
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
legend(h_lgnd1,"bat "+bats);
legend(h_lgnd2,"bat "+bats);
xlabel('Position (m)')
ylabel('Velocity (m/s)')
linkaxes(pnl.de.axis, 'xy')
% plotbrowser
figname = 'flight_patterns';
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
    plot(repelem(x,2) , ylimits, '-', 'color', 0.9.*[1 1 1], 'LineWidth',0.5)
    text(x, ylimits(2)+0.02*diff(ylimits), LM(ii_LM).name, 'Rotation', 45, 'FontSize',6)
end

end

