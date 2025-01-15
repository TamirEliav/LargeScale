%% replay bump vs landmarks/behavior
close all
clear
clc

%%
out_dir = 'E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\cell_revision';
directions_clrs = {[0    0.3843    0.7451];[ 0.5216    0.2471         0]};

%%
% events = load('L:\processed_data_structs\replay_events.mat');
coverage = load('E:\Tamir\work\PROJECTS\LargeScale\paper_replay\data_prepared_for_figures\replay_coverage.mat');

%%
bats_to_analyze = [2382 184];
exp_IX = find(ismember(coverage.T.bat_num, bats_to_analyze));
% exp_IX = 1:height(coverage.T);

%%
clear exp_all
for ii_exp = 1:length(exp_IX)
    exp_ID = coverage.T.exp_ID{exp_IX(ii_exp)};
    exp = exp_load_data(exp_ID,'details','uturns','flight');
    exp_all(ii_exp) = exp;
end
details_all = [exp_all.details];
uturns_all = [exp_all.uturns];
flights_all = [exp_all.flight];
speed_traj_all = cat(1,flights_all.speed_traj);

exp_ID = 'b2382_d190623';
exp = exp_load_data(exp_ID,'LM','rest','flight');
exp.LM(ismember({exp.LM.name},{'ball','enter'}))=[];
exp.LM(end+1).pos_proj = exp.rest.balls_loc(1);
exp.LM(end).name = 'landing-platform ';
exp.LM(end+1).pos_proj = exp.rest.balls_loc(2);
exp.LM(end).name = 'landing-platform ';

%%
y = squeeze(coverage.coverage_all(exp_IX,1:2,:,:));
y = squeeze(mean(y,[1 2 3]));
[pks,locs] = findpeaks(y,'MinPeakProminence',0.1);
x = interp1([1 length(y)], exp.rest.balls_loc, 1:length(y), 'linear');
locs_m = interp1([1 length(y)], exp.rest.balls_loc, locs, 'linear');

%% fit gaussian
xlimits = [0 130];
% Parameters
x0 = locs_m(2);
win = 30;
local_idx = abs(x - x0) <= win/2;

% Local data
x_local = x(local_idx)';
y_local = y(local_idx);

% Fit a Gaussian
gaussEqn = 'a*exp(-((x-b)^2)/(2*c^2))+d';
startPoints = [max(y_local), x0, win/2 median(y)];
fitResult = fit(x_local, y_local, gaussEqn, 'Start', startPoints);

% Extract fit parameters
a = fitResult.a; % Amplitude
b = fitResult.b; % Mean
c = fitResult.c; % Standard deviation
d = fitResult.d; % baseline

exp_func = @(x) builtin('exp', x);
[F,XI] = ksdensity(b,x,"Bandwidth",c);
F = F./max(F)*a+d;

%%
% fig1 = figure(Units="centimeters",Position=[5 5 40 20]);
fig1 = figure(WindowState="maximized");
htl = tiledlayout(5,1,'TileSpacing','tight','Padding','tight');

nexttile
hold on
for ii_exp = 1:size(speed_traj_all,1)
    for ii_dir = 1:2
        traj = speed_traj_all(ii_exp,ii_dir);
        plot(traj.bins_centers, traj.vel_median,'Color',directions_clrs{ii_dir});
    end
end
ylim([-8 8])
title('Speed trajectories (avg per session)')

nexttile
hold on
for ii_exp = 1:size(speed_traj_all,1)
    FE = flights_all(ii_exp).FE;
    [FE([FE.direction]==-1).direction] = deal(2);
    for ii_dir = 1:2
        TF = [FE.direction]==ii_dir;
        plot([FE(TF).pos], [FE(TF).vel], '.', "Color",directions_clrs{ii_dir},'MarkerSize',2);
    end
end
ylim([-8 8])
title('Speed trajectories (all individual flight epochs)')

nexttile
hold on
plot(x,y,'-k')
% plot(x, a * exp_func(-((x - b).^2) / (2 .* c^2))+d, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Gaussian Fit');
plot(XI,F,'r--')
for ii_LM = 1:length(exp.LM)
    xline(exp.LM(ii_LM).pos_proj, 'k-', exp.LM(ii_LM).name)
end
xlim(xlimits)
ylabel('Replay probability density')
% yyaxis right

nexttile
hold on
uturns_pos = [uturns_all.pos];
uturns_dir = [uturns_all.dir];
for ii_dir = 1:2
    histogram(uturns_pos(uturns_dir==ii_dir),BinWidth=1,DisplayStyle="stairs",EdgeColor=directions_clrs{ii_dir});
end
ylabel('counts')
title('uturns position histogram (colors=flight directions)')
for ii_LM = 1:length(exp.LM)
    xline(exp.LM(ii_LM).pos_proj, 'k-', exp.LM(ii_LM).name)
end
xlim(xlimits)

nexttile
hold on
ksdensity([exp.LM.pos_proj],x,"Bandwidth",c);
for ii_LM = 1:length(exp.LM)
    xline(exp.LM(ii_LM).pos_proj, 'k-', exp.LM(ii_LM).name)
end
xlim(xlimits)
xlabel('Position (m)');
title('Landmaks positions convolved with gaussian kernel of width fitted to bump')

linkaxes(htl.Children, 'x')

filename = fullfile(out_dir, 'landmarks_replay_over_representation.pdf');
exportgraphics(fig1, filename)






%%





















%%














%%
