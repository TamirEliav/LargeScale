function exp_calc_speed_traj(exp_ID)

%% load exp data
exp = exp_load_data(exp_ID,'pos','flight');
prm = PARAMS_GetAll();

%%
bins_edges = prm.FR_map.bin_limits(1) : prm.FR_map.bin_size : prm.FR_map.bin_limits(2);
bins_centers = (bins_edges(1:end-1) + bins_edges(2:end))/2;
nBins = length(bins_centers);
FE = exp.flight.FE;
FE = FE([FE.distance]>prm.flight.full_min_distance); % take only full filghts
directions = [1 -1];
speed_traj = struct();
for ii_dir = 1:2
    dir_IX = [FE.direction] == directions(ii_dir);
    flights = FE(dir_IX);

    pos_all = [flights.pos];
    vel_all = [flights.vel];
    speed_all = abs(vel_all);
    
    subs = discretize([flights.pos], bins_edges);
    vel_mean   = accumarray(subs', vel_all', [nBins 1], @mean,   nan);
    vel_median = accumarray(subs', vel_all', [nBins 1], @median, nan);
    vel_std    = accumarray(subs', vel_all', [nBins 1], @std,    nan);

    speed_cv.raw = std(speed_all) / mean(speed_all);
    speed_cv.across_pos = nanstd(vel_median) / abs(nanmean(vel_median));
    IX = get_data_in_ti(pos_all, prm.fields.valid_speed_pos);
    speed_cv.raw_high_speed = std(speed_all(IX)) / mean(speed_all(IX));
    IX = get_data_in_ti(bins_centers, prm.fields.valid_speed_pos);
    speed_cv.across_pos_high_speed =  nanstd(vel_median(IX)) / abs(nanmean(vel_median(IX)));
    
    m = nanmedian(vel_median);
    m1 = interp1(bins_centers, vel_median, prm.fields.valid_speed_pos(1)); % ball 1
    m2 = interp1(bins_centers, vel_median, prm.fields.valid_speed_pos(2)); % ball 2
    
    speed_traj(ii_dir).bins_centers = bins_centers;
    speed_traj(ii_dir).vel_mean = vel_mean;
    speed_traj(ii_dir).vel_median = vel_median;
    speed_traj(ii_dir).vel_median_median = m;
    speed_traj(ii_dir).vel_low_speed_edge = [m1 m2];
    speed_traj(ii_dir).vel_low_speed_edge_prc = [m1 m2]./m;
    speed_traj(ii_dir).speed_cv = speed_cv;
end

%% save updated flight struct
flight = exp.flight;
flight.speed_traj = speed_traj;
file_name = fullfile('L:\Analysis\Results\exp\flight',[exp_ID '_exp_flight']);
save(file_name,'flight');


%%


