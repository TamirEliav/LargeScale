function exp_calc_speed_traj(exp_ID)

%% load exp data
exp = exp_load_data(exp_ID,'pos','flight');
prm = PARAMS_GetAll();

%%
bins_edges = prm.FR_map.bin_limits(1) : prm.FR_map.bin_size : prm.FR_map.bin_limits(2);
bins_centers = (bins_edges(1:end-1) + bins_edges(2:end))/2;
nBins = length(bins_centers);
flight = exp.flight;
directions = [1 -1];
speed_traj = struct();
for ii_dir = 1:2
    IX = [flight.FE.direction] == directions(ii_dir);
    flights = flight.FE(IX);
    IX = discretize([flights.pos], bins_edges);
    vel_mean = accumarray(IX', [flights.vel]',[nBins 1], @mean,nan);
    vel_median = accumarray(IX', [flights.vel]',[nBins 1], @median,nan);
    m = nanmedian(vel_median);
    m1 = interp1(bins_centers, vel_median, prm.fields.valid_speed_pos(1));
    m2 = interp1(bins_centers, vel_median, prm.fields.valid_speed_pos(2));
    
    speed_traj(ii_dir).bins_centers = bins_centers;
    speed_traj(ii_dir).vel_mean = vel_mean;
    speed_traj(ii_dir).vel_median = vel_median;
    speed_traj(ii_dir).vel_median_median = m;
    speed_traj(ii_dir).vel_low_speed_edge = [m1 m2];
    speed_traj(ii_dir).vel_low_speed_edge_prc = [m1 m2]./m;
end

%% save updated flight struct
flight.speed_traj = speed_traj;
file_name = fullfile('L:\Analysis\Results\exp\flight',[exp_ID '_exp_flight']);
save(file_name,'flight');


%%


