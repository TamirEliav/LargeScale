function cell_calc_time_stability(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','spikes');
exp = exp_load_data(cell.details.exp_ID,'details');
prm = PARAMS_GetAll();

%%
ts_limits = exp.details.session_ts([1 end]);
edges = ts_limits(1) : (prm.RecStability.BinSize*1e6*60) : ts_limits(end);
N = histcounts(cell.spikes.ts,edges);
FR = N ./ (prm.RecStability.BinSize*60);
bin_centers = (edges(1:end-1)+edges(1:end-1)) ./ 2;

%% create struct
RecStability.BinSize = prm.RecStability.BinSize;
RecStability.bin_centers = bin_centers;
RecStability.FR = FR;

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\RecStability',[cell_ID '_cell_RecStability']);
save(filename, 'RecStability');



end