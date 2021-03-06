function cell_calc_time_stability(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','spikes');
exp = exp_load_data(cell.details.exp_ID,'details');
prm = PARAMS_GetAll();

%%
[N,edges] = histcounts(cell.spikes.ts,...
                       'BinWidth',prm.RecStability.BinSize*1e6*60,...
                       'Normalization','count');
FR = N ./ (prm.RecStability.BinSize*60);
bin_centers = edges2centers(edges);

%% create struct
RecStability.BinSize = prm.RecStability.BinSize;
RecStability.bin_centers = bin_centers;
RecStability.FR = FR;

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\RecStability',[cell_ID '_cell_RecStability']);
save(filename, 'RecStability');



end