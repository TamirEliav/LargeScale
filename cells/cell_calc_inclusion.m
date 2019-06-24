function cell_calc_inclusion(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID, 'stats');
prm = PARAMS_GetAll();

%%
for ii_dir = 1:2
    conditions = [];
    conditions(1) = cell.stats.dir(ii_dir).num_full_flights >= prm.inclusion.min_full_flights;
    conditions(2) = cell.stats.dir(ii_dir).spikes_num_air   >= prm.inclusion.min_spikes_air;
    inclusion(ii_dir).conditions = conditions;
    inclusion(ii_dir).TF = all(conditions);
    inclusion(ii_dir).pyr = cell.stats.all.meanFR_all < prm.inclusion.interneuron_FR_thr;
end

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\inclusion',[cell_ID '_cell_inclusion']);
save(filename, 'inclusion');


end