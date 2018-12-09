function cell_calc_FR_map_shuffles(cell_ID)

%% load cell data
cell = cell_load_data(cell_ID,'details','FE');
exp = exp_load_data(cell.details.exp_ID,'details','position');
prm = PARAMS_GetAll();

%%
for ii_dir = 1:2 
    FE_PSTH_shuffle(ii_dir) = FE_compute_PSTH_shuffle(cell.FE{ii_dir},...
                                                      prm.FR_map.shuffles_num,...
                                                      prm.FR_map.shuffles_max_shift);
end

shuffle = FE_PSTH_shuffle;


%% save data to file
filename = fullfile('L:\Analysis\Results\cells\shuffle',[cell_ID '_cell_shuffle']);
save(filename, 'shuffle');
    


end