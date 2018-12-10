%% create all cells spikes data
clear
clc

%% load cells summary and choose cells
cells_t = DS_get_cells_summary();
cells_t(~ismember(cells_t.bat, [79,148,34,9861,2289] ),:) = [];
cells_t(~strcmp(cells_t.brain_area, 'dCA1'),:)=[];

%%
for ii_cell = 1:height(cells_t)
    %%
    cell_ID = cells_t.cell_ID{ii_cell};
    fprintf('cell %d/%d %s\n', ii_cell, height(cells_t), cell_ID);
    
    %%
    try
        
    cell_create_details(cell_ID);
    cell_create_spikes_data(cell_ID);
    cell_calc_time_stability(cell_ID);
    cell_create_flight_data(cell_ID);
    cell_calc_FR_map(cell_ID);
    cell_calc_FR_map_shuffles(cell_ID);
    cell_calc_Ipos(cell_ID);
    cell_calc_fields(cell_ID);
    cell_calc_significant(cell_ID);
    cell_plot_map_fields(cell_ID);
    
    catch err
        disp(err)
    end
    
    close all
    
end
