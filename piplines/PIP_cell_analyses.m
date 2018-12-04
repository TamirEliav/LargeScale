%% create all cells spikes data
clear
clc

%% load cells summary and choose cells
cells = DS_get_cells_summary();
cells(cells.bat~=9861,:)=[];
% cells(~strcmp(cells.brain_area, 'dCA1'),:)=[];

%%
for ii_cell = 1:height(cells)
    %%
    cell_ID = cells.cell_ID{ii_cell};
    fprintf('cell %d/%d %s\n', ii_cell, height(cells), cell_ID);
    
    %%
    try
        
%     cell_create_details(cell_ID);
%     cell_create_spikes_data(cell_ID);
%     cell_calc_time_stability(cell_ID);
%     cell_create_flight_data(cell_ID);
%     cell_calc_FR_map(cell_ID);
%     cell_calc_Ipos(cell_ID);
    cell_calc_fields(cell_ID);
    cell_plot_map_fields(cell_ID);
    
    catch err
        disp(err)
    end
    
    close all
    
end
