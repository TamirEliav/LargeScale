%% update all cells spikes data (add position/velocity/flight#/...)
clear
clc

%%
cells_t = DS_get_cells_summary();
for ii_cell = 1:height(cells_t)
    %%
    cell_ID = cells_t.cell_ID{ii_cell};
    fprintf('%d/%d: %s\n', ii_cell, height(cells_t), cell_ID);
    
    %%
    cell_update_spikes_data(cell_ID);
end