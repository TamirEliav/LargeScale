%% create all cells spikes data
clear
clc

%%
cells_t = DS_get_cells_summary();
for ii_cell = 1:height(cells_t)
    cell_ID = cells_t.cell_ID{ii_cell};
    fprintf('cell %d/%d %s\n', ii_cell, height(cells_t), cell_ID);
    cell_create_spikes_data(cell_ID);
end
