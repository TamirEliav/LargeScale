%% create all cells details
clear
clc

%%
cells_t = DS_get_cells_summary();
for ii_cell = 1:height(cells_t)
    cell_ID = cells_t.cell_ID{ii_cell};
    cell_create_details(cell_ID);
end