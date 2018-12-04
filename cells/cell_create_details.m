function cell_create_details(cell_ID)

%% load cells summary
cells_t = DS_get_cells_summary();
details = table2struct( cells_t(cell_ID,:) );

%% add some details...
details.cell_ID = cell_ID;
details.TT_ID = regexprep(cell_ID, '_SS\d+','');
details.exp_ID = regexprep(cell_ID, '_TT\d+_SS\d+','');

%% save cell details
filename_cell_details = ['L:\Analysis\Results\cells\details\' cell_ID '_cell_details' ];
save(filename_cell_details, 'details');

end