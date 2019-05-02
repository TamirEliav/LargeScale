function cell_create_details(cell_ID)

%% load cells summary
cells_t = DS_get_cells_summary();
details = table2struct( cells_t(cell_ID,:) );

%% add some details...
details.cell_ID = cell_ID;
details.TT_ID = regexprep(cell_ID, '_SS\d+','');
details.exp_ID = regexprep(cell_ID, '_TT\d+_SS\d+','');
if isempty(details.stable_ts)
    details.stable_ts = [];
else
    details.stable_ts = eval(details.stable_ts);
end
exp = exp_load_data(details.exp_ID,'details');
details.brain_area = exp.details.TT_loc{details.TT};

%% save cell details
filename_cell_details = ['L:\Analysis\Results\cells\details\' cell_ID '_cell_details' ];
save(filename_cell_details, 'details');

end