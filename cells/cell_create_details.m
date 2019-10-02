function cell_create_details(cell_ID)

%% load cells summary
cells_t = DS_get_cells_summary();
details = table2struct( cells_t(cell_ID,:) );

%% add some details...
details.cell_ID = cell_ID;
details.TT_day_ID = regexprep(cell_ID, '_SS\d+','');
details.TT_ID = regexprep(details.TT_day_ID, '_d\d+','');
details.exp_ID = regexprep(cell_ID, '_TT\d+_SS\d+','');
if isempty(details.stable_ts)
    details.stable_ts = [];
else
    details.stable_ts = eval(details.stable_ts);
end
exp = exp_load_data(details.exp_ID,'details');
details.brain_area = exp.details.TT_loc{details.TT};
details.depth = exp.details.depth(details.TT);

%% add TT position 
TT_pos_t = load('L:\TTs_position\atlas_coordinates_Tamir_sep2019.mat');
TT_pos_t = struct2table(TT_pos_t.atlas_coordinates_Tamir);
str = sprintf('bat_%04d_TT%d',details.bat, details.TT);
TT_IX = find(contains(TT_pos_t.name,str));
if isempty(TT_IX)
    details.TT_pos = [nan nan];
else
    details.TT_pos = TT_pos_t(TT_IX,:).proximo_distance;
    details.TT_pos = details.TT_pos{:}';
end

%% save cell details
filename_cell_details = ['L:\Analysis\Results\cells\details\' cell_ID '_cell_details' ];
save(filename_cell_details, 'details');

end