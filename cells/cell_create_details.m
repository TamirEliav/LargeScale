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
% TT_pos_t = load('L:\TTs_position\atlas_coordinates_Tamir_sep2019.mat');
% TT_pos_t = load('L:\TTs_position\atlas_coordinates_Tamir_Nov2019.mat');
TT_pos_t = load('L:\TTs_position\atlas_coordinates_Tamir_Sep2020.mat');
TT_pos_t = struct2table(TT_pos_t.atlas_coordinates_Tamir);
str = sprintf('bat_%04d_TT%d',details.bat, details.TT);
TT_IX = find(contains(TT_pos_t.name,str));
if isempty(TT_IX)
    details.TT_pos = [nan nan];
    details.TT_pos_proximodistal_prc = nan;
    details.TT_pos_longitudinal_prc = nan;
else
    TT_pos_t = TT_pos_t(TT_IX,:);
    details.TT_pos = TT_pos_t.proximo_distance;
    details.TT_pos = details.TT_pos{:}';
%     details.TT_pos_ML = TT_pos_t.ML_position;
%     details.TT_pos_rel_AP = TT_pos_t.rel_AP_position;
%     details.TT_pos_rel_DV = TT_pos_t.DV_position;
%     details.TT_pos_slide = TT_pos_t.slide;
    details.TT_pos_proximodistal_prc = TT_pos_t.proximo_distal_TT_precentage;
    details.TT_pos_longitudinal_prc = TT_pos_t.longitudinal_TT_precentage;
end
details.TT_pos_prc = [details.TT_pos_proximodistal_prc details.TT_pos_longitudinal_prc];

%% save cell details
filename_cell_details = ['L:\Analysis\Results\cells\details\' cell_ID '_cell_details' ];
save(filename_cell_details, 'details');

end