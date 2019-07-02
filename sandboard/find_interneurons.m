%%
clear
clc

%% load cells
cells_t = DS_get_cells_summary();
cells_t(~ismember(cells_t.bat, [79,148,34,9861,2289] ),:) = [];
prm = PARAMS_GetAll();

%% filter cells - brain region / sorting quality
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.details];
cells(~contains({cells.brain_area}, 'CA1')) = [];
cells(~ismember([cells.ClusterQuality], [2])) = [];
% cells(cellfun(@isempty, {cells.stable_ts})) = [];
cells_t = cells_t({cells.cell_ID},:);

%% filter cells - mean FR
cells = cellfun(@(c)(cell_load_data(c,'stats')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.stats];
cells = [cells.all];
cells_t([cells.meanFR_all]<prm.inclusion.interneuron_FR_thr,:) = [];

%% disp final cells table
clear cells
% cells_t
whos cells_t

%% create folder with those cells figures
dir_IN = 'L:\Analysis\Results\cells\figures';
dir_OUT = 'L:\Analysis\Results\cells\interneurons';
tmpl_list = cells_t.cell_ID;
ext = 'tif';
util_copy_files_by_template(dir_IN, dir_OUT, tmpl_list, ext)








%%




%%
