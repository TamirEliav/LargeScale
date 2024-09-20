%% 
clear
clc

%%
prm=PARAMS_GetAll();
cells_t = DS_get_cells_summary();

cells_t(~ismember(cells_t.bat, [79,148,34,9861,2289] ),:) = [];

cells_to_exlude = [50 362 390 413 416 319];
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.details];
cells(~contains({cells.brain_area}, {'CA1'})) = [];
cells(~ismember([cells.ClusterQuality], [2])) = [];
cells(ismember([cells.cell_num], cells_to_exlude)) = [];
cells_t = cells_t({cells.cell_ID},:);

%%
cells = cellfun(@(c)(cell_load_data(c,'signif','meanFR','FR_map','details','fields')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
meanFR = [cells.meanFR];
isPyr = [meanFR.all]<=prm.inclusion.interneuron_FR_thr;
cells_pyr = cells(isPyr);
cells_int = cells(~isPyr);

%%
data = struct;

%% Interneurons
maps = cat(1,cells_int.FR_map);
maps = [maps.all];
details = [cells_int.details];
bat_num = [details.bat]';
bat_num = [bat_num bat_num];
cell_IDs = string({details.cell_ID})';
cell_IDs = [cell_IDs cell_IDs];

data.bin_centers = maps(1).bin_centers;
data.INT_maps = cat(1,maps.PSTH);
data.INT_bat_num = [bat_num(:)];
data.INT_cell_IDs = [cell_IDs(:)];

%% Pyramidal
signif = cat(1,cells_pyr.signif);
signif = arrayfun(@(x)(x.TF),signif);
maps = cat(1,cells_pyr.FR_map);
maps = maps(signif);
maps = [maps.all];
fields_per_dir = cat(1,cells_pyr(:).fields);
fields_per_dir = fields_per_dir(signif);
details = [cells_pyr.details];
bat_num = [details.bat]';
bat_num = [bat_num bat_num];
bat_num = bat_num(signif);
cell_IDs = string({details.cell_ID})';
cell_IDs = [cell_IDs cell_IDs];
cell_IDs = cell_IDs(signif);

data.PYT_maps = cat(1,maps.PSTH);
data.PYT_maps01 = zeros(size(data.PYT_maps));
data.PYT_bat_num = bat_num;
data.PYR_cell_IDs = [cell_IDs(:)];
for ii_map = 1:length(maps)
    fields = fields_per_dir{ii_map};
    fields([fields.in_low_speed_area]) = [];
    fields_per_dir{ii_map} = fields;
    fields_edges = cat(1,fields.edges_prc);
    IX = get_data_in_ti(data.bin_centers,fields_edges);
    data.PYT_maps01(ii_map,IX) = 1;
end

%% save data for Nischal
dir_out = 'L:\DATA_for_people\Nischal';
save(fullfile(dir_out,'data_PYR_INT'),'data')

%% copy interneurons figures to a seperate folder for Nachum
% details = [cells_int.details];
% INT_cell_IDs = {details.cell_ID}';
% dir_IN = 'L:\Analysis\Results\cells\figures';
% dir_OUT = 'L:\DATA_for_people\Nischal\INT_figures';
% tmpl_list = INT_cell_IDs;
% ext = 'tif';
% util_copy_files_by_template(dir_IN, dir_OUT, tmpl_list, ext)











%%

