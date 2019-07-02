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
cells_t([cells.meanFR_all]>prm.inclusion.interneuron_FR_thr,:) = [];

%% disp final cells table
% cells_t
whos cells_t


%% load cells FR maps
cells = cellfun(@(c)(cell_load_data(c,'fields')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
% sdf = arrayfun(@(c)([c.fields{:}]), cells, 'UniformOutput',0);

%% find cells with field with less than 10 spikes
list_cells_less_10_spikes = {};
for ii_cell = 1:length(cells)
    found_field = 0;
    for ii_dir = 1:2
        fields = cells(ii_cell).fields{ii_dir};
        if isempty(fields)
            continue;
        end
        nSpikes = cellfun(@length, {fields.spikes_ts});
        if any(nSpikes <= 10)
            found_field = 1;
        end
    end
    if found_field 
        list_cells_less_10_spikes{end+1} = cells_t.cell_ID{ii_cell};
    end
end

%% create folder with those cells figures
dir_IN = 'L:\Analysis\Results\cells\figures';
dir_OUT = 'L:\Analysis\Results\cells\figures\with_field_less_than_10_spikes';
tmpl_list = list_cells_less_10_spikes;
ext = 'tif';
util_copy_files_by_template(dir_IN, dir_OUT, tmpl_list, ext)








%%




%%
