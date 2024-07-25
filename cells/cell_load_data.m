function cell_data = cell_load_data(cell_ID, varargin)

%% options to load
field_options = {
    'details',
    'spikes',
    'RecStability',
    'FE',
    'FR_map',
    'Ipos',
    'fields',
    'shuffle',
    'signif',
    'meanFR',
    'time_AC',
    'stats',
    'inclusion',
    'cluster_control',
    'cluster_quality',
    'replay_FR_map',
    };

%% parse input
load_all = nargin==1;
if load_all
    fields2load = field_options;
else
	fields2load = varargin;
end

%% allow cell_num as input instead of cell ID
if isnumeric(cell_ID)
    cell_num = cell_ID;
    cells_t = DS_get_cells_summary();
    IX = find(cells_t.cell_num == cell_num);
    if isempty(IX)
        error(sprintf('cell number %d does not exist!',field_name, cell_num))
    end
    cell_ID = cells_t.cell_ID{IX};
end
if ~ischar(cell_ID)
    error("wrong cell_ID input: "+cell_ID)
end

%%
cell_data=struct();
for ii_field = 1:length(fields2load)
    field_name = fields2load{ii_field};
    switch field_name
        case field_options
            file2load = fullfile('L:\Analysis\Results\cells', field_name, [cell_ID '_cell_' field_name '.mat']);
            if exist(file2load,'file')
                field_data = load(file2load);
                cell_data.(field_name) = field_data.(field_name);
            elseif ~load_all
                % raise error only if the user asked for this field
                % specifically, otherwise just load all the available
                % fields
                error(sprintf('field %s data does not exist for cell %s',field_name, cell_ID))
            end
        otherwise
            error(sprintf('field name %s not recognized',field_name));
    end
end

end







