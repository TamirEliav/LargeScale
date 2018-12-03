function cell_calc_fields(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','FE','FR_map');
exp = exp_load_data(cell.details.exp_ID,'details');
prm = PARAMS_GetAll();

%%
fields_per_dir = {};
for ii_dir = 1:2
    
    %% arrange data
    FR_map = cell.FR_map(ii_dir).all;
    FE = cell.FE{ii_dir};
    
    %% detect peaks
    fields = fields_detect_peaks(FR_map,prm);
    
    %% TODO: remove fields with not enough spikes
    
    %% Field width/edges
    [widths, edges] = fields_calc_width_edges(FR_map, fields, prm.fields.width_href);
    [fields(:).width_href]   = disperse(widths);
    [fields(:).edges_href]   = disperse(edges');
    
    %% remove overlapping, lower fields
    fields = fields_remove_overlaps(FR_map, fields, prm);
    
    %% limit field size if neighboring fields overlap (in field width and not in ref height for overlap)
    fields = fields_limit_width_by_neighbors(FR_map, fields, prm);
    
    %% get spikes in field
    fields = fields_add_spikes_data(FE, FR_map, fields, prm);
    
    %% remove unstable fields
    fields = fields_remove_unstable(fields, prm);
    
    %% Remove non significant fields based on shuffling analysis.
    

    %%
    [~,IX] = sort([fields.loc],'ascend');
    fields = fields(IX);
    fields_per_dir{ii_dir} = fields;
    
end

%% save cell fields
fields = fields_per_dir; % change name before saving
filename_cell_fields= ['L:\Analysis\Results\cells\fields\' cell_ID '_cell_fields' ];
save(filename_cell_fields, 'fields');

end


%%



%%
