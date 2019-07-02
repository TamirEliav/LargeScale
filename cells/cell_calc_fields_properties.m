function cell_calc_fields_properties(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','FE','FR_map','fields');
exp = exp_load_data(cell.details.exp_ID,'details','flight');
prm = PARAMS_GetAll();

%%
directions = [-1 1];
for ii_dir = 1:2
    for ii_field = 1:length(cell.fields{ii_dir})
        field = cell.fields{ii_dir}(ii_field);
        % add velocity info
        direction = directions(ii_dir);
        vel = interp1(exp.flight.speed_traj(ii_dir).bins_centers, ...
                      exp.flight.speed_traj(ii_dir).vel_median,   ...
                      field.loc);
        vel2 = median(field.spikes_vel);
        % mark fields that are entirely in the low speed area
        IX = get_data_in_ti(field.edges_href, prm.fields.valid_speed_pos);
        in_low_speed_area = isempty(IX);
        % update struct
        cell.fields{ii_dir}(ii_field).direction = direction;
        cell.fields{ii_dir}(ii_field).vel = vel;
        cell.fields{ii_dir}(ii_field).vel2 = vel2;
        cell.fields{ii_dir}(ii_field).in_low_speed_area = in_low_speed_area;
    end
end

%% save updated fields results
fields = cell.fields;
filename_cell_fields= ['L:\Analysis\Results\cells\fields\' cell_ID '_cell_fields'];
save(filename_cell_fields, 'fields');


end