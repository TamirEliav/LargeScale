function cell_calc_stats(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID);
exp = exp_load_data(cell.details.exp_ID,'details');
prm = PARAMS_GetAll();

%% create struct
FR_map = cell.FR_map;
FE = cell.FE;
fields = cell.fields;
fields_all = [];
for ii_dir = 1:2
    if isempty(fields{ii_dir})
        continue
    end
    fields_to_add = fields{ii_dir};
    % TODO: workaround to solve the problem that sometimes I don't have the
    % field 'overlap_edges'... maybe change that in 'cell_calc_fields'...
    if isfield(fields_to_add,'overlap_edges')
        fields_to_add = rmfield(fields_to_add,'overlap_edges');
    end
    fields_all = [fields_all fields_to_add];
end

%%
stats_per_dir = repelem(struct(),2);
for ii_dir = 1:2

    stats_per_dir(ii_dir).spikes_num_air = sum([FE{ii_dir}.num_spikes]);
    stats_per_dir(ii_dir).total_distance = sum([FE{ii_dir}.distance]) .* 1e-3; % in km
    stats_per_dir(ii_dir).time_in_air    = sum([FE{ii_dir}.duration]) ./ 60; % in minutes

    stats_per_dir(ii_dir).SI_bits_spike  = FR_map(ii_dir).all.SI_bits_spike;
    stats_per_dir(ii_dir).SI_bits_sec    = FR_map(ii_dir).all.SI_bits_sec;
    stats_per_dir(ii_dir).sparsity       = FR_map(ii_dir).all.sparsity;
    stats_per_dir(ii_dir).AC_width       = FR_map(ii_dir).all.AC.width;
    stats_per_dir(ii_dir).corr_odd_even  = FR_map(ii_dir).corr_odd_even.rho;
    stats_per_dir(ii_dir).corr_begin_end = FR_map(ii_dir).corr_begin_end.rho;
    stats_per_dir(ii_dir).corr_all_full  = FR_map(ii_dir).corr_all_full.rho;
    stats_per_dir(ii_dir).peakFR         = max(FR_map(ii_dir).all.PSTH);
    stats_per_dir(ii_dir).map_signif     = cell.signif(ii_dir).TF;
    stats_per_dir(ii_dir).peakIpos       = max(cell.Ipos.data(ii_dir).Ipos);
    
    stats_per_dir(ii_dir).field_num      = length(fields{ii_dir});
    if length(fields{ii_dir}) > 0 % only if there are fields
        stats_per_dir(ii_dir).field_largest  = max([fields{ii_dir}.width_prc]);
        stats_per_dir(ii_dir).field_smallest = min([fields{ii_dir}.width_prc]);
        stats_per_dir(ii_dir).field_ratio_LS = max([fields{ii_dir}.width_prc]) / min([fields{ii_dir}.width_prc]);
        stats_per_dir(ii_dir).field_CV       = nanstd([fields{ii_dir}.width_prc]) / nanmean([fields{ii_dir}.width_prc]);
        stats_per_dir(ii_dir).spikes_num_field = length([fields{ii_dir}.spikes_ts]);
    else
        stats_per_dir(ii_dir).field_largest  = nan;
        stats_per_dir(ii_dir).field_smallest = nan;
        stats_per_dir(ii_dir).field_ratio_LS = nan;
        stats_per_dir(ii_dir).field_CV       = nan;
        stats_per_dir(ii_dir).spikes_num_field = 0;
    end
    stats_per_dir(ii_dir).spikes_prc_field = 100.* stats_per_dir(ii_dir).spikes_num_field / stats_per_dir(ii_dir).spikes_num_air;
end

%%
stats_all = struct();

stats_all.IsoDist = cell.spikes.Isolation_dis;
stats_all.L_Ratio = cell.spikes.L_Ratio;
stats_all.meanFR_all = cell.meanFR.all_sessions;
stats_all.meanFR_flight = cell.meanFR.in_flight;

stats_all.spikes_num_air = sum([stats_per_dir.spikes_num_air]);
stats_all.total_distance = sum([stats_per_dir.total_distance]);
stats_all.time_in_air    = sum([stats_per_dir.time_in_air]);
stats_all.spikes_num_field = sum([stats_per_dir.spikes_num_field]);
stats_all.spikes_prc_field = 100.* stats_all.spikes_num_field / stats_all.spikes_num_air;

stats_all.field_num      = length(fields_all);
if length(fields_all) > 0 % only if there are fields
    stats_all.field_largest  = max([fields_all.width_prc]);
    stats_all.field_smallest = min([fields_all.width_prc]);
    stats_all.field_ratio_LS = max([fields_all.width_prc]) / min([fields_all.width_prc]);
    stats_all.field_CV       = nanstd([fields_all.width_prc]) / nanmean([fields_all.width_prc]);
else
    stats_all.field_largest  = nan;
    stats_all.field_smallest = nan;
    stats_all.field_ratio_LS = nan;
    stats_all.field_CV       = nan;
end

%% combine
stats.all = stats_all;
stats.dir = stats_per_dir;

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\stats',[cell_ID '_cell_stats']);
save(filename, 'stats');



end





