function cell_calc_stats(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID);
prm = PARAMS_GetAll();
% load LM data
exp_ID = 'b2289_d180615';
exp = exp_load_data(exp_ID,'LM');
LM = exp.LM;
LM( contains({LM.name},{'ball','enter'}) ) = [];

%% create struct
FR_map = cell.FR_map;
FE = cell.FE;
FE_all = [FE{:}];
fields = cell.fields;
fields_all = [];
for ii_dir = 1:2
    if isempty(fields{ii_dir})
        continue
    end
    fields_to_add = fields{ii_dir};
    % workaround to solve the problem that sometimes I don't have the
    % field 'overlap_edges'... maybe change that in 'cell_calc_fields'...
    if isfield(fields_to_add,'overlap_edges')
        fields_to_add = rmfield(fields_to_add,'overlap_edges');
    end
    fields_all = [fields_all fields_to_add];
end

%% by dir
stats_per_dir = repelem(struct(),2);
for ii_dir = 1:2

    stats_per_dir(ii_dir).spikes_num_air = sum([FE{ii_dir}.num_spikes]);
    stats_per_dir(ii_dir).total_distance = sum([FE{ii_dir}.distance]) .* 1e-3; % in km
    stats_per_dir(ii_dir).time_in_air    = sum([FE{ii_dir}.duration]) ./ 60; % in minutes
    stats_per_dir(ii_dir).num_flights      = length([FE{ii_dir}]);
    stats_per_dir(ii_dir).num_full_flights = sum([FE{ii_dir}.distance] > prm.flight.full_min_distance );

    stats_per_dir(ii_dir).SI_bits_spike  = FR_map(ii_dir).all.SI_bits_spike;
    stats_per_dir(ii_dir).SI_bits_sec    = FR_map(ii_dir).all.SI_bits_sec;
    stats_per_dir(ii_dir).sparsity       = FR_map(ii_dir).all.sparsity;
    stats_per_dir(ii_dir).AC_width       = FR_map(ii_dir).all.AC.width;
    stats_per_dir(ii_dir).corr_odd_even  = FR_map(ii_dir).corr_odd_even.rho;
    stats_per_dir(ii_dir).corr_begin_end = FR_map(ii_dir).corr_begin_end.rho;
%     stats_per_dir(ii_dir).corr_all_full  = FR_map(ii_dir).corr_all_full.rho;
    stats_per_dir(ii_dir).corr_all_full  = nan; % ignore
    stats_per_dir(ii_dir).peakFR         = max(FR_map(ii_dir).all.PSTH);
    stats_per_dir(ii_dir).peakIpos       = max(cell.Ipos.data(ii_dir).Ipos);
    stats_per_dir(ii_dir).map_signif     = cell.signif(ii_dir).TF;
    stats_per_dir(ii_dir).map_signif_shuffle = cell.signif(ii_dir).SI_shuffle_signif;
    
    % fields related
    stats_per_dir(ii_dir).field_num_include_near_balls = length(fields{ii_dir});
    if length(fields{ii_dir}) > 0
        stats_per_dir(ii_dir).spikes_num_field = length([fields{ii_dir}.spikes_ts]);
        valid_speed = ~[fields{ii_dir}.in_low_speed_area];
    else
        stats_per_dir(ii_dir).spikes_num_field = 0;
        valid_speed = [];
    end
    % for scale stats, consider only fields outside the low speed area 
    fields_valid_speed = fields{ii_dir}(valid_speed);
    stats_per_dir(ii_dir).field_num = length(fields_valid_speed);
    switch length(fields_valid_speed)
        case 0 
            stats_per_dir(ii_dir).field_largest     = nan;
            stats_per_dir(ii_dir).field_smallest    = nan;
            stats_per_dir(ii_dir).field_ratio_LS    = nan;
            stats_per_dir(ii_dir).field_CV          = nan;
            stats_per_dir(ii_dir).field_size_mean   = nan;
            stats_per_dir(ii_dir).field_size_median = nan;
            stats_per_dir(ii_dir).field_largest_vel   = nan;
            stats_per_dir(ii_dir).field_largest_vel2  = nan;
            stats_per_dir(ii_dir).field_smallest_vel  = nan;
            stats_per_dir(ii_dir).field_smallest_vel2 = nan;
            stats_per_dir(ii_dir).field_ratio_LS_vel  = nan;
            stats_per_dir(ii_dir).field_ratio_LS_vel2 = nan;
        case 1
            stats_per_dir(ii_dir).field_largest     = nan;
            stats_per_dir(ii_dir).field_smallest    = nan;
            stats_per_dir(ii_dir).field_ratio_LS    = nan;
            stats_per_dir(ii_dir).field_CV          = nan;
            stats_per_dir(ii_dir).field_size_mean   = nanmean([fields_valid_speed.width_prc]);
            stats_per_dir(ii_dir).field_size_median = nanmedian([fields_valid_speed.width_prc]);
            stats_per_dir(ii_dir).field_largest_vel   = nan;
            stats_per_dir(ii_dir).field_largest_vel2  = nan;
            stats_per_dir(ii_dir).field_smallest_vel  = nan;
            stats_per_dir(ii_dir).field_smallest_vel2 = nan;
            stats_per_dir(ii_dir).field_ratio_LS_vel  = nan;
            stats_per_dir(ii_dir).field_ratio_LS_vel2 = nan;
        otherwise % #fields>=2
            [fields_valid_speed.LM_nearest_by_peak]  = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields_valid_speed.loc]  , 'nearest'));
            stats_per_dir(ii_dir).field_largest     = max([fields_valid_speed.width_prc]);
            stats_per_dir(ii_dir).field_smallest    = min([fields_valid_speed.width_prc]);
            stats_per_dir(ii_dir).field_ratio_LS    = max([fields_valid_speed.width_prc]) / min([fields_valid_speed.width_prc]);
            stats_per_dir(ii_dir).field_CV          = nanstd([fields_valid_speed.width_prc]) / nanmean([fields_valid_speed.width_prc]);
            stats_per_dir(ii_dir).field_size_mean   = nanmean([fields_valid_speed.width_prc]);
            stats_per_dir(ii_dir).field_size_median = nanmedian([fields_valid_speed.width_prc]);
            stats_per_dir(ii_dir).field_largest_vel   = max([fields_valid_speed.vel]);
            stats_per_dir(ii_dir).field_largest_vel2  = max([fields_valid_speed.vel2]);
            stats_per_dir(ii_dir).field_smallest_vel  = min([fields_valid_speed.vel]);
            stats_per_dir(ii_dir).field_smallest_vel2 = min([fields_valid_speed.vel2]);
            stats_per_dir(ii_dir).field_ratio_LS_vel  = stats_per_dir(ii_dir).field_largest_vel  / stats_per_dir(ii_dir).field_smallest_vel;
            stats_per_dir(ii_dir).field_ratio_LS_vel2 = stats_per_dir(ii_dir).field_largest_vel2 / stats_per_dir(ii_dir).field_smallest_vel2;
    end
    stats_per_dir(ii_dir).spikes_prc_field = 100.* stats_per_dir(ii_dir).spikes_num_field / stats_per_dir(ii_dir).spikes_num_air;
end

%% all
stats_all = struct();

stats_all.IsoDist = cell.spikes.Isolation_dis;
stats_all.L_Ratio = cell.spikes.L_Ratio;
stats_all.meanFR_all     = cell.meanFR.all;
stats_all.meanFR_flight  = cell.meanFR.in_flight;
stats_all.num_flights      = length(FE_all);
stats_all.num_full_flights = sum([FE_all.distance] > prm.flight.full_min_distance );

stats_all.spikes_num_air = sum([stats_per_dir.spikes_num_air]);
stats_all.total_distance = sum([stats_per_dir.total_distance]);
stats_all.time_in_air    = sum([stats_per_dir.time_in_air]);
stats_all.spikes_num_field = sum([stats_per_dir.spikes_num_field]);
stats_all.spikes_prc_field = 100.* stats_all.spikes_num_field / stats_all.spikes_num_air;

% fields related
stats_all.field_num_include_near_balls = length(fields_all);
if length(fields_all) > 0
    valid_speed = ~[fields_all.in_low_speed_area];
else
    valid_speed = [];
end
% for scale stats, consider only fields outside the low speed area 
fields_valid_speed = fields_all(valid_speed);
stats_all.field_num = length(fields_valid_speed);
stats_all.field_largest     = nan;
stats_all.field_smallest    = nan;
stats_all.field_ratio_LS    = nan;
stats_all.field_CV          = nan;
stats_all.field_size_mean   = nan;
stats_all.field_size_median = nan;
stats_all.field_largest_vel   = nan;
stats_all.field_largest_vel2  = nan;
stats_all.field_smallest_vel  = nan;
stats_all.field_smallest_vel2 = nan;
stats_all.field_ratio_LS_vel  = nan;
stats_all.field_ratio_LS_vel2 = nan;
stats_all.field_largest_dist2LM  = nan;
stats_all.field_smallest_dist2LM = nan;
switch length(fields_valid_speed)
    case 0
        % do nothing
    case 1
        stats_all.field_size_mean   = nanmean([fields_valid_speed.width_prc]);
        stats_all.field_size_median = nanmedian([fields_valid_speed.width_prc]);
    otherwise % #fields>=2
        [fields_valid_speed.LM_nearest_by_peak]  = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields_valid_speed.loc]  , 'nearest'));
        [stats_all.field_largest,  L_IX] = max([fields_valid_speed.width_prc]);
        [stats_all.field_smallest, S_IX] = min([fields_valid_speed.width_prc]);
        stats_all.field_ratio_LS    = max([fields_valid_speed.width_prc]) / min([fields_valid_speed.width_prc]);
        stats_all.field_CV          = nanstd([fields_valid_speed.width_prc]) / nanmean([fields_valid_speed.width_prc]);
        stats_all.field_size_mean   = nanmean([fields_valid_speed.width_prc]);
        stats_all.field_size_median = nanmedian([fields_valid_speed.width_prc]);
        stats_all.field_largest_vel   = fields_valid_speed(L_IX).vel;
        stats_all.field_largest_vel2  = fields_valid_speed(L_IX).vel2;
        stats_all.field_smallest_vel  = fields_valid_speed(S_IX).vel;
        stats_all.field_smallest_vel2 = fields_valid_speed(S_IX).vel2;
        stats_all.field_ratio_LS_vel  = stats_all.field_largest_vel  / stats_all.field_smallest_vel;
        stats_all.field_ratio_LS_vel2 = stats_all.field_largest_vel2 / stats_all.field_smallest_vel2;
        stats_all.field_largest_dist2LM  = abs(fields_valid_speed(L_IX).loc - fields_valid_speed(L_IX).LM_nearest_by_peak);
        stats_all.field_smallest_dist2LM = abs(fields_valid_speed(S_IX).loc - fields_valid_speed(S_IX).LM_nearest_by_peak);
end

%% combine
stats.all = stats_all;
stats.dir = stats_per_dir;

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\stats',[cell.details.cell_ID '_cell_stats']);
save(filename, 'stats');



end






