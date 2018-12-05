function cell_calc_FR_map_shuffles(cell_ID)

%% load cell data
cell = cell_load_data(cell_ID,'details','FE');
exp = exp_load_data(cell.details.exp_ID,'details','position');
prm = PARAMS_GetAll();

%%
rng(0);
shuffle_by_dir = {};
for ii_dir = 1:2
    
    %%
    FE_PSTH_all_shuffles = {};
    for ii_shuffle = 1:prm.FR_map.shuffles_num
        %% create FE from real data with shuffled spikes ts (circular shuffle per flight!)
        FE = cell.FE{ii_dir};
        ti = cellfun(@(x,y)([x;y]),{FE.start_ts},{FE.end_ts},'UniformOutput',false);
        shuffles_shifts = rand(1,length(ti)) .* prm.FR_map.shuffles_max_shift*1e6;
        [FE.spikes_ts] = disperse(cellfun(@shuffle_circular_shift, {FE.spikes_ts}, ti, num2cell(shuffles_shifts),'UniformOutput',false));
        % calc new spikes positions / velocity
        spikes_pos = interp1([FE.ts], [FE.pos], [FE.spikes_ts]);
        spikes_vel = interp1([FE.ts], [FE.vel], [FE.spikes_ts]);
        FE_spikes_IX = mat2cell(1:length([FE.spikes_ts]), 1, cellfun(@length, {FE.spikes_ts}) );
        [FE.spikes_pos] = disperse(cellfun(@(x)(spikes_pos(x)), FE_spikes_IX, 'UniformOutput',false));
        [FE.spikes_vel] = disperse(cellfun(@(x)(spikes_vel(x)), FE_spikes_IX, 'UniformOutput',false));
        
        %% calc FR map
        FE_PSTH = FE_compute_PSTH(FE);
        FE_PSTH = rmfield(FE_PSTH, {'spike_density','time_spent','bin_centers','bin_edges'});
        FE_PSTH_all_shuffles{ii_shuffle} = FE_PSTH;
    end
    shuffle_by_dir{ii_dir}.FE_PSTH = [FE_PSTH_all_shuffles{:}];
    shuffle_by_dir{ii_dir}.shuffles_num = prm.FR_map.shuffles_num;
    shuffle_by_dir{ii_dir}.shuffle_max_shift = prm.FR_map.shuffles_max_shift;
    
end

shuffle = shuffle_by_dir;


%% save data to file
filename = fullfile('L:\Analysis\Results\cells\shuffle',[cell_ID '_cell_shuffle']);
save(filename, 'shuffle');
    


end