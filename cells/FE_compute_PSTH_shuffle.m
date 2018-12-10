function FE_PSTH_shuffle = FE_compute_PSTH_shuffle(FE_data,n_shuffles, shuffles_max_shift)

rng(0);
FE_PSTH_all_shuffles = {};
parfor ii_shuffle = 1:n_shuffles
    %% create FE from real data with shuffled spikes ts (circular shuffle per flight!)
    FE = FE_data;
    ti = cellfun(@(x,y)([x;y]),{FE.start_ts},{FE.end_ts},'UniformOutput',false);
    shuffles_shifts = rand(1,length(ti)) .* shuffles_max_shift*1e6;
    [FE.spikes_ts] = disperse(cellfun(@shuffle_circular_shift, {FE.spikes_ts}, ti, num2cell(shuffles_shifts),'UniformOutput',false));
    % calc new (shuffled) spikes positions
    spikes_pos = interp1([FE.ts], [FE.pos], [FE.spikes_ts]);
    FE_spikes_IX = mat2cell(1:length([FE.spikes_ts]), 1, cellfun(@length, {FE.spikes_ts}) );
    [FE.spikes_pos] = disperse(cellfun(@(x)(spikes_pos(x)), FE_spikes_IX, 'UniformOutput',false));

    %% calc FR map
    FE_PSTH = FE_compute_PSTH(FE);
    FE_PSTH = rmfield(FE_PSTH, {'spike_density','time_spent','bin_centers','bin_edges'});
    FE_PSTH_all_shuffles{ii_shuffle} = FE_PSTH;
end
FE_PSTH_shuffle.FE_PSTH = [FE_PSTH_all_shuffles{:}];
FE_PSTH_shuffle.shuffles_num = n_shuffles;
FE_PSTH_shuffle.shuffles_max_shift = shuffles_max_shift;

end





