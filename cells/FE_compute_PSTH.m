function FE_PSTH = FE_compute_PSTH(FE)

%%
if isempty(FE)
    pos_fs = 100;
else
    pos_fs = 1e6/median(diff([FE.ts]));
end

%%
prm = PARAMS_GetAll();

%%
pos = [FE.pos];
spikes_pos = [FE.spikes_pos];
bin_size = prm.FR_map.bin_size;
bin_limits = prm.FR_map.bin_limits;
bin_edges = bin_limits(1):bin_size:bin_limits(end);
bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))./2;
min_time_spent_per_meter = prm.FR_map.min_time_spent_per_meter;
min_time_spent_per_bin = prm.FR_map.min_time_spent_per_meter .* bin_size;
ker_SD = prm.FR_map.ker_SD;
[PSTH,spike_density,time_spent] = computePSTH(pos,pos_fs,spikes_pos,bin_edges,min_time_spent_per_bin,ker_SD);
[SI_bits_spike, SI_bits_sec] = computeSI(PSTH,time_spent);
sparsity = cumputeSparsity(PSTH,time_spent);

%%
FE_PSTH.PSTH = PSTH;
FE_PSTH.spike_density = spike_density;
FE_PSTH.time_spent = time_spent;
FE_PSTH.min_time_spent_per_meter = min_time_spent_per_meter;
FE_PSTH.bin_size = bin_size;
FE_PSTH.bin_edges = bin_edges;
FE_PSTH.bin_centers = bin_centers;
FE_PSTH.SI_bits_spike = SI_bits_spike;
FE_PSTH.SI_bits_sec = SI_bits_sec;
FE_PSTH.sparsity  = sparsity;

end
