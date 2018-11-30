function cell_calc_FR_map(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','FE');
exp = exp_load_data(cell.details.exp_ID,'details');
prm = PARAMS_GetAll();

%%
for ii_dir = 1:2
    
    %%
    FE = cell.FE{ii_dir};
    FE_full = FE([FE.distance]>100); % TODO: decide what we consider to be a full flight! and use param for 100...
    FE_odd  = FE(1:2:end);
    FE_even = FE(2:2:end);
    FE_begin = FE(1 : round(length(FE)/2)        );
    FE_end   = FE(    round(length(FE)/2)+1 : end);

    FR_map(ii_dir).all   = FE_compute_PSTH(FE);
    FR_map(ii_dir).full  = FE_compute_PSTH(FE_full);
    FR_map(ii_dir).odd   = FE_compute_PSTH(FE_odd);
    FR_map(ii_dir).even  = FE_compute_PSTH(FE_even);
    FR_map(ii_dir).begin = FE_compute_PSTH(FE_begin);
    FR_map(ii_dir).end   = FE_compute_PSTH(FE_end);
    
    % calc correlations between partial subsets
    FR_map(ii_dir).corr_all_full  = FE_PSTH_compute_corr(FR_map(ii_dir).all,   FR_map(ii_dir).full);
    FR_map(ii_dir).corr_odd_even  = FE_PSTH_compute_corr(FR_map(ii_dir).odd,   FR_map(ii_dir).even);
    FR_map(ii_dir).corr_begin_end = FE_PSTH_compute_corr(FR_map(ii_dir).begin, FR_map(ii_dir).end);
    
end

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\FR_map',[cell_ID '_cell_FR_map']);
save(filename, 'FR_map');



end


%%
function FE_PSTH = FE_compute_PSTH(FE)
% TODO: take params from the general params function!!!!!!
    pos = [FE.pos];
    pos_fs = 1e6/min(diff([FE.ts]));
    spikes_pos = [FE.spikes_pos];
    bin_size = 0.1;
    bin_limits = [0 200];
    bin_edges = bin_limits(1):bin_size:bin_limits(end);
    bin_centers = (bin_edges(1:end-1)+bin_edges(1:end-1))./2;
    min_time_spent_per_meter = 0.75;
    min_time_spent_per_bin = min_time_spent_per_meter .* bin_size;
    [PSTH,spike_density,time_spent] = computePSTH(pos,pos_fs,spikes_pos,bin_edges,min_time_spent_per_bin);
    [SI_bits_spike, SI_bits_sec] = computeSI(PSTH,time_spent);
    sparsity = cumputeSparsity(PSTH,time_spent);

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
    
    FE_PSTH.AC = FE_PSTH_compute_AC(FE_PSTH);
end

%%
function PSTH_corr = FE_PSTH_compute_corr(PSTH1,PSTH2)
    [rho,pval] = corr(PSTH1.PSTH', PSTH2.PSTH', 'rows','complete');
    PSTH_corr.rho = rho;
    PSTH_corr.pval = pval;
end

%%
function PSTH_AC = FE_PSTH_compute_AC(PSTH)
    x = PSTH.PSTH;
    x(isnan(x)) = 0;
    [c, lags] = xcorr(x,'coeff');
    lags = lags .* PSTH.bin_size;
    
    [~,~,w,~] = findpeaks(c,'SortStr', 'descend');
    AC_width = w(1) .* PSTH.bin_size;
    
    PSTH_AC.c = c;
    PSTH_AC.lags = lags;
    PSTH_AC.width = AC_width;
    
end










%%
