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

    FR_map(ii_dir).all   = FE_PSTH_compute_AC(FE_compute_PSTH(FE));
    FR_map(ii_dir).full  = FE_PSTH_compute_AC(FE_compute_PSTH(FE_full));
    FR_map(ii_dir).odd   = FE_PSTH_compute_AC(FE_compute_PSTH(FE_odd));
    FR_map(ii_dir).even  = FE_PSTH_compute_AC(FE_compute_PSTH(FE_even));
    FR_map(ii_dir).begin = FE_PSTH_compute_AC(FE_compute_PSTH(FE_begin));
    FR_map(ii_dir).end   = FE_PSTH_compute_AC(FE_compute_PSTH(FE_end));
    
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
function PSTH_corr = FE_PSTH_compute_corr(PSTH1,PSTH2)
    [rho,pval] = corr(PSTH1.PSTH', PSTH2.PSTH', 'rows','complete');
    PSTH_corr.rho = rho;
    PSTH_corr.pval = pval;
end

%%
function FE_PSTH = FE_PSTH_compute_AC(FE_PSTH)
    x = FE_PSTH.PSTH;
    x(isnan(x)) = 0;
    [c, lags] = xcorr(x,'coeff');
    lags = lags .* FE_PSTH.bin_size;
    
    [~,~,w,~] = findpeaks(c,'SortStr', 'descend');
    AC_width = w(1) .* FE_PSTH.bin_size;
    
    FE_PSTH.AC.c = c;
    FE_PSTH.AC.lags = lags;
    FE_PSTH.AC.width = AC_width;
end










%%