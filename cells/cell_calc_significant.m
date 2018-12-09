function cell_calc_significant(cell_ID)


%% load data
cell = cell_load_data(cell_ID,'details','FR_map','fields','shuffle');
prm = PARAMS_GetAll();

%%
signif = struct();
for ii_dir = 1:2
    TF = true;
    TF = TF & ~isempty(cell.fields{ii_dir});                                % at least one signif field
    TF = TF & cell.FR_map(ii_dir).all.SI_bits_spike > prm.signif.SI_thr;    % SI > thr
    TF = TF & cell.FR_map(ii_dir).all.SI_bits_spike > ...                   % SI > shuffle
        prctile([cell.shuffle(ii_dir).FE_PSTH.SI_bits_spike],prm.signif.SI_thr_shuffle);
    
    signif(ii_dir).TF = TF;
    signif(ii_dir).SI_thr = prm.signif.SI_thr;
    signif(ii_dir).SI_thr_shuffle = prm.signif.SI_thr_shuffle;
    signif(ii_dir).odd_even_FR_map_corr_thr = prm.signif.odd_even_FR_map_corr_thr;
end

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\signif',[cell_ID '_cell_signif']);
save(filename, 'signif');
    













end