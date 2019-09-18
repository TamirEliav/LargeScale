function cell_calc_significant(cell_ID)


%% load data
cell = cell_load_data(cell_ID,'details','FR_map','fields','shuffle');
prm = PARAMS_GetAll();

%%
signif = struct();
for ii_dir = 1:2
    has_valid_field   = ~isempty(cell.fields{ii_dir});                                % at least one signif field
    SI_thr_signif     = cell.FR_map(ii_dir).all.SI_bits_spike > prm.signif.SI_thr;    % SI > thr
    SI_shuffle_signif = cell.FR_map(ii_dir).all.SI_bits_spike > ...                   % SI > shuffle
                        prctile([cell.shuffle(ii_dir).FE_PSTH.SI_bits_spike],prm.signif.SI_thr_shuffle);
    TF = true;
    TF = TF & has_valid_field;
    TF = TF & SI_thr_signif;
    TF = TF & SI_shuffle_signif;
    % TODO: what about inclusion criteria? should it go here?!
    
    signif(ii_dir).TF = TF;
    signif(ii_dir).SI_thr = prm.signif.SI_thr;
    signif(ii_dir).SI_thr_shuffle = prm.signif.SI_thr_shuffle;
    signif(ii_dir).SI_thr_signif = SI_thr_signif;
    signif(ii_dir).SI_shuffle_signif = SI_shuffle_signif;
end

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\signif',[cell_ID '_cell_signif']);
save(filename, 'signif');
    













end