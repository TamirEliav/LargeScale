function cell_calc_significant(cell_ID)


%% load data
cell = cell_load_data(cell_ID,'details','FR_map','fields','shuffle','FE');
prm = PARAMS_GetAll();

%%
signif = struct();
for ii_dir = 1:2
    % at least one signif field (in valid location)
    if isempty(cell.fields{ii_dir})
        has_valid_field = 0;
    elseif sum(~[cell.fields{ii_dir}.in_low_speed_area]) == 0
        has_valid_field = 0;
    else
        has_valid_field = 1;
    end
    % SI > thr
    SI_thr_signif = cell.FR_map(ii_dir).all.SI_bits_spike > prm.signif.SI_thr;
    % SI > shuffle
    SI_shuffle_signif = cell.FR_map(ii_dir).all.SI_bits_spike > ...                   
                        prctile([cell.shuffle(ii_dir).FE_PSTH.SI_bits_spike],prm.signif.SI_thr_shuffle);
    % has enough full flights
    has_min_flights = length(cell.FE{ii_dir}) >= prm.inclusion.min_full_flights;
    % has enough spikes during flight
    has_min_spikes = length(cell.FE{ii_dir}) >= prm.inclusion.min_spikes_air;
    
    %% apply all conditions
    TF = true;
    TF = TF & SI_thr_signif;
    TF = TF & SI_shuffle_signif;
    TF = TF & has_valid_field;
    TF = TF & has_min_flights; % this is actually an inclusion criteria
    TF = TF & has_min_spikes; % this is actually an inclusion criteria
    
    signif(ii_dir).TF = TF;
    signif(ii_dir).SI_thr = prm.signif.SI_thr;
    signif(ii_dir).SI_thr_shuffle = prm.signif.SI_thr_shuffle;
    signif(ii_dir).SI_thr_signif = SI_thr_signif;
    signif(ii_dir).SI_shuffle_signif = SI_shuffle_signif;
    signif(ii_dir).has_valid_field = has_valid_field;
    signif(ii_dir).has_min_flights = has_min_flights;
    signif(ii_dir).has_min_spikes = has_min_spikes;
    
end

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\signif',[cell_ID '_cell_signif']);
save(filename, 'signif');
    













end