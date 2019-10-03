%%
clear
clc

%%
feature = 'SI';
% feature = 'fields_ratio';
% feature = 'fields_smallest';
% feature = 'fields_largest';
% feature = 'fields_number';

dir_IN = 'L:\Analysis\Results\cells\figures\';
dir_OUT = fullfile( dir_IN , feature );
pattern = '*_map_fields.tif';
% dir_OUT = 'L:\Analysis\Results\cells\figures\good_cluster_CA1\by_SI';
% dir_OUT = 'L:\Analysis\Results\cells\figures\good_cluster_CA1\by_field_ratio';
% dir_OUT = 'L:\Analysis\Results\cells\figures\good_cluster_CA1\by_nSpikesAir';
% dir_OUT = 'L:\Analysis\Results\cells\figures\good_cluster_CA1\by_nSpikesAir_dir2';
% dir_OUT = 'L:\Analysis\Results\cells\figures\good_cluster_CA1\by_meanFR_air';
% dir_OUT = 'L:\Analysis\Results\cells\figures\good_cluster_CA1\by_meanFR_all';
% dir_OUT = 'L:\Analysis\Results\cells\figures\good_cluster_CA1\by_nFullFlights_dir1';
mkdir(dir_OUT);

%% 
prm = PARAMS_GetAll();
files = dir(fullfile(dir_IN,pattern));
for ii_file = 1:length(files)
    %%
    file = files(ii_file);
    cell_ID = DS_extract_cell_ID(file.name);
    cell = cell_load_data(cell_ID, 'stats','signif');
    
    
    %% skip cells...
    if cell.stats.all.meanFR_all > prm.inclusion.interneuron_FR_thr
        continue;
    end
    if ~any([cell.signif.TF])
        continue;
    end
    
    %% get feature
    switch feature
        case 'SI'
            SI = max([cell.stats.dir.SI_bits_spike].*[cell.signif.TF]);
            file_new = sprintf('SI_%.2f__%s',SI,file.name);
        case 'fields_smallest'
            fields_smallest = cell.stats.all.field_smallest;
            file_new = sprintf('fields_smallest_%.2f__%s', fields_smallest, file.name);
        case 'fields_largest'
            fields_largest = cell.stats.all.field_largest;
            file_new = sprintf('fields_largest_%.2f__%s', fields_largest, file.name);
        case 'fields_ratio'
            fields_ratio = cell.stats.all.field_ratio_LS;
            file_new = sprintf('fields_ratio_%.2f__%s', fields_ratio,file.name);
        case 'fields_number'
            fields_number = cell.stats.all.field_num;
            file_new = sprintf('fields_number_%.2f__%s', fields_number,file.name);
        case 'num_spikes_air'
            num_spikes_air = min([cell.stats.dir.spikes_num_air]);
            file_new = sprintf('nSpikesAir_%.2f__%s',num_spikes_air,file.name);
        case 'num_spikes_air_dir1'
            num_spikes_air = cell.stats.dir(2).spikes_num_air;
            file_new = sprintf('nSpikesAir_%.2f__%s',num_spikes_air,file.name);
        case 'meanFR_air'
            meanFR_air = cell.stats.all.meanFR_flight;
            file_new = sprintf('meanFR_air_%.2f__%s',meanFR_air ,file.name);
        case 'meanFR_all'
            meanFR_all = cell.stats.all.meanFR_all;
            file_new = sprintf('meanFR_all_%.2f__%s',meanFR_all ,file.name);
        case 'nFullFlights_dir2'
            nFlights = cell.stats.dir(2).num_full_flights;
            file_new = sprintf('nFullFlights_%.2f__%s',nFlights  ,file.name);
        otherwise
            error('not supported')
    end
    
    %% copy file with new prefix
    file_old = fullfile(dir_IN,  file.name);
    file_new = fullfile(dir_OUT, file_new);
    copyfile(file_old, file_new);
end




