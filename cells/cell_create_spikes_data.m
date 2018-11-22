function cell_create_spikes_data(cell_ID)

%% TODO: TEMP - load data in mat files from Shir
% filename_data_shir = ['L:\Analysis\from_Shir_20181011\Wild_cells_rearranged\', 'cell_data_' cell_ID '.mat'];
% load(filename_data_shir, 'spikes');

%% load cell details
cell_data = cell_load_data(cell_ID,'details');

%% load cell spikes data from NTT file of the relevant TT
spike_sorting_dir = 'L:\Analysis\pre_proc\SpikeSorting';
NTT_file = fullfile(...
    spike_sorting_dir,...
    sprintf('%04d',cell_data.details.bat),...
    datestr(cell_data.details.date,'yyyymmdd'),...
    'spikes_NTT',...
    ['spikes_' cell_data.details.TT_ID '.NTT']);

[Timestamps, CellNumbers, Samples, Header] = ...
     Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 1, [] );

%% create spikes structure
IX = find(CellNumbers==cell_data.details.SS);
spikes.ts = Timestamps(IX);
spikes.waveforms = Samples(:,:,IX);
spikes.NTT_file = NTT_file;
 
% add some further information about the spikes
[spikes.L_Ratio,spikes.Isolation_dis,spikes.features_space] = ...
    cell_calc_cluster_quality(Samples, IX);

%% save spikes data in mat file
filename_cell_spikes = ['L:\Analysis\Results\cells\spikes\' cell_ID '_cell_spikes'];
save(filename_cell_spikes, 'spikes');


end