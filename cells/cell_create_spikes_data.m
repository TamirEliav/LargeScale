function cell_create_spikes_data(cell_ID)

%% load cell/exp details
cell=cell_load_data(cell_ID,'details');
exp=exp_load_data(cell.details.exp_ID, 'path');

%% load cell spikes data from NTT file of the relevant TT
% spike_sorting_dir = 'L:\Analysis\pre_proc\SpikeSorting';
% NTT_file = fullfile(...
%     spike_sorting_dir,...
%     sprintf('%04d',cell.details.bat),...
%     datestr(cell.details.date,'yyyymmdd'),...
%     'spikes_NTT',...
%     ['spikes_' cell.details.TT_ID '.NTT']);

NTT_file = fullfile(exp.path.spikes_sorting,...
                    ['spikes_' cell.details.TT_day_ID '.NTT']);
if isempty(cell.details.stable_ts)
    [Timestamps, CellNumbers, Samples, Header] = ...
         Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 1, [] );
else
    [Timestamps, CellNumbers, Samples, Header] = ...
         Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 4, cell.details.stable_ts );
end

%% parse header
ADBitVolts = sscanf(Header{contains(Header,'ADBitVolts')},'-ADBitVolts %f %f %f %f');

%% convert bits to uVolts
Samples = Samples .* ADBitVolts' .* 1e6;

%% create spikes structure
IX = find(CellNumbers==cell.details.unit);
spikes.ts = Timestamps(IX);
spikes.waveforms = Samples(:,:,IX);
spikes.NTT_file = NTT_file;
 
% add cluster quality info
[spikes.L_Ratio, spikes.Isolation_dis, features_space] = ...
    cell_calc_cluster_quality(Samples, IX);
% spikes.features_space = features_space(IX);

%% save spikes data in mat file
filename_cell_spikes = ['L:\Analysis\Results\cells\spikes\' cell_ID '_cell_spikes'];
save(filename_cell_spikes, 'spikes');


end