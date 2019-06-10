function cell_calc_mean_FR(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','spikes','FE');
exp = exp_load_data(cell.details.exp_ID,'details');
prm = PARAMS_GetAll();

%% mean FR - all stable recordings
% note, I take the first/last spike as the start/end of time interval as:
% 1) it should correspond to the stability ts interval
% 2) it shouldn't make a big difference in case the cell was not firing from
% the very beggining of the recordings to the very end of it (and if this
% is the case we also can't be sure we were actually recording the cell...)
nSpikes = length(cell.spikes.ts);
T = range(cell.spikes.ts) * 1e-6;
mean_FR_all = nSpikes / T;

%% mean FR in sessions
ti = exp.details.session_ts;
[IX, IX_per_ti] = get_data_in_ti(cell.spikes.ts, exp.details.session_ts);
num_spikes_per_session = cellfun(@length,IX_per_ti);
mean_FR_per_session = num_spikes_per_session' ./ (diff(ti,1,2)*1e-6);
mean_FR_all_sessions = sum(num_spikes_per_session) / sum(diff(ti,1,2)*1e-6);
% TODO: calc mean FR only for the part in the behave session where the cell
% was stable! (after I have the definition of stable timestamps)

%% mean FR in flight only
total_time_in_flight   = sum(    [cell.FE{1}.duration   cell.FE{2}.duration] );
total_spikes_in_flight = length( [cell.FE{1}.spikes_ts  cell.FE{2}.spikes_ts] );
mean_FR_in_flight = total_spikes_in_flight / total_time_in_flight;

%% create struct
meanFR.all = mean_FR_all;
meanFR.per_session = mean_FR_per_session;
meanFR.all_sessions = mean_FR_all_sessions;
meanFR.in_flight = mean_FR_in_flight;

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\meanFR',[cell_ID '_cell_meanFR']);
save(filename, 'meanFR');



end