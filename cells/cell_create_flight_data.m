function cell_create_flight_data(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','spikes');
exp = exp_load_data(cell.details.exp_ID,'details','flight','pos');
prm = PARAMS_GetAll();

%% calc pos/vel for each spikes
cell.spikes.pos = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos,       cell.spikes.ts, 'linear');
cell.spikes.vel = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.vel_csaps, cell.spikes.ts, 'linear');

%% take only full flights
IX = find([exp.flight.FE.distance] > prm.flight.full_min_distance);
FE_full = exp.flight.FE(IX);

%% for cell that are partially stable (for part of the total recording), take only the elevant part!!
if ~isempty(cell.details.stable_ts)
    FE_ts = [FE_full.start_ts; FE_full.end_ts]';
    valid_FE = all(FE_ts > cell.details.stable_ts(1)  & FE_ts < cell.details.stable_ts(2),2);
    FE_full(~valid_FE) = []; % remove invalid FE
end

%%
FE_by_dir = {};
directions = [1 -1];
for ii_dir = 1:length(directions)
    %%
    FE_dir_IX = find([FE_full.direction]==directions(ii_dir));
    FE = FE_full(FE_dir_IX);
    
    %% add spikes data to FE struct
    ti = [FE.start_ts; FE.end_ts]';
    [~, spikes_IX_per_ti] = get_data_in_ti(cell.spikes.ts, ti);
    spikes_ts_by_epoch    = cellfun(@(x)(cell.spikes.ts(x)),            spikes_IX_per_ti, 'UniformOutput',false);
    spikes_wvfrm_by_epoch = cellfun(@(x)(cell.spikes.waveforms(:,:,x)), spikes_IX_per_ti, 'UniformOutput',false);
    spikes_pos_by_epoch   = cellfun(@(x)(cell.spikes.pos(x)),           spikes_IX_per_ti, 'UniformOutput',false);
    spikes_vel_by_epoch   = cellfun(@(x)(cell.spikes.vel(x)),           spikes_IX_per_ti, 'UniformOutput',false);
    num_spikes_by_epoch = cellfun(@length, spikes_ts_by_epoch);
    
    [FE(:).spikes_ts] = disperse(spikes_ts_by_epoch);
    [FE(:).spikes_pos] = disperse(spikes_pos_by_epoch);
    [FE(:).spikes_vel] = disperse(spikes_vel_by_epoch);
    [FE(:).spikes_wvfrm] = disperse(spikes_wvfrm_by_epoch);
    [FE(:).num_spikes] = disperse(num_spikes_by_epoch);
    
    %%
    FE_by_dir{ii_dir} = FE;
end

%% change name
FE = FE_by_dir;

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\FE',[cell_ID '_cell_FE']);
save(filename, 'FE');



end



