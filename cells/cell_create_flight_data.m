function cell_create_flight_data(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','spikes');
exp = exp_load_data(cell.details.exp_ID,'details','flight','position');
FE = exp.flight.FE;
prm = PARAMS_GetAll();

%% calc pos/vel for each spikes (TODO: maybe that should be done in a seperate function)
cell.spikes.pos = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos_csaps, cell.spikes.ts, 'linear');
cell.spikes.vel = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.vel_csaps, cell.spikes.ts, 'linear');

%% todo: for cell that are partially stable (for part of the total recording), take only the elevant part!!

%%
directions = [1 -1];
for ii_dir = 1:length(directions)
    %%
    FE_dir_IX = find([FE.direction]==directions(ii_dir));
    FE_dir = FE(FE_dir_IX);
    
    %%
    ti = [FE_dir.start_ts; FE_dir.end_ts]';
    [spikes_ts, spikes_IX_per_ti] = get_data_in_ti(cell.spikes.ts, ti);
    spikes_pos = interp1([FE_dir.ts], [FE_dir.pos], spikes_ts, 'linear');
    spikes_vel = interp1([FE_dir.ts], [FE_dir.vel], spikes_ts, 'linear');
    spikes_ts_by_epoch    = cellfun(@(x)(cell.spikes.ts(x)),            spikes_IX_per_ti, 'UniformOutput',false);
    spikes_wvfrm_by_epoch = cellfun(@(x)(cell.spikes.waveforms(:,:,x)), spikes_IX_per_ti, 'UniformOutput',false);
    spikes_pos_by_epoch   = cellfun(@(x)(cell.spikes.pos(x)),           spikes_IX_per_ti, 'UniformOutput',false);
    spikes_vel_by_epoch   = cellfun(@(x)(cell.spikes.vel(x)),           spikes_IX_per_ti, 'UniformOutput',false);
    num_spikes_by_epoch = cellfun(@length, spikes_ts_by_epoch);
    
    [FE_dir(:).spikes_ts] = disperse(spikes_ts_by_epoch);
    [FE_dir(:).spikes_pos] = disperse(spikes_pos_by_epoch);
    [FE_dir(:).spikes_vel] = disperse(spikes_vel_by_epoch);
    [FE_dir(:).spikes_wvfrm] = disperse(spikes_wvfrm_by_epoch);
    [FE_dir(:).num_spikes] = disperse(num_spikes_by_epoch);
    
    %%
    data_FE{ii_dir} = FE_dir;
end

%% change name
FE = data_FE;

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\FE',[cell_ID '_cell_FE']);
save(filename, 'FE');



end



