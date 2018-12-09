function POS_detect_flight(exp_ID)

%% load exp data
exp = exp_load_data(exp_ID,'position');
prm = PARAMS_GetAll();

%% arrange relevant data
session_ti = exp_get_sessions_ti(exp_ID,'Behave');
IX = get_data_in_ti(exp.pos.proc_1D.ts, session_ti);
pos.ts = exp.pos.proc_1D.ts(IX);
pos.pos = exp.pos.proc_1D.pos(IX);
pos.vel = exp.pos.proc_1D.vel_csaps(IX);
pos.fs = exp.pos.proc_1D.fs;

%% 1. extract basic epochs (low speed thr crossing)
xthr_IX = find( abs(pos.vel) > prm.flight.speed_low_thr );
start_IX = [xthr_IX(1) xthr_IX( find(diff(xthr_IX)>1)+1 )               ];
end_IX   = [           xthr_IX( find(diff(xthr_IX)>1)   )  xthr_IX(end) ];
start_ts = pos.ts(start_IX)';
end_ts = pos.ts(end_IX)';
duration = (end_ts - start_ts) * 1e-6;
direction = sign(pos.pos(end_IX) - pos.pos(start_IX))';
distance = abs(pos.pos(end_IX) - pos.pos(start_IX))';

% create struct array of flight epochs
FE = repelem(struct,length(start_IX));
[FE(:).start_IX] = disperse(start_IX);
[FE(:).end_IX] = disperse(end_IX);
[FE(:).start_ts] = disperse(start_ts);
[FE(:).end_ts] = disperse(end_ts);
[FE(:).duration] = disperse(duration);
[FE(:).direction] = disperse(direction);
[FE(:).distance] = disperse(distance);

%% add ts/position/velocity of all samples per epoch
ti = [FE.start_ts;FE.end_ts]';
[~, IX_per_ti] = get_data_in_ti(pos.ts,ti);
ts_by_epoch = cellfun(@(x)(pos.ts(x)), IX_per_ti, 'UniformOutput',false);
pos_by_epoch = cellfun(@(x)(pos.pos(x)), IX_per_ti, 'UniformOutput',false);
vel_by_epoch = cellfun(@(x)(pos.vel(x)), IX_per_ti, 'UniformOutput',false);

[FE(:).ts] = disperse(ts_by_epoch);
[FE(:).pos] = disperse(pos_by_epoch);
[FE(:).vel] = disperse(vel_by_epoch);

%% 2. remove epochs without high speed
xthr_IX = find( abs(pos.vel) > prm.flight.speed_high_thr );
xthr_ts = pos.ts(xthr_IX);
ti = [FE.start_ts;FE.end_ts]';
[~, IX_per_ti] = get_data_in_ti(xthr_ts, ti);
[FE(:).duration_high_speed] = disperse( cellfun(@length, IX_per_ti) ./ pos.fs );
IX = find( [FE.duration_high_speed] < prm.flight.high_speed_min_duration );
FE(IX) = [];

%% TODO: remove epochs with extreme speed (or that should be already done in the position processing?!)

%% 3. divide to full/partial epochs
% TODO: decide if to do that here or when calculating the FR map...

%% 4. assign epoch numbers
[FE.number] = disperse(1:length(FE));


%% create flight struct
flight.FE = FE;
flight.speed_low_thr = prm.flight.speed_low_thr;
flight.speed_high_thr = prm.flight.speed_high_thr;
flight.high_speed_min_duration = prm.flight.high_speed_min_duration;

%% save flight struct
file_name = fullfile('L:\Analysis\Results\exp\flight',[exp_ID '_exp_flight']);
save(file_name,'flight');
    

end



