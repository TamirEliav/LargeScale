function POS_detect_flight(exp_ID)

%% load exp data
exp = exp_load_data(exp_ID);
prm = PARAMS_GetAll();

%% TODO: need to re-check everything after I have the position processing of interpolation and upsampling (+adding NANs!!)

%% arrange relevant data
session_ti = exp_get_sessions_ti(exp_ID,'Behave');
ts = exp.pos.proc_1D.ts;
IX = get_data_in_ti(ts, session_ti);
pos.ts = exp.pos.proc_1D.ts(IX);
pos.pos = exp.pos.proc_1D.pos_csaps(IX);
pos.vel = exp.pos.proc_1D.vel_csaps(IX);
pos.fs = exp.pos.proc_1D.fs;

%% detect low speed thr crossing 
xthr = abs(pos.vel) > prm.flight.speed_low_thr;
xthr_IX = find(xthr);
start_IX = [xthr_IX(1) xthr_IX( find(diff(xthr_IX)>1)+1 )               ];
end_IX   = [           xthr_IX( find(diff(xthr_IX)>1)   )  xthr_IX(end) ];
subs = zeros(size(pos.ts));
subs(start_IX) = 1;
subs = cumsum(subs);
subs(~xthr) = nan;
figure
plot(subs,'.')

%% check for high speed thr
min_high_speed_samples = prm.flight.min_duration_high_speed * pos.fs;
IX = ~isnan(subs);
vals = ones(size(subs));
high_speed_num_samples = accumarray(subs(IX)', 1, [], @sum);
high_speed_duration = high_speed_num_samples ./ pos.fs;
HSD=interp1(unique(subs(~isnan(subs))), high_speed_duration, subs,'nearest');

%%
figure
hold on
% scatter(pos.ts, subs, 1, jet(length(subs)), 'filled')
scatter(pos.ts, subs, 10, HSD, 'filled')
colorbar
colormap cool
caxis([min(HSD) max(HSD)])

%% 1. extract basic epochs (low speed thr crossing)
xthr_IX = find( abs(pos.vel) > prm.flight.speed_low_thr );
start_IX = [xthr_IX(1) xthr_IX( find(diff(xthr_IX)>1)+1 )               ];
end_IX   = [           xthr_IX( find(diff(xthr_IX)>1)   )  xthr_IX(end) ];
start_ts = ts(start_IX)';
end_ts = ts(end_IX)';
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

%% 2. remove invalid epochs


%% 3. divide to full/partial epochs

%%
flight.start_IX = [xthr_IX(1) xthr_IX( find(diff(xthr_IX)>1)+1 )               ];
flight.end_IX   = [           xthr_IX( find(diff(xthr_IX)>1)   )  xthr_IX(end) ];
flight.start_ts = ts(flight.start_IX)';
flight.end_ts = ts(flight.end_IX)';
flight.durations = (flight.end_ts - flight.start_ts) * 1e-6;
flight.direction = sign(pos.pos(flight.end_IX) - pos.pos(flight.start_IX))';
flight.distance = abs(pos.pos(flight.end_IX) - pos.pos(flight.start_IX))';

% remove invalid flights
IX = find( flight.durations<duration_min | flight.durations>duration_max );
flight.start_IX(IX) = [];
flight.end_IX(IX) = [];
flight.start_ts(IX) = [];
flight.end_ts(IX) = [];
flight.durations(IX) = [];
flight.direction(IX) = [];
flight.distance(IX) = [];

%% add params
flight.duration_min = duration_min;
flight.duration_max = duration_max;
flight.speed_thr = speed_thr;


%% export the results
switch nargout 
    case 1
        % no need to do anything
    case 0
        % save flights data to file
        mkdir(exp_path.position);
        file_name = fullfile(exp_path.position,['flight_' exp_ID]);
        save(file_name,'flight');
end
    

end