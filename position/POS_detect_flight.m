function POS_detect_flight(exp_ID)

%% params
prm = PARAMS_GetAll();

%% load exp data
exp = exp_load_data(exp_ID);

%% TODO: need to re-check everything after I have the position processing of interpolation and upsampling (+adding NANs!!)

%% get relevant data
session_ti = exp_get_sessions_ti(exp_ID,'Behave');
ts = exp.pos.proc_1D.ts;
IX = get_data_in_ti(ts, session_ti);
pos.ts = exp.pos.proc_1D.ts(IX);
pos.pos = exp.pos.proc_1D.pos_csaps(IX);
pos.vel = exp.pos.proc_1D.vel_csaps(IX);

%% detect speed thr crossing 
xthr_IX = find( abs(pos.vel) > prm.flight.speed_low_thr);
start_IX = [xthr_IX(1) xthr_IX( find(diff(xthr_IX)>1)+1 )               ];
end_IX   = [           xthr_IX( find(diff(xthr_IX)>1)   )  xthr_IX(end) ];
subs = zeros(size(pos.ts));
subs(start_IX) = 1;
subs = cumsum(subs);

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