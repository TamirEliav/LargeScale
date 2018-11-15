function flight = POS_detect_flight(exp_ID)

%%
duration_min = 10; % seconds
duration_max = 40; % seconds
speed_thr = 4; % m/s
% merge_thr = 0; % not implemented yet, still need to elimintate nans in the position data by interpolation (TODO)

%% read exp info
[exp_path exp_info] = DS_get_path(exp_ID);

%% load pre-processed position data
pos = POS_load(exp_ID);

%%
ts = pos.ts_nlg_usec;
xthr_IX = find( abs(pos.vel_linearized_csaps) > speed_thr );
xthr_ts = ts(xthr_IX);

flight.exp_ID = exp_ID;
% flights.ephocs
% epochs = struct();
% epochs
flight.start_IX = [xthr_IX(1) xthr_IX( find(diff(xthr_IX)>1)+1 )               ];
flight.end_IX   = [           xthr_IX( find(diff(xthr_IX)>1)   )  xthr_IX(end) ];
flight.start_ts = ts(flight.start_IX)';
flight.end_ts = ts(flight.end_IX)';
flight.durations = (flight.end_ts - flight.start_ts) * 1e-6;
flight.direction = sign(pos.pos_linearized(flight.end_IX) - pos.pos_linearized(flight.start_IX))';
flight.distance = abs(pos.pos_linearized(flight.end_IX) - pos.pos_linearized(flight.start_IX))';

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