function pos_holes_stats = exp_calc_position_holes_stats(exp_ID)

%%
disp(exp_ID)

%% load exp data
exp = exp_load_data(exp_ID,'details','pos','flight');
prm = PARAMS_GetAll();

%%
FE = [exp.flight.FE];
FE([FE.distance]<prm.flight.full_min_distance) = []; % remove short flights
FE_ts = [[FE.start_ts];[FE.end_ts]]';
IX = get_data_in_ti(exp.pos.raw.ts_nlg_usec', FE_ts);
rawdata_fraction = length(IX)/exp.pos.raw.fs / sum([FE.duration]);

%% create results struct
pos_holes_stats = struct();
pos_holes_stats.rawdata_fraction = rawdata_fraction;


end




