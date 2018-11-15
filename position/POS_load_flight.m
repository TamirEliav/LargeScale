function flight = POS_load_flight(exp_ID)

%% get exp info
[exp_path exp_info] = DS_get_path(exp_ID);

%%
file_name = fullfile(exp_path.position,['flight_' exp_ID]);
load(file_name)

end