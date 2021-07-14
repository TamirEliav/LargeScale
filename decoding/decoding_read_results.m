%%
clear
clc

%%
exp_ID = 'b9861_d180527';
dir_IN = 'D:\sequences\seq_uri_eden\decoded';
filename = fullfile(dir_IN, [exp_ID '_sleep_res.h5']);
exp = exp_load_data(exp_ID, 'details','path','PE','MUA');

%% read
pos = h5read(filename,'/position');
posterior = h5read(filename,'/acausal_posterior');
state = h5read(filename,'/state');
time = h5read(filename,'/time');

%%
posterior_state = squeeze(sum(posterior,1));
posterior_pos = squeeze(sum(posterior,2));
[~,MAP_pos_IX] = max(posterior_pos,[],1);
[~,MAP_state_IX] = max(posterior_state,[],1);
MAP_pos = pos(MAP_pos_IX);
MAP_state = state(MAP_state_IX);

%%
figure
clear hax;
hax(1)=subplot(311)
plot(exp.MUA.t, exp.MUA.zFR)
subplot(312)
plot(time,MAP_pos,'.')
