function PRE_run_exp(exp_ID)

%%
% clear all
% exp_ID = 'b0148_d170711';
% exp_ID = 'b0079_d160926';

%% read exp info
exp=exp_load_data(exp_ID,'details','path');

%% Position (bsp related)
bsp_extract_data(exp.path.bsp)

%% Neural (Spikes+LFP)
Nlg2Nlx(exp.path.raw) % TODO: insert params from here, rather than change them inside Nlg2Nlx function, also consider to remove some params (like ref channel...). this code should ONLY be a reader/parser code!
PRE_filter_CSCs(exp_ID);
PRE_detect_spikes(exp_ID);
% TODO: add some analysis AFTER spike sorting (separation matrices/FR/AC/xcorr/ISI/...)

%% Audio
nlg2nlx_audio(exp.path.audio);

%% sync position/neural/audio
PRE_sync_bsp_to_nlg(exp.path.bsp, exp.path.nlx, exp.path.sync);
PRE_sync_audio_to_nlg(exp.path.audio, exp.path.nlx, exp.path.sync);

%% Position pre-process
POS_pre_process(exp_ID);
% Motion_pre_process(exp_ID); % TODO: calc flights ephocs - consider putting this somewhere else in the code, later...

%% flight ryhthm
% PRE_calc_flight_rhythm(exp_ID);

%% LFP
LFP_filter_theta(exp_ID);

%% show event to fill the excel record
PRE_event_disp(exp_ID);

%% spikes pre-process
% TODO: calc basic spike train metrics (like AC,xcorr,ISI,FR in each session, bursty, clustering 
% TODO: assign position data to the spikes


end

