function decoding_seq_quantify_add_info(exp_ID, epoch_type, params_opt , event_type)
arguments
    %% 
    exp_ID = 'b0184_d191208'
    epoch_type {mustBeMember(epoch_type,{'sleep','rest','flight'})} = 'rest'
    params_opt = 11;
    event_type {mustBeMember(event_type,{'PE','posterior','ripples','MUA'})} = 'posterior'
end

%% load data
exp = exp_load_data(exp_ID, 'details','path','rest','flight');
dir_IN_OUT = 'F:\sequences\events_quantification';
filename = fullfile(dir_IN_OUT, sprintf('%s_events_%s_dec_prm_%d_%s',exp_ID,epoch_type,params_opt,event_type));
load(filename);

%% epochs ts/name
switch epoch_type
    case 'rest'
        epochs_ti = exp.rest.ti;
        epochs_rest_ball_num = [exp.rest.events.ball_num];
    case 'sleep'
        sleep_session_TF = contains(exp.details.session_names, 'sleep','IgnoreCase',true);
        epochs_ti = exp.details.session_ts(sleep_session_TF,:);
        epochs_rest_ball_num = nan(1,size(epochs_ti,1));
    case 'flight'
        FE = [exp.flight.FE];
        FE([FE.distance]<100)=[];
        epochs_ti = [FE.start_ts; FE.end_ts]';
        epochs_rest_ball_num = nan(1,size(epochs_ti,1));
end
n_epochs = size(epochs_ti,1);
epochs_durations = range(epochs_ti,2) * 1e-6; % in seconds

%%
[events.epoch_num] = disperse(interp1(epochs_ti(:,1), 1:n_epochs, [events.peak_ts], 'previous','extrap'));
[events.epoch_duration] = disperse(epochs_durations([events.epoch_num]));
[events.time_from_epoch_start] = disperse(abs([events.peak_ts] - epochs_ti([events.epoch_num],1)') .* 1e-6); % in seconds
[events.time_to_epoch_end] =     disperse(abs([events.peak_ts] - epochs_ti([events.epoch_num],2)') .* 1e-6) ; % in seconds
[events.rest_ball_num] = disperse(epochs_rest_ball_num([events.epoch_num]));

%% save all events to mat file
save(filename, 'events','params');


end
