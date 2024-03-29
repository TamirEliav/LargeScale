function decoding_seq_quantify(decode, event_type)
arguments
    %% 
    decode
    event_type {mustBeMember(event_type,{'PE','posterior','ripples','MUA'})} = 'posterior'
end
%%
exp_ID = decode.exp_ID;
epoch_type = decode.epoch_type;
params_opt = decode.params_opt;

%% folders
dir_OUT = 'F:\sequences\events_quantification';
mkdir(dir_OUT);

%% load data
exp = exp_load_data(exp_ID, 'details','path','rest','ripples','MUA','PE','pos','flight','flight_6m');
events = decoding_load_events(exp_ID, epoch_type, params_opt, event_type);

%%
if isfield(exp,'flight_6m')
    FE = [exp.flight_6m.FE];
end
if isfield(exp,'flight')
    FE = [exp.flight.FE];
end

%% params & consts
radius = 5;
RW_speed = median(abs([FE.vel]));

%% calc
clear seqs_model seqs_bayes
for ii_event = 1:length(events)
    %%
    event = events(ii_event);
    IX = [event.start_IX:event.end_IX];
    prob_model = squeeze(decode.posterior(:,event.state_num,IX));
%     prob_bayes = squeeze(decode.likelihood(:,event.state_num,IX));
    time = decode.time(IX);
    pos = decode.pos;
    state = decode.state(event.state_num);
    res_model = decoding_seq_quantify_event(prob_model,'pos',pos,'time',time,'radius',radius,'state',state,'RW_speed',RW_speed,'time_units_factor',1e-6);
%     res_bayes = decoding_seq_quantify_event(prob_bayes,'pos',pos,'time',time,'radius',radius,'state',state,'RW_speed',RW_speed,'time_units_factor',1e-6);
    seqs_model(ii_event) = res_model;
%     seqs_bayes(ii_event) = res_bayes;
end
% add seq quantification to events struct
if length(events)~=0
    [events.seq_model] = disperse(seqs_model);
%     [events.seq_bayes] = disperse(seqs_bayes);
else
    [events.seq_model] = struct();
%     [events.seq_bayes] = struct();
end

%% save all events to mat file
params = struct();
params.decode = decode.params;
params.radius = radius;
params.RW_speed = RW_speed;
filename = fullfile(dir_OUT, sprintf('%s_events_%s_dec_prm_%d_%s',exp_ID,epoch_type,params_opt,event_type));
save(filename, 'events','params');


end








%%
