function [events, params]= decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type)
arguments
    %% 
    exp_ID = 'b9861_d180526'
    epoch_type {mustBeMember(epoch_type,{'sleep','rest','flight'})} = 'sleep'
    params_opt = 11;
    event_type {mustBeMember(event_type,{'PE','posterior','ripples','MUA'})} = 'posterior'
end

dir_IN = 'F:\sequences\events_quantification';
filename = fullfile(dir_IN, sprintf('%s_events_%s_dec_prm_%d_%s',exp_ID,epoch_type,params_opt,event_type));
load(filename);
