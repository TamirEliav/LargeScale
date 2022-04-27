function events = decoding_load_events(exp_ID, epoch_type, params_opt, event_type)
arguments
    %% 
    exp_ID = 'b9861_d180526'
    epoch_type {mustBeMember(epoch_type,{'sleep','rest','flight'})} = 'sleep'
    params_opt = 11;
    event_type {mustBeMember(event_type,{'PE','posterior','ripples','MUA'})} = 'posterior'
end

%% load exp data
exp = exp_load_data(exp_ID, 'rest', 'PE');

%%
switch event_type
    case 'PE'
        %%
        PE_to_use = exp.PE.thr;
        PE_ti = [PE_to_use.start_ts; PE_to_use.end_ts]';
        switch epoch_type
            case 'rest'
                rest_ti = exp.rest.ti;
                IX = any(PE_ti>shiftdim(rest_ti(:,1),-2) & PE_ti<shiftdim(rest_ti(:,2),-2), [2 3]);
                PE_to_use = PE_to_use(IX);
            case 'sleep'
                sleep_ti = exp_get_sessions_ti(exp_ID, 'Sleep1','Sleep2');
                sleep_ti(any(isnan(sleep_ti),2),:) = []; % remove nan in case of missing sessions
                IX = any(PE_ti>shiftdim(sleep_ti(:,1),-2) & PE_ti<shiftdim(sleep_ti(:,2),-2), [2 3]);
                PE_to_use = PE_to_use(IX);
            case 'flight'
                error('Not supported!')
        end
        events = PE_to_use;
    case 'posterior'
        %%
        dir_IN = 'F:\sequences\posterior_events';
        filename = fullfile(dir_IN, sprintf('%s_posterior_events_%s_dec_prm_%d',exp_ID,epoch_type,params_opt));
        load(filename,'events');
end

