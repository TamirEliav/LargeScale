%%
clear
clc

%% choose bats / sessions
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
clear exp_list
groupsummary(T,'bat_num')
if exist('bats_to_include','var')
    T = groupfilter(T,"bat_num",@(x)ismember(x,bats_to_include),'bat_num');
end
bats = unique(T.bat_num);

%% load data
exps = [];
for ii_exp = 1:height(T)
    % load exp data
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details');
    session_IX = contains(exp.details.session_names,'sleep','IgnoreCase',true);
    sleep_ti = exp.details.session_ts(session_IX,:);
    epoch_type = 'sleep';
    params_opt = 11;
    [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    % apply inclusion criteria 
    seqs = [events.seq_model];
    [~, TF] = decoding_apply_seq_inclusion_criteria(seqs);
    events(~TF)=[];
    exps(ii_exp).events = events;
    exps(ii_exp).sleep_ti = sleep_ti;
    exps(ii_exp).details = exp.details;
end

%%
fig = figure;
fig.WindowState = 'maximized';
for ii_exp = 1:length(exps)
    % prepare data 
    exp = exps(ii_exp);
    events = exp.events;
    % plot
    clf
    hold on
    for ii_event = 1:length(events)
        event = events(ii_event);
        x1 = event.seq_model.start_pos;
        x2 = event.seq_model.end_pos;
        y = event.peak_ts;
        plot([x1 x2],[y y])
    end
    pause
    sgtitle(exp.details.exp_ID,'interpreter','none');
end
% close(fig);





%%