%%
clear
clc

%% replay - 82 days inclusion list
[exp_list,exp_t] = decoding_get_inclusion_list();
exp_t = exp_t(exp_list,:);

%%
exps = cellfun(@(exp_ID) exp_load_data(exp_ID,'rest','uturns') ,exp_list)   ;
uturns_all = [exps.uturns];

%%
histogram([uturns_all.pos])


%%
figure
hold on
for ii_exp = 1:length(exp_list)
    %%
    exp_ID = exp_list{ii_exp}
    exp = exp_load_data(exp_ID,'rest','uturns');
    epoch_type = 'rest';
    params_opt = 11;
    event_type = 'posterior';
    [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
    [~, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
    events(~TF) = [];
    if isempty(exp.uturns.ts) || isempty(events)
        continue;
    end
    seqs = [events.seq_model];

   %% super-naive
   x = exp.uturns.pos;
   replay_ts = [events.peak_ts];
   matches = arrayfun(@(u) find(replay_ts > u, 1, 'first'), exp.uturns.ts, 'UniformOutput', false);
   invalid_ix = cellfun(@isempty,matches);
   x(invalid_ix)=[];
   ix = [matches{:}];
   y = [seqs(ix).middle_pos];
   plot(x,y,'ok');

end




%%
