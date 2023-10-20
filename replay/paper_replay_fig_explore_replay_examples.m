%%
clear
clc

%% load replay seq events from all included bats 
[exp_list,exp_T] = decoding_get_inclusion_list();
exp_T = exp_T(exp_list,:);

%%
params_opt = 11;
epoch_type = 'sleep';
sleep_events = load_all_seq_data(exp_list,epoch_type,params_opt);
epoch_type = 'rest';
rest_events = load_all_seq_data(exp_list,epoch_type,params_opt);
sleep_seqs = [sleep_events.seq_model];
rest_seqs = [rest_events.seq_model];
[sleep_seqs.epoch_type] = disperse(repelem({'sleep'},length(sleep_seqs)));
[rest_seqs.epoch_type] = disperse(repelem({'rest'},length(rest_seqs)));
seqs_all = [sleep_seqs rest_seqs]

%% gplots
[h,ax,bigax] = plot_seqs_features(seqs_all, {'duration','compression','score','distance_norm','middle_pos_norm','speed'}, 'epoch_type');

%%
figure
tiledlayout('flow')

nexttile
hold on
histogram([sleep_seqs.compression],'DisplayStyle','stairs','Normalization','pdf')
histogram([rest_seqs.compression],'DisplayStyle','stairs','Normalization','pdf')
xlabel('compression')
legend('sleep','rest')

nexttile
hold on
histogram([sleep_seqs.distance],'DisplayStyle','stairs','Normalization','pdf')
histogram([rest_seqs.distance],'DisplayStyle','stairs','Normalization','pdf')
xlabel('distance')
legend('sleep','rest')

nexttile
hold on
histogram([sleep_seqs.duration],'DisplayStyle','stairs','Normalization','pdf')
histogram([rest_seqs.duration],'DisplayStyle','stairs','Normalization','pdf')
xlabel('Duration (s)')
legend('sleep','rest')

nexttile
hold on
histogram([sleep_seqs.speed],'DisplayStyle','stairs','Normalization','pdf')
histogram([rest_seqs.speed],'DisplayStyle','stairs','Normalization','pdf')
xlabel('Speed (m/s)')
legend('sleep','rest')

%% stats
[~,pval_speed] = kstest2([sleep_seqs.speed],[rest_seqs.speed])
[~,pval_duration] = kstest2([sleep_seqs.duration],[rest_seqs.duration])
[~,pval_compression] = kstest2([sleep_seqs.compression],[rest_seqs.compression])
[~,pval_distance] = kstest2([sleep_seqs.distance],[rest_seqs.distance])

%%
% seqs_all = containers.Map({'sleep','rest'},{sleep_seqs,rest_seqs})

%%
function events_all = load_all_seq_data(exp_list,epoch_type,params_opt)
events_all=[];
for ii_exp = 1:length(exp_list)
    exp_ID = exp_list{ii_exp};
    exp = exp_load_data(exp_ID,'details','path');
    [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    if isempty(events)
        continue;
    end
    seqs = [events.seq_model];

    %% apply inclusion criteria
    [seqs, TF] = decoding_apply_seq_inclusion_criteria(seqs);
    if isempty(seqs)
        continue;
    end
    events(~TF) = [];
    events_all = [events_all events];
end
end

%% create gplotmatrix
function [h,ax,bigax] = plot_seqs_features(seqs_all, feature_list, grp_feat)
M = zeros(length(seqs_all),length(feature_list));
for ii_feat = 1:length(feature_list)
    feat_name = feature_list{ii_feat};
    M(:,ii_feat) = [seqs_all.(feat_name)];
end
G = categorical({seqs_all.(grp_feat)});
fig = figure;
fig.WindowState = 'maximized';
[h,ax,bigax] = gplotmatrix(M,[],G,[],[],[],[],'grpbars',feature_list);
end





