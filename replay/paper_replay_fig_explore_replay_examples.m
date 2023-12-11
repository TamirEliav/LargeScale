%%
clear
clc

%% load replay seq events from all included bats 
[exp_list,exp_T] = decoding_get_inclusion_list();
exp_T = exp_T(exp_list,:);

%%
params_opt = 11;
epoch_types = {'sleep','rest'};
epoch_type = 'sleep';
sleep_events = load_all_seq_data(exp_list,epoch_type,params_opt);
epoch_type = 'rest';
rest_events = load_all_seq_data(exp_list,epoch_type,params_opt);
sleep_seqs = [sleep_events.seq_model];
rest_seqs = [rest_events.seq_model];
[sleep_seqs.epoch_type] = disperse(repelem({'sleep'},length(sleep_seqs)));
[rest_seqs.epoch_type] = disperse(repelem({'rest'},length(rest_seqs)));
seqs_all = [sleep_seqs rest_seqs];
seqs_per_epoch = {sleep_seqs,rest_seqs};
events_per_epoch = {sleep_events,rest_events};

%% filter seqs
midrange = [20 80];
filters  = [
    struct('fn','score','range',[0.8 1],'type','val');
    struct('fn','duration','range',midrange,'type','prc');
    struct('fn','distance_norm','range',midrange,'type','prc');
    struct('fn','compression','range',midrange,'type','prc');
    struct('fn','confidence_HPD','range',[0 1],'type','val');
    struct('fn','middle_pos_norm','range',[0.1 0.9],'type','val');
    ]

filters  = [
    struct('fn','score','range',[0.8 1],'type','val');
    struct('fn','duration','range',[0.1 0.3],'type','val');
    struct('fn','distance_norm','range',[0.05 0.10],'type','val');
    struct('fn','compression','range',[5 15],'type','val');
    struct('fn','confidence_HPD','range',[0 1],'type','val');
    struct('fn','middle_pos_norm','range',[0.1 0.9],'type','val');
    ]

sum(filter_seqs(sleep_seqs,filters))
sum(filter_seqs(rest_seqs,filters))

%%
win_s = 0.5;
dir_OUT = 'L:\paper_replay\figures\Fig_replay_examples\auto';
for ii_epoch_type = 1:length(epoch_types)
    epoch_type = epoch_types{ii_epoch_type}
    seqs = seqs_per_epoch{ii_epoch_type};
    events = events_per_epoch{ii_epoch_type};
    [TF,seqs] = filter_seqs(seqs,filters);
    sum(TF)
    events = events(TF);
    [~, sorted_IX] = sort([seqs.score],'descend');
    seqs = seqs(sorted_IX);
    events = events(sorted_IX);
    res_dir = fullfile(dir_OUT,epoch_type+"_"+string(datetime('now','Format','yyyyMMdd_HHmmSS')));
    mkdir(res_dir);
    writetable(struct2table(filters), fullfile(res_dir,'filters.csv'));
    writetable(struct2table(events), fullfile(res_dir,'events.csv'));
    writetable(struct2table(seqs), fullfile(res_dir,'seqs.csv'));
    for ii_seq = 1:length(events)
        event = events(ii_seq);
        paper_replay_fig_single_replay_example( ...
            event.exp_ID, ...
            event.epoch_type, ...
            params_opt, ...
            event.num, ...
            win_s, ...
            'res_dir',res_dir, ...
            'filename_prefix',sprintf('%.3d_',ii_seq));
    end
end

%% run over manually-selected list and plot all examples
replay_examples_list_filename = "L:\Analysis\Code\inclusion_lists\replay_examples.xlsx";
ex_list = table2struct(readtable(replay_examples_list_filename));
win_s = 0.5;
dir_OUT = 'L:\paper_replay\figures\Fig_replay_examples\manual';
res_dir = fullfile(dir_OUT,string(datetime('now','Format','yyyyMMdd_HHmmSS')));
mkdir(res_dir);
for ii_ex = 1:length(ex_list)
    ii_ex
    ex = ex_list(ii_ex);
    addFieldsToWorkspace(ex)
%     if ~any(ismember(bat,[194 2289 9861]))
%         continue
%     end
    if ~strcmp(epoch_type,'sleep')
        continue
    end
    paper_replay_fig_single_replay_example( ...
            exp_ID, ...
            epoch_type, ...
            params_opt, ...
            event_num, ...
            win_s, ...
            'res_dir',res_dir, ...
            'filename_prefix',sprintf('%.3d_',ii_ex), ...
            'title_str_prefix',sprintf('%.3d_',ii_ex));
end

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





%% -------------------------------------------------------------------------
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
    [events.exp_ID] = disperse(repelem({exp_ID},length(events)));
    [events.epoch_type] = disperse(repelem({epoch_type},length(events)));
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


%% filter seqs
function [TF,seqs] = filter_seqs(seqs,filters)

%%
TF = true(size(seqs));
for ii_filt = 1:length(filters)
    filt = filters(ii_filt);
    x = [seqs.(filt.fn)];
    switch filt.type
        case 'val'
            range_vals = filt.range;
        case 'prc'
            range_vals = prctile(x,filt.range);
    end
    TF = TF & (x>range_vals(1) & x<range_vals(2));
end
seqs = seqs(TF);

end








%%



