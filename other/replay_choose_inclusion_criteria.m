%%
close all
clear cll
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
events_rest = {};
events_sleep = {};
for ii_exp = 1:height(T)
    % load exp data
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details');
    params_opt = 11;
    epoch_type = 'rest';
    events_rest{ii_exp} = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    epoch_type = 'sleep';
    events_sleep{ii_exp} = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
end
events_rest = [events_rest{:}];
events_sleep = [events_sleep{:}];

%% apply inclusion criteria 
seqs_rest = [events_rest.seq_model];
seqs_sleep = [events_sleep.seq_model];
% [seqs, TF] = decoding_apply_seq_inclusion_criteria(seqs);

%%
close all
fig=figure;
fig.WindowState = 'maximized';
tiledlayout('flow')
nexttile
plot([seqs_rest.score],[seqs_rest.distance],'.')
xlim([0 1])
ylim([0 30])
xlabel('Score')
ylabel('Distance (m)')
yline(3)
nexttile
plot([seqs_sleep.score],[seqs_sleep.distance],'.')
xlim([0 1])
ylim([0 30])
xlabel('Score')
ylabel('Distance (m)')
yline(3)
nexttile
x = [seqs_rest.score];
histogram(x([seqs_rest.distance]>3),linspace(0,1,50));
hold on
histogram(x([seqs_rest.distance]<3),linspace(0,1,50));
xlabel('Score')
title('rest')
legend("dist>3","dist<3")
nexttile
x = [seqs_sleep.score];
histogram(x([seqs_sleep.distance]>3),linspace(0,1,50));
histogram(x([seqs_sleep.distance]>3),linspace(0,1,50));
hold on
histogram(x([seqs_sleep.distance]<3),linspace(0,1,50));
xlabel('Score')
title('sleep')
legend("dist>3","dist<3")
