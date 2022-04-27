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

%% load data (sleep)
events = {};
for ii_exp = 1:height(T)
    % load exp data
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details');
    epoch_type = 'sleep';
    params_opt = 11;
    [events_session, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    events{ii_exp} = events_session;
end
events = [events{:}];
seqs = [events.seq_model];
[seqs_sleep, TF] = decoding_apply_seq_inclusion_criteria(seqs);
clear events seqs

%% load data (rest)
events = {};
for ii_exp = 1:height(T)
    % load exp data
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details');
    epoch_type = 'rest';
    params_opt = 11;
    [events_session, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    events{ii_exp} = events_session;
end
events = [events{:}];
seqs = [events.seq_model];
[seqs_rest, TF] = decoding_apply_seq_inclusion_criteria(seqs);
clear events seqs

%%
fig=figure;
fig.Units='centimeters';
fig.Position=[10 5 15 8];
tiledlayout(2,3,'TileSpacing','loose')
% ----- sleep -------
nexttile
x = [seqs_sleep.duration]; h=histogram(x,'DisplayStyle','stairs','EdgeColor','k','BinWidth',0.05,'Normalization','count','LineWidth',2);
xlabel('Duration (s)'); ylabel('Count'); hax=gca; hax.YLim = [0 1.05*max(h.Values)];
h=text(-.8,0.5,{'Sleep';'replay'},'Units','normalized','HorizontalAlignment','center','FontWeight','bold','FontSize',14);
nexttile
x = [seqs_sleep.compression]; h=histogram(x,'DisplayStyle','stairs','EdgeColor','k','BinWidth',0.5,'Normalization','count','LineWidth',2);
xlabel('Compression-ratio'); ylabel('Count'); hax=gca; hax.YLim = [0 1.05*max(h.Values)];
nexttile
x = [seqs_sleep.distance]; h=histogram(x,'DisplayStyle','stairs','EdgeColor','k','BinWidth',0.5,'Normalization','count','LineWidth',2);
xlabel('Distance (m)'); ylabel('Count'); hax=gca; hax.YLim = [0 1.05*max(h.Values)];
% ----- rest -------
nexttile
x = [seqs_rest.duration]; h=histogram(x,'DisplayStyle','stairs','EdgeColor','k','BinWidth',0.05,'Normalization','count','LineWidth',2);
xlabel('Duration (s)'); ylabel('Count'); hax=gca; hax.YLim = [0 1.05*max(h.Values)];
text(-.8,0.5,{'Awake';'replay'},'Units','normalized','HorizontalAlignment','center','FontWeight','bold','FontSize',14);
nexttile
x = [seqs_rest.compression]; h=histogram(x,'DisplayStyle','stairs','EdgeColor','k','BinWidth',0.5,'Normalization','count','LineWidth',2);
xlabel('Compression-ratio'); ylabel('Count'); hax=gca; hax.YLim = [0 1.05*max(h.Values)];
nexttile
x = [seqs_rest.distance]; h=histogram(x,'DisplayStyle','stairs','EdgeColor','k','BinWidth',0.5,'Normalization','count','LineWidth',2);
xlabel('Distance (m)'); ylabel('Count'); hax=gca; hax.YLim = [0 1.05*max(h.Values)];
saveas(fig,'E:\Tamir\Conferences_Meetings_Visits\2022\04_Technion\figures\replay_sleep_rest_hist','pdf')






