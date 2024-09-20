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
T = table();
T.bat_num = exp_t.bat_num;
T.uturns_count = arrayfun(@(x)length(x.pos),uturns_all)';

%%
figure
subplot(211)
histogram([uturns_all.pos],20)
subplot(212)
histogram([uturns_all.pos_norm],20)


%% arrange data
uturn_pos_all = [];
uturn_pos_norm_all = [];
replay_pos_all = [];
replay_pos_norm_all = [];
replay_gTakeLand_all = [];
time_from_uturn_all = [];
bat_num_all = [];
for ii_exp = 1:length(exp_list)
    %%
    exp_ID = exp_list{ii_exp};
    exp = exp_load_data(exp_ID,'details','rest','uturns');
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
    
    %% 
    uturn_pos = exp.uturns.pos;
    uturn_pos_norm = exp.uturns.pos_norm;
    replay_ts = [events.peak_ts];
    matches = arrayfun(@(u) find(replay_ts > u, 1, 'first'), exp.uturns.ts, 'UniformOutput', false);
    invalid_ix = cellfun(@isempty,matches);
    uturn_pos(invalid_ix)=[];
    uturn_pos_norm(invalid_ix)=[];
    ix = [matches{:}];
    replay_pos = [seqs(ix).middle_pos];
    replay_pos_norm = [seqs(ix).middle_pos_norm];
    replay_gTakeLand = classify_replay_landing_takeoff_other(seqs(ix));
    time_from_uturn = [events(ix).peak_ts] - exp.uturns.ts(~invalid_ix);
    time_from_uturn = time_from_uturn*1e-6;

    uturn_pos_all = [uturn_pos_all; uturn_pos'];
    uturn_pos_norm_all = [uturn_pos_norm_all; uturn_pos_norm'];
    replay_pos_all = [replay_pos_all; replay_pos'];
    replay_pos_norm_all = [replay_pos_norm_all; replay_pos_norm'];
    replay_gTakeLand_all = [replay_gTakeLand_all; replay_gTakeLand];
    time_from_uturn_all = [time_from_uturn_all; time_from_uturn'];
    bat_num_all = [bat_num_all; exp.details.batNum.*ones(size(uturn_pos))'];
end


%%
bat_index_all = findgroups(bat_num_all);
data_t = table(uturn_pos_all, uturn_pos_norm_all, replay_pos_all, replay_pos_norm_all, time_from_uturn_all, bat_num_all, bat_index_all);

%%
figure
gplotmatrix(table2array(data_t),[],bat_num_all);

%%
opt = 5;
TF = true(size(uturn_pos_norm_all));
TF = TF & uturn_pos_norm_all>0.05 & uturn_pos_norm_all<0.95;
TF = TF & replay_gTakeLand_all=="Mid-air";
max_timediff = 30;
TF = TF & time_from_uturn_all<max_timediff;
X = uturn_pos_norm_all(TF);
Y = replay_pos_norm_all(TF);

fig=figure;
fig.Units = 'centimeters';
fig.OuterPosition = [5 5 25 20];
axes('position', [0.1 0.1 0.4 0.7],'Units','normalized');
plot(X, Y, '.')
xlabel('uturn position (norm)')
ylabel('replay position (norm)')
axis equal
axis tight
h=refline(1,0);
h.Color = [1 1 1]*0.5;

axes('position', [0.55 0.1 0.4 0.7],'Units','normalized');
histogram2(X,Y,20,'DisplayStyle','tile')
axis equal
axis tight

axes('position', [0.25 0.8 0. 0.15],'Units','normalized');
axis off
[r,r_pval] = corr(X,Y,'type','Pearson');
[rho,rho_pval] = corr(X,Y,'type','Spearman');
% lm = fitlm(X,Y)
str = {
    sprintf('r = %.2g, p = %.2g',r,r_pval);
    sprintf('rho = %.2g, p = %.2g',rho,rho_pval);
    '';
    'midair replay';
    'uturn position 5-95%',
    sprintf('uturn-replay time diff < %d s',max_timediff)};
text(.05,.5,str,'units','normalized','HorizontalAlignment','center')

dir_out = 'L:\paper_replay\figures\uturns';
filename = "uturns_scatter_opt_"+ opt;
file_out = fullfile(dir_out,filename);
saveas(fig,file_out,'jpg');



%%








%%

