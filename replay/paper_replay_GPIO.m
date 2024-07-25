%%
clear
clc

%%
exp_list_2bats = {
    'b2299_d191202',
    'b2299_d191203',
    'b2299_d191204',
    'b2299_d191205',
    'b2299_d191208',
    'b2299_d191209',
    'b2299_d191210',
    'b2299_d191213',
    };

%%
clear res_all
balls_locs = [];
for ii_exp = 1:length(exp_list_2bats)
    exp_ID = exp_list_2bats{ii_exp};
    exp = exp_load_data(exp_ID,'rest');
    balls_locs(ii_exp,:) = exp.rest.balls_loc;
    filename = fullfile('F:\sequences\decoded_figs\flight\conf_mat\', sprintf('%s_flight_decoding_opt_4.mat',exp_ID));
    res = load(filename);
    res_all(ii_exp) = res.res;
%     res.res.mean_err_prob
%     res.res.pos_err_median_prc
end
exp_list_mean_dec_acc = 100-[res_all.mean_err_prob]'.*100;
exp_list_median_pos_err_prc = [res_all.pos_err_median_prc]'.*100;
T_res = table(exp_list_2bats, exp_list_mean_dec_acc,exp_list_median_pos_err_prc)
% what distance from the balls do we 
fig_3_replay_inc_range = [25 115];
balls_locs - fig_3_replay_inc_range
mean(abs(balls_locs - fig_3_replay_inc_range),'all')

%%
x = [seqs_pooled.env_size];
thr = 150;
figure
histogram(x)
mean(x(x<thr))
mean(x(x>thr))
median(x(x<thr))
median(x(x>thr))
std(x(x<thr))
std(x(x>thr))

%%
L = [];
for ii_exp = 1:height(T)
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'rest');
    L(ii_exp) = diff(exp.rest.balls_loc);
end


%%
for ii_ex = 1:length(replay_examples)
    ex = replay_examples(ii_ex);
    addFieldsToWorkspace(ex);
    exp = exp_load_data(exp_ID,'ripples','details');
    exp.details.TT_loc{exp.ripples.stats.best_TT}
end


%% calc median rest pause duration
exp_list = decoding_get_inclusion_list();
rest_events_all = [];
for ii_exp = 1:length(exp_list)
    exp_ID = exp_list{ii_exp};
    exp = exp_load_data(exp_ID,'rest');
    events = exp.rest.events;
    rest_events_all = [rest_events_all events];
end
prctile([rest_events_all.duration],[25 50 75])
length([rest_events_all.duration])
figure
histogram([rest_events_all.duration])

%% calc median distance from flight end and landing ball
exp_list = decoding_get_inclusion_list();
dist_2_ball = [];
for ii_exp = 1:length(exp_list)
    exp_ID = exp_list{ii_exp};
    exp = exp_load_data(exp_ID,'rest','flight');
    FEs = exp.flight.FE;
    FEs([FEs.distance]<100) = [];
    FEs_end_pos = arrayfun(@(fe)fe.pos(end),FEs);
    ball_IX = interp1([1 -1],[2 1], [FEs.direction]);
    d = exp.rest.balls_loc(ball_IX) - FEs_end_pos;
%     d = abs(d);
    d = d.* [FEs.direction];
    dist_2_ball = [dist_2_ball d];
end

%% load all replays, and add state prob trace
epoch_types = {'sleep','rest'};
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
nExp = length(exp_list);
nEpochTypes = length(epoch_types);
events_all_per_session = cell(nEpochTypes, nExp);
% figure
%%
for ii_exp = 1:nExp
    for ii_epoch_type = 1:nEpochTypes
        
        %% load events
        exp_ID = exp_list{ii_exp};
        exp = exp_load_data(exp_ID,'details');
        epoch_type = epoch_types{ii_epoch_type};
        fprintf('%d: %s (%s)\n',ii_exp,exp_ID, epoch_type);
        params_opt = 11;
        event_type = 'posterior';
        [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
        [~, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
        events(~TF) = [];
        if isempty(events)
            continue;
        end
        %% load decoding
        decode = decoding_load_data(exp_ID,epoch_type,params_opt);
        %% extract prob trace
%         figure
%         title([exp_ID ' - ' epoch_type],'Interpreter','none');
%         hold on
        for ii_event = 1:length(events)
            event = events(ii_event);
            IX = event.start_IX:event.end_IX;
            prob = decode.posterior_state(event.state_num,IX);
%             plot(prob)
            events(ii_event).prob = prob;
        end
        if ~isfield(events,'rest_ball_num')
            [events.rest_ball_num] = disperse(nan(1,length(events)));
        end

        %%
        events_all_per_session{ii_epoch_type,ii_exp} = events;
        
    end
end

%% add number of ripples per replay
load("L:\processed_data_structs\replay_events.mat");
nEpochTypes = length(epoch_types);
parfor ii_exp = 1:length(exp_list)
    exp_ID = exp_list{ii_exp};
    exp = exp_load_data(exp_ID,'ripples');
    for ii_epoch_type = 1:nEpochTypes
        epoch_type = epoch_types{ii_epoch_type};
        fprintf('%d: %s (%s)\n',ii_exp,exp_ID, epoch_type);
        events = events_all_per_session{ii_epoch_type,ii_exp};
        if isempty(events)
            continue
        end
        seqs = [events.seq_model];
        replays_ti = [[seqs.start_ts];[seqs.end_ts]]';
        ripples_ts = [exp.ripples.events.peak_ts];
        [~, IX_per_ti, ~, ~] = get_data_in_ti(ripples_ts,replays_ti);
        num_ripples = cellfun(@length,IX_per_ti);
        [events.num_ripples] = disperse(num_ripples);
        events_all_per_session{ii_epoch_type,ii_exp} = events;
    end
end

%% save the events to mat file
res_dir = 'L:\processed_data_structs';
filename = 'replay_events2';
file_OUT = fullfile(res_dir, filename);
% save(file_OUT,'events_all_per_session','T','epoch_types','exp_list','params_opt','event_type');

%% calc percentage of replays that has a drop below 0.5 (=merged occurred) 
load("L:\processed_data_structs\replay_events.mat");
events = [events_all_per_session{:}];
seqs = [events.seq_model];
N = 100;
t = linspace(0,1,N);
M = zeros(length(events),N);
for ii_event = 1:length(events)
    event = events(ii_event);
    prob = event.prob;
    L = length(prob);
    M(ii_event,:) = interp1(1:L,prob,linspace(1,L,N));
end
sdf=mean(any(M<0.5,2));

%%
x = linspace(0,1,N);
figure
h=tiledlayout("flow",'TileSpacing','tight');
nexttile
hold on
clear hsel
hsel(1)=shadedErrorBar(x,M,{@mean,@nansem},'lineprops',{'r-'});
% nexttile
hsel(2)=shadedErrorBar(x,M,{@median,@nansem});
legend([hsel.mainLine],'mean','median')
xlabel('Replay time (norm)')
ylabel('State probability')
nexttile
prcs = [1 10:10:90 99];
prcs = flip(prcs);
cmap = cool(length(prcs));
% cmap = flip(cmap);
hax=gca;
hax.ColorOrder = cmap;
hold on
for ii_prc = 1:length(prcs)
    prc = prcs(ii_prc);
    y = prctile(M,prc,1);
    plot(x,y);
end
legend(prcs+"th percentile ")
xlabel('Replay time (norm)')
ylabel('State probability')
linkaxes(h.Children,'xy');
xlim([0 1])
ylim([0 1])

%% load all replays
epoch_types = {'sleep','rest'};
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
nExp = length(exp_list);
nEpochTypes = length(epoch_types);
events_all_per_session = cell(nEpochTypes, nExp);
for ii_exp = 1:nExp
    for ii_epoch_type = 1:nEpochTypes
        exp_ID = exp_list{ii_exp};
        exp = exp_load_data(exp_ID,'details');
        epoch_type = epoch_types{ii_epoch_type};
        fprintf('%d: %s (%s)\n',ii_exp,exp_ID, epoch_type);
        params_opt = 11;
        event_type = 'posterior';
        [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
        [~, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
        events(~TF) = [];
        if isempty(events)
            continue;
        end
        if ~isfield(events,'rest_ball_num')
            [events.rest_ball_num] = disperse(nan(1,length(events)));
        end
        events_all_per_session{ii_epoch_type,ii_exp} = events;
    end
end

%% calc percentage of replays > 20% of env size
events = [events_all_per_session{:}];
seqs = [events.seq_model];
thr = 0.2;
m = mean([seqs.distance_norm]>thr);
figure
% histogram([seqs.distance_norm],'Normalization','cdf',DisplayStyle='stairs')
ecdf([seqs.distance_norm]);
hax=gca;
h=hax.Children;
h.LineWidth=4;
xline(thr,'--')
yline(1-m,'--')
xlabel('Replay distance (norm)')
ylabel('CDF')


%% 
clear
clc
close all
exp_list = decoding_get_inclusion_list();
for ii_exp = 17:length(exp_list)
    exp_ID = exp_list{ii_exp};
    decoding_calc_session_seqs_spikes_corr(exp_ID);
    close all
end


%%
figure
clear panels
panels(1) = axes('Units','normalized','Position',[0.1 0.5 0.8 0.2])
panels(2) = axes('Units','normalized','Position',[0.1 0.5 0.8 0.2])
axes(panels(1))
imagesc(randn(100))
axes(panels(2))
plot(randn(1,100),randn(1,100),'ok')
panels(2).Color = 'None';

%%
