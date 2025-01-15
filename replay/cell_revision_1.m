%%
close all
clear
clc

%%
out_dir = 'E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\cell_revision';

%% load data
short_tunnel_data = load('L:\processed_data_structs\Shir_replay_events_15m.mat');
long_tunnel_data = load('L:\processed_data_structs\replay_events.mat');
short_tunnel_data.traj = load('E:\Tamir\work\PROJECTS\LargeScale\processed_data_structs\speed_traj_15m.mat');

%%
seqs_all = {};
% long tunnel
events = [long_tunnel_data.events_all_per_session{1,:}];
seqs = [events.seq_model];
seqs_all{1,1} = seqs;
events = [long_tunnel_data.events_all_per_session{2,:}];
seqs = [events.seq_model];
seqs_all{1,2} = seqs;
% short tunnel
seqs_all{2,1} = [short_tunnel_data.events_15m.sleep.seq_model];
seqs_all{2,2} = [short_tunnel_data.events_15m.rest.seq_model];

%%
epoch_types = long_tunnel_data.epoch_types;
tunnel_types = {'(200m or 130m)';'(15m)'};
epoch_type_clrs = {[.6 .1 .8],[.1 .8 .1]};
fig1 = figure(WindowState="maximized");
tiledlayout(2,4,'TileIndexing','columnmajor')
for ii_tunnel_type = 1:2
    for ii_epoch_type = 1:2
        seqs = seqs_all{ii_tunnel_type, ii_epoch_type};
        nexttile
        plot([seqs.duration],[seqs.compression],'.','Color',epoch_type_clrs{ii_epoch_type});
        xlabel('duration')
        ylabel('compression')
        title({epoch_types{ii_epoch_type};tunnel_types{ii_tunnel_type}})
        nexttile
        histogram([seqs.compression].*[seqs.duration])
        xlabel({'replay behavioral duration (s)';'[replay duration * compression]'});
        ylabel('counts')
    end
end
filename = fullfile(out_dir, "duration_vs_compression.pdf")
exportgraphics(fig1,filename)

%%
fig2 = figure;
axes(Units="centimeters",Position=[2 2 10 5 ])
for ii_tunnel_type = 1:2
    for ii_epoch_type = 1:2
        seqs = seqs_all{ii_tunnel_type, ii_epoch_type};
        line_styles = {'-';':'};
        histogram([seqs.compression].*[seqs.duration],'DisplayStyle','stairs','EdgeColor',epoch_type_clrs{ii_epoch_type},'Normalization','pdf','LineWidth',2,'LineStyle',line_styles{ii_tunnel_type},'DisplayName',""+epoch_types{ii_epoch_type}+" "+tunnel_types{ii_tunnel_type})
        hold on
        xlabel({'Replay behavioral duration (s)';'[Replay duration x compression ratio]'});
        ylabel('Probability density function')
    end
end
legend(Location="northeast",Box="off")
box off
filename = fullfile(out_dir, "replay_behavioral_duration.pdf");
% exportgraphics(fig2,filename)
saveas(fig2,filename,'pdf')

%%
params_opt = 11; % decoding opt
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
bats = unique(T.bat_num);
bats_clr_map = containers.Map(num2cell([bats;2295]),{'r','g','b','k','m','c',[0.8 0.8 0],[0.5 0.5 0.5]});
events_all_per_session = {};
speed_mean_per_session = nan(length(epoch_types),height(T));
speed_median_per_session = nan(length(epoch_types),height(T));
speed_std_per_session = nan(length(epoch_types),height(T));
distance_mean_per_session = nan(length(epoch_types),height(T));
flight_speed_all = [];
exp_arena_all = {};
for ii_exp = 1:height(T)
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details','flight');
    flight_speed_all(ii_exp,:) = [exp.flight.speed_traj.vel_median_median];
    exp_arena_all{ii_exp} = exp.details.recordingArena;
    for ii_epoch_type = 1:length(epoch_types)
        epoch_type = epoch_types{ii_epoch_type};
        event_type = 'posterior';
        [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
        [seqs, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
        events(~TF) = [];
        if isempty(events)
            continue;
        end
        seqs_map_IX = [seqs.state_direction];
        seqs_map_IX(seqs_map_IX == -1) = 2;
        [seqs.session_flight_speed] = disperse(flight_speed_all(ii_exp,seqs_map_IX));
        [events.seq_model] = disperse(seqs);
        [events.recordingArena] = deal(exp.details.recordingArena);
        event_struct = events(1,[]);
        [events(2:end).prev_event] = disperse(events(1:end-1));
        events(1).prev_event = event_struct;
        events_all_per_session{ii_epoch_type}{ii_exp} = events;
        speed_mean_per_session(ii_epoch_type,ii_exp) = mean([seqs.speed]);
        speed_median_per_session(ii_epoch_type,ii_exp) = median([seqs.speed]);
        speed_std_per_session(ii_epoch_type,ii_exp) = std([seqs.speed]);
        distance_mean_per_session(ii_epoch_type,ii_exp) = mean([seqs.distance]);
    end
end
flight_speed_all = mean(abs(flight_speed_all),2);
% events_all = {[events_all_per_session{1}{:}], [events_all_per_session{2}{:}]};
% seqs_all = {[events_all{1}.seq_model], [events_all{2}.seq_model]};
% gArena = {categorical({events_all{1}.recordingArena}), categorical({events_all{2}.recordingArena}) };

%% add 15m data
n = length(flight_speed_all);
for ii_exp = 1:size(short_tunnel_data.traj.speed_traj_all_exp,2)
    % flight speed
    avg_flight_speed = mean(abs([short_tunnel_data.traj.speed_traj_all_exp{1,ii_exp}.vel_median_median]));
    exp_ID = short_tunnel_data.traj.speed_traj_all_exp{2,ii_exp};
    bat_num = 2295;
    flight_speed_all(n+ii_exp) = avg_flight_speed;

    % replay distance
    for ii_epoch_type = 1:2
        epoch_type = epoch_types{ii_epoch_type};
        events = short_tunnel_data.events_15m.(epoch_type);
        seqs = [events.seq_model];
        TF = ismember({events.exp_ID},exp_ID);
        seqs(~TF)=[];
        distance_mean_per_session(ii_epoch_type,n+ii_exp)= mean([seqs.distance]);
    end

    % table info
    T.bat_num(n+ii_exp) = bat_num;
    T.exp_ID{n+ii_exp} = exp_ID;
end
bats_clr_map(bat_num) = 0.5*[1 1 1];
bats = bats_clr_map.keys;
bats = [bats{:}];

%%
fig3 = figure(Units="centimeters",Position=[5 5 20 20]);
tiledlayout("flow");
for ii_epoch_type = 1:2
    nexttile
    X = flight_speed_all;
    Y = distance_mean_per_session(ii_epoch_type,:)';
    plot(X,Y, 'o','Color',epoch_type_clrs{ii_epoch_type})
    xlabel('Avg. flight speed (m/s)')
    ylabel('Avg. replay distance (m)')
    hold on
    [r,pval] = corr(X,Y,'rows','pairwise');
    text(0.1,0.95, sprintf('r=%.2f, P=%.2g',r,pval),'Units','normalized')
    lsline
    nexttile
    hold on
    for ii_bat = 1:length(bats)
        bat_num = bats(ii_bat);
        IX = T.bat_num==bat_num;
        X_bat = prctile(X(IX),[50 25 75]);
        Y_bat = prctile(Y(IX),[50 25 75]);
%         X_bat(isnan(X_bat)) = [];
%         Y_bat(isnan(Y_bat)) = [];
        lw = 1.5;
        c = bats_clr_map(bat_num);
        plot(X_bat(1)*[1 1], Y_bat([2 3]),'LineWidth',lw,'Color',c,'DisplayName',""+bat_num);
        plot(X_bat([2 3]), Y_bat([1])*[1 1],'LineWidth',lw,'Color',c,'DisplayName',"");
    end
    legend(NumColumns=1,Location="bestoutside",Position=[0.9 0.8 0.1 0.15]);
end
filename = fullfile(out_dir, "replay_distance_vs_flight_speed.pdf");
exportgraphics(fig3,filename)

%%























%%
