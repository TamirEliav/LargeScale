%% trial-by-trial speed correlation replay vs flight
clear
clc
close all

%%
out_dir = 'E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\cell_revision';

%%
load('L:\processed_data_structs\replay_events.mat');

%%
flight_speed_pos_range = [5 10];
clear exp_all
for ii_exp = 1:height(T)
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details','rest','flight');
    for ii_fe = 1:length(exp.flight.FE)
        fe = exp.flight.FE(ii_fe);
        switch fe.direction
            case 1
                pos_range = exp.rest.balls_loc(2) - flight_speed_pos_range;
            case -1
                pos_range = exp.rest.balls_loc(1) + flight_speed_pos_range;
        end
        pos_range  = sort(pos_range);
        TF = [fe.pos]>pos_range(1) & [fe.pos]<pos_range(2);
        vel = [fe.vel];
        exp.flight.FE(ii_fe).speed_before_landing = median(vel(TF));
    end
    exp_all(ii_exp) = exp;
end

%%
bats = unique(T.bat_num);
bats_colors = [1 0 0; 0 1 0; 0 0 1; 0 0 0; 1 0 1; 0 1 1; 0.8 0.8 0];

%%
fig1 = figure(WindowState="maximized");
htl = tiledlayout("flow");
x_per_sessions = {};
y_per_sessions = {};
x_norm_per_sessions = {};
y_norm_per_sessions = {};
c_all_sessions = {};
ii_epoch_type = find(epoch_types=="rest");
for ii_exp = 1:height(T)
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_all(ii_exp);
    events = events_all_per_session{ii_epoch_type,ii_exp};
    if isempty(events)
        continue;
    end
    FE_num = interp1([exp.flight.FE.start_ts],1:length(exp.flight.FE), [events.start_ts], "previous",'extrap');
    [events.FE_num] = disperse(FE_num);
    events(isnan(FE_num)) = [];
    if isempty(events)
        continue;
    end
    FE_speed_before_landing = [exp.flight.FE.speed_before_landing];
    speed_before_landing = FE_speed_before_landing([events.FE_num]);
    seqs = [events.seq_model];
    [seqs.speed_before_landing] = disperse(speed_before_landing);
    TakeLand_thr = 0.05;
    gTakeLand = classify_replay_landing_takeoff_other(seqs, TakeLand_thr);
    gPastFuture = categorical( ([events.rest_ball_num] == 1 & [seqs.state_direction] == 1) | ...
                               ([events.rest_ball_num] == 2 & [seqs.state_direction] == -1), ...
                               [false true],["Past","Future"])';
    gForRev = categorical([seqs.forward],[true false],["Forward","Reverse"])';
    TF = gTakeLand == 'Landing' & gPastFuture == 'Past';
    seqs(~TF)=[];
    if isempty(seqs)
        continue;
    end
    x = abs([seqs.speed_before_landing])';
    y = [seqs.speed]';
%     y = [seqs.compression]';
    x_per_sessions{ii_exp} = x;
    y_per_sessions{ii_exp} = y;
    avg_speed_before_landing = median(abs(FE_speed_before_landing),'omitnan');
    x_norm_per_sessions{ii_exp} = x ./ avg_speed_before_landing;
    y_norm_per_sessions{ii_exp} = y ./ avg_speed_before_landing;
    c = find(exp.details.batNum==bats);
    c_all_sessions{ii_exp} = repelem(c,length(x));

    if length(seqs)<10
        continue
    end
    nexttile
    hold on
    plot(x,y, 'o');
%     plot(x,y, 'o','Color',bats_colors(c,:));
%     hax=gca;
%     hax.Colormap = bats_colors;
    lsline
    [r,pval] = corr(x,y,'type','Pearson');
    text(0.05,.9,sprintf('r=%.2g',r),Units="normalized");
    text(0.05,.8,sprintf('P=%.2g',pval),Units="normalized");
    xlabel('Speed before landing (m/s)')
    ylabel('Replay speed (m/s)')
    title(exp.details.exp_ID,Interpreter="none");
end

nexttile
hold on
title('pooled')
x = cat(1,x_per_sessions{:});
y = cat(1,y_per_sessions{:});
plot(x,y, 'o')
lsline
[r,pval] = corr(x,y,'type','Pearson');
text(0.05,.9,sprintf('r=%.2g',r),Units="normalized");
text(0.05,.8,sprintf('P=%.2g',pval),Units="normalized");
xlabel('Speed before landing (m/s)')
ylabel('Replay speed (m/s)')

nexttile
hold on
title('pooled + normalized')
x = cat(1,x_norm_per_sessions{:});
y = cat(1,y_norm_per_sessions{:});
plot(x,y, 'o')
lsline
[r,pval] = corr(x,y,'type','Pearson');
text(0.05,.9,sprintf('r=%.2g',r),Units="normalized");
text(0.05,.8,sprintf('P=%.2g',pval),Units="normalized");
xlabel({'Speed before landing';'norm. to median speed before landing'});
ylabel({'Replay speed';'norm. to median speed before landing'})

fig1.Children.Title.String ='Replay speed vs flight speed, for past landing replays';

filename = fullfile(out_dir, 'trial-by-trial_speed_correlation_replay_vs_flight.pdf');
exportgraphics(fig1, filename)





%%
















%%
