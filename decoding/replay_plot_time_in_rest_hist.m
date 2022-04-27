%%
close all
clear cll
clc

%% params
grp_opt = 11; % opt 5 is the best option
TakeLand_thr = 0.05;
time_bin_width = 1;
max_time = 25;
t0_lw = 2;
normalize_by_TakeLand_fraction = true;

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
events = {};
for ii_exp = 1:height(T)
    % load exp data
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details');
%     epoch_type = 'sleep';
    epoch_type = 'rest';
    params_opt = 11;
    [events_session, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    events{ii_exp} = events_session;
end
T.nEvents = cellfun(@length,events)'; % note this is without filtering sequence by features!
sortrows( groupsummary(T,'bat_num',["median","mean","max","sum"],"nEvents"),"sum_nEvents", 'descend')
events = [events{:}];

%% apply inclusion criteria 
seqs = [events.seq_model];
[seqs, TF] = decoding_apply_seq_inclusion_criteria(seqs);
events(~TF)=[];

%% classify events
gEarlyLate = categorical([events.time_from_epoch_start]<1,[false,true],["Late","Early"])';
gForRev = categorical([seqs.forward],[true false],["Forward","Reverse"])';
gTakeLand = classify_replay_landing_takeoff_other(seqs, TakeLand_thr);
gPastFuture = categorical( ([events.rest_ball_num] == 1 & [seqs.state_direction] == 1) | ...
                           ([events.rest_ball_num] == 2 & [seqs.state_direction] == -1), ...
    [false true],["Past","Future"])';
gRestBall = categorical( [events.rest_ball_num], [1 2], ["Ball 1","Ball 2"])';

%% choose grouping
switch grp_opt
    case 1
        g = gPastFuture;
    case 2
        g = gForRev;
    case 3
        g = gTakeLand;
    case 4
        g = gPastFuture .* gForRev;
    case 5
        g = gTakeLand .* gPastFuture;
%         G = categories(g);
%         g = mergecats(g,G(contains(G,"Mid-air")),"Mid-air");
    case 6
        g = gTakeLand .* gForRev;
%         G = categories(g);
%         g = mergecats(g,G(contains(G,"Mid-air")),"Mid-air");
    case 7
        g = gTakeLand .* gPastFuture .* gForRev;
        G = categories(g);
        g = mergecats(g,G(contains(G,"Mid-air")),"Mid-air");
    case 8
        g = gRestBall;
    case 9
        g = gRestBall .* gTakeLand;
    case 10
        g = gRestBall .* gPastFuture;
    case 11
        g = gRestBall .* gForRev;
end
G = categories(g);
nG = length(G);

%% count rest time on the balls
exps = cellfun(@(exp_ID)(exp_load_data(exp_ID,'details','rest')),T.exp_ID);
rest = [exps.rest];
rest_durations = range(cat(1,rest.ti),2)*1e-6;
[N,EDGES] = histcounts(rest_durations,'BinWidth',time_bin_width,'BinLimits',[0 max(rest_durations)]);
N_rest = flip(cumsum(flip(N)));
CENTERS = edges2centers(EDGES);

%% count events during rest time
t1 = [events.time_from_epoch_start];
t2 = [events.time_to_epoch_end];
N_events_from_rest_start = zeros(nG,length(CENTERS));
N_events_from_rest_end = zeros(nG,length(CENTERS));
for ii_g = 1:nG
    N_events_from_rest_start(ii_g,:) = histcounts(t1(g==G{ii_g}),EDGES);
    N_events_from_rest_end(ii_g,:) = histcounts(t2(g==G{ii_g}),EDGES);
end

%% calc event rate
events_rate_from_rest_start = (N_events_from_rest_start./(N_rest.*time_bin_width));
events_rate_from_rest_end = (N_events_from_rest_end./(N_rest.*time_bin_width));
events_rate_ylabel_str = 'Events rate (events/s)';

%% in case of takeoff/landing classification, normalize the rate by the relative space
if normalize_by_TakeLand_fraction
    TakeLand_fraction = TakeLand_thr;
    Midair_fraction = 1-2*TakeLand_thr;
    TakeLand_G_IX = contains(G,{'Takeoff','Landing'});
    midair_G_IX = contains(G,'Mid-air');
    events_rate_from_rest_start(TakeLand_G_IX, :) = events_rate_from_rest_start(TakeLand_G_IX, :) ./ TakeLand_fraction;
    events_rate_from_rest_end(TakeLand_G_IX, :) = events_rate_from_rest_end(TakeLand_G_IX, :) ./ TakeLand_fraction;
    events_rate_from_rest_start(midair_G_IX, :) = events_rate_from_rest_start(midair_G_IX, :) ./ Midair_fraction;
    events_rate_from_rest_end(midair_G_IX, :) = events_rate_from_rest_end(midair_G_IX, :) ./ Midair_fraction;
    events_rate_ylabel_str = {'Events rate (events/s)';'adjusted to Takeoff/Landing/mid-air fractions'};
end

%%
fig=figure;
fig.WindowState = 'maximized';
tiledlayout(2,2,'TileSpacing','compact');
nexttile
stairs(CENTERS, events_rate_from_rest_start', 'LineWidth',2);
xline(0,'LineWidth',t0_lw,'Color','k')
hax=gca;
hax.XLim = [0 max_time];
hax.YLim(1) = 0;
xlabel('Time from rest start (s)')
ylabel(events_rate_ylabel_str);
nexttile
hold on
stairs(-CENTERS, events_rate_from_rest_end','LineWidth',2);
xline(0,'LineWidth',t0_lw,'Color','k')
hax=gca;
hax.XLim = [-max_time 0];
hax.YLim(1) = 0;
xlabel('Time from rest end (s)')
ylabel(events_rate_ylabel_str);
nexttile
hold on
stairs(CENTERS, N_events_from_rest_start','LineWidth',2);
xline(0,'LineWidth',t0_lw,'Color','k')
hax=gca;
hax.XLim = [0 max_time];
hax.YLim(1) = 0;
xlabel('Time from rest start (s)')
ylabel('Events count');
nexttile
hold on
stairs(-CENTERS, N_events_from_rest_end','LineWidth',2);
xline(0,'LineWidth',t0_lw,'Color','k')
hax=gca;
hax.XLim = [-max_time 0];
hax.YLim(1) = 0;
xlabel('Time from rest end (s)')
ylabel('Events count');
legend(G,'Location','bestoutside','Units','normalized','position',[0.25 0.01 0.5 0.05],'NumColumns',3);

%%
sgtitle({'Reverse replay immediately after landing but no forward replay before takeoff';...
         sprintf('bin width = %.2fs',time_bin_width)});
dir_out = 'F:\sequences\figures\timing_of_awake_replay_in_the_rest_epoch';
filename = sprintf('rest_replay_timing_hist_group_opt_%d_bin_%dms',grp_opt,round(time_bin_width*1e3));
if normalize_by_TakeLand_fraction
    filename = [filename '_norm_TakeLandMidair_fractions'];
end
filename = fullfile(dir_out,filename);

saveas(fig,filename,'jpg');

%%
disp('done!')





