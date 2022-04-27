%%
exps_IDs = exp_t.exp_ID;

%%
exps = cellfun(@(exp_ID)( exp_load_data(exp_ID,'details','PE') ), exps_IDs);
PEs = {exps.PE};
PE_all = [exps.PE];
nPE = cellfun(@length, PEs);
mean_FR  = cellfun(@(PE)(mean([PE.peak_FR])), PEs);
mean_zFR  = cellfun(@(PE)(mean([PE.peak_zFR])), PEs);
order_feature = nPE .* mean_zFR;
[~,ranks]  = ismember(order_feature,flip(unique(order_feature)));
% [~,sort_IX] = sort(order_feature,'descend');
[exps.nPE] = disperse(nPE);
[exps.mean_zFR] = disperse(mean_zFR);
[exps.order_feature] = disperse(order_feature);
[exps.rank] = disperse(ranks);

%%
figure
tiledlayout('flow')
nexttile
h=scatter(nPE,mean_FR,20,order_feature,'filled');
h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('ranks',ranks);
h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('feature value',order_feature);
% lsline
xlabel('No. of events')
ylabel('Mean firing rate (Hz)')
colorbar

nexttile
h=scatter(nPE,mean_zFR,20,order_feature,'filled');
h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('ranks',ranks);
h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('feature value',order_feature);
% lsline
xlabel('No. of events')
ylabel('Mean firing rate (z)')
colorbar

%% sort exps
% exps=exps(sort_IX);

T1 = table( arrayfun(@(exp)(exp.details.exp_ID), exps,'UniformOutput',0),...
            arrayfun(@(exp)(exp.details.batNum), exps,'UniformOutput',1),...
            datetime(arrayfun(@(exp)(exp.details.date), exps,'UniformOutput',1),'Format','yyMMdd'),...
            arrayfun(@(exp)(exp.nPE), exps,'UniformOutput',1),...
            arrayfun(@(exp)(exp.mean_zFR), exps,'UniformOutput',1),...
            arrayfun(@(exp)(exp.order_feature), exps,'UniformOutput',1),...
            arrayfun(@(exp)(exp.rank), exps,'UniformOutput',1),...
            'VariableNames',{'exp_ID','batNum','date','nPE','mean_zFR','order_feature','rank'});
T1.Properties.RowNames = string(datestr(T1.date,'yymmdd'));

%%
cells_table_filename = "L:\Analysis\Code\inclusion_lists\cells_summary.xlsx";
cells_t = readtable(cells_table_filename);
T2 = groupcounts(cells_t,'date');
T2.Properties.VariableNames(2) = {'nCells'};
T2.date = datetime(T2.date,'Format','yyMMdd');
T2.Properties.RowNames = string(datestr(T2.date,'yymmdd'));

%%
T = outerjoin(T2,T1,'Keys','Row');
T = sortrows(T,'nCells','descend');
clc
T

%% creating replay examples
decoding_plot_event('b9861_d180527', 'sleep', 11, 824, 0.5)
decoding_plot_event('b9861_d180527', 'sleep', 11, 681, 0.5)
decoding_plot_event('b9861_d180527', 'sleep', 11, 503, 0.5)
decoding_plot_event('b9861_d180526', 'sleep', 11, 27,  0.5)
decoding_plot_event('b9861_d180526', 'sleep', 11, 42,  0.5)
decoding_plot_event('b9861_d180526', 'sleep', 11, 113, 0.5)
decoding_plot_event('b9861_d180526', 'sleep', 11, 113, 1)
decoding_plot_event('b9861_d180526', 'sleep', 11, 113, 1.5)
decoding_plot_event('b9861_d180526', 'sleep', 11, 113, 2)
decoding_plot_event('b9861_d180526', 'sleep', 13, 113, 1)
decoding_plot_event('b9861_d180526', 'sleep', 11, 132, 0.5)
decoding_plot_event('b9861_d180526', 'sleep', 11, 132, 1)
decoding_plot_event('b0184_d191201', 'sleep', 11, 136, 0.5)
decoding_plot_event('b0184_d191130', 'sleep', 11, 53,  0.5)
decoding_plot_event('b0184_d191130', 'sleep', 11, 333,  0.5)
decoding_plot_event('b0184_d191130', 'sleep', 11, 350,  0.5)

%%
decoding_plot_event('b9861_d180526', 'sleep', 11, 'posterior', 49, [-.45 1.12]);

%%
decoding_plot_event('b9861_d180526', 'sleep', 14, 113, 1);
decoding_plot_event('b9861_d180527', 'sleep', 11, 824, 0.42);
decoding_plot_event('b0184_d191130', 'sleep', 11, 53,  [-.32 .11]);

%%
% decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',49);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',98);

%%
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',15);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',22);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',23);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',25);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',44);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',49);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',50);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',51);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',52);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',55);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',57);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',60);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',70);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',75);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',84);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',90);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',98);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',105);

%%
clc
X1 = [3 1 10 8 16];
u = [1 1 1 1 1];
% X2 = conv(X1,u,'full');
X2 = imfilter(X1,u,'symmetric','conv','same');
% X2 = imfilter(X1,u,'replicate','conv','same');
disp(X1)
disp(X2)


%%
clear; clc
exp_ID = 'b0184_d191130';
exp = exp_load_data(exp_ID, 'details', 'path');
dir_IN = exp.path.spikes_raw;
limits_ts = [];
TT = exp.details.refCh(1);
ch = exp.details.refCh(2);
csc_file = sprintf('spikes_%s_TT%d_ch%d.ncs',exp_ID,TT,ch);
csc_file = fullfile(dir_IN,csc_file);
csc_ref_ch = Nlx_csc_read(csc_file, limits_ts);
data = zeros([size(exp.details.activeChannels),length(csc_ref_ch)]);
%%
tic
for TT = 1:exp.details.numTT
    parfor ch = 1:4
        fprintf('TT%d_ch%d\n',TT,ch);
        % load ch data
        csc_file = sprintf('spikes_%s_TT%d_ch%d.ncs',exp_ID,TT,ch);
        csc_file = fullfile(dir_IN,csc_file);
        if ~exist(csc_file,'file')
            continue;
        end
        data(TT,ch,:) = Nlx_csc_read(csc_file,limits_ts);
    end
end
toc

%%
sdf=median(data,[1 2]);


%% select best TT for dislaying ripple traces
events_all = ripples.all;
IX = [events_all.peak_IX];
IX = IX+[-10:10]';
pripple_TT_at_peaks = pripple_TT(IX,:);
mean(pripple_TT_at_peaks)
n_contrib_events = zeros(1,nTT);
contrib_events = false(length(events_all),nTT);
for TT=1:nTT
    events_TT = ripples.by_TT{TT};
    if isempty(events_TT)
        continue;
    end
    t1 = [events_all.peak_ts];
    t2 = [events_TT.peak_ts];
    tdiff = abs(t1-t2');
    thr = 100; % usec
    TF = tdiff < thr;
    fprintf('TT %d: %d\n',TT,sum(any(TF)))
    contrib_events(:,TT) = any(TF);
    n_contrib_events(TT) = sum(any(TF));
end
figure
subplot(131);bar(1:nTT,n_contrib_events); ylabel('num events')
subplot(132);bar(1:nTT,mean(pripple_TT_at_peaks)); ylabel('pripple')
subplot(133);bar(1:nTT,n_contrib_events.*mean(pripple_TT_at_peaks)); ylabel('product')


%%
IX = [events_all.peak_IX];
IX = IX+[-200:200]';
% signal = LFP2;
% signal = sum(ripple,3);
signal = pripple_TT;
sdf = signal(IX,:);
sdf = reshape(sdf,size(IX,1),size(IX,2),nTT);

%%
close all
for TT=TT_to_use
    figure
    sdf2=squeeze(sdf(:,:,TT));
    sdf3=sdf2;
%     sdf3 = sdf3+[1:nEvents].*25;
%     plot(sdf3)
    imagesc(sdf3')
    title("TT"+TT)
end

%%
% sdf4= squeeze(sum(sdf,1));
sdf4 = squeeze(max(sdf,[],1));
% sdf4(isnan(sdf4))=0;
sdf4 = sdf4(:,TT_to_use);
[IDX C] = kmeans(sdf4,2);

[coeff,score,latent,tsquared,explained,mu] = pca(sdf4);
[~,sorted_IX] = sort(score(:,2));

% figure
% scatter(score(:,1),score(:,2),5,IDX,'filled')

figure
imagesc(coeff(:,[1 2]))
yticklabels(TT_to_use)
colormap bone
colorbar

figure
imagesc(contrib_events(sorted_IX,TT_to_use))

figure
subplot(121)
imagesc(squeeze(sdf(:,sorted_IX,1))')
% imagesc(squeeze(sdf(:,:,1))')
subplot(122)
imagesc(squeeze(sdf(:,sorted_IX,12))')
% imagesc(squeeze(sdf(:,:,12))')

%% 09/01/2022 - copy flight decoding results per inclusion criteria
% inc_type = 'manual';
inc_type = 'auto';
inc_list_filename = 'F:\sequences\inclusion\exp_inc_list.xlsx';
T = readtable(inc_list_filename, 'ReadRowNames',1,'ReadVariableNames',1);
dir_IN = 'F:\sequences\decoded_figs\flight\conf_mat';
main_dir_OUT = 'F:\sequences\inclusion\';
main_dir_OUT = fullfile(main_dir_OUT,inc_type);
switch inc_type
    case 'manual'
        [group, id] = findgroups(T.included_manual);
    case 'auto'
        [group, id] = findgroups(T.included_auto);
end
for g = 1:length(id)
    IX=group==g;
    tmpl_list = T.exp_ID(IX);
    ext = 'jpg';
    dir_OUT = fullfile(main_dir_OUT,"grp_"+id(g));
    mkdir(dir_OUT);
    util_copy_files_by_template(dir_IN, dir_OUT, tmpl_list, ext);
end
%%
T2 = T(T.included_auto,:);
figure
tiledlayout('flow')
nexttile
histogram(T2.pos_err_median_m,20)
text(0.75,0.9,"n="+length(x),'units','normalized');
xlabel('Median position error (m)')
ylabel('Counts')
nexttile
histogram(T2.pos_err_median_prc,20)
text(0.75,0.9,"n="+length(x),'units','normalized');
xlabel('Median position error (%)')
ylabel('Counts')
nexttile
histogram(T2.mean_err_prob,20)
text(0.75,0.9,"n="+length(x),'units','normalized');
xlabel('Error probability')
ylabel('Counts')
nexttile
histogram(T2.err_prob_by_predicted_max,20)
text(0.75,0.9,"n="+length(x),'units','normalized');
xlabel('Max predicted error prob.')
ylabel('Counts')
saveas(gcf,'F:\sequences\inclusion\inc_features_hist','jpg');

%% check replay compression value of examples in my thesis figure
compressions = [];
for ii_example = 1:size(panel_A,1)
    %% load data
    exp_ID = examples_list(ii_example).exp_ID;
    epoch_type = examples_list(ii_example).epoch_type;
    event_num = examples_list(ii_example).event_num;
%     decode = decoding_load_data(exp_ID, 'sleep', 11);
%     exp = exp_load_data(exp_ID,'details','path','MUA','ripples');
%     TT = exp.ripples.stats.best_TT;
%     [LFP.signal, LFP.ts, LFP.fs, LFP.params] = LFP_load(exp_ID,TT,'band','ripple');
%     LFP.avg_signal = nanmean(LFP.signal,[2 3]);
    events = decoding_load_events_quantification(exp_ID,epoch_type,11,"posterior");
    event = events(event_num);
    compressions(ii_example) = event.seq_model.compression;
end

%% 26/02/2022 - check sleep session durations
session_names = [];
session_durations = [];
for ii_exp = 1:height(T)
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details');
    session_names = [session_names; string(exp.details.session_names)];
    ti = exp_get_sessions_ti(exp_ID,'Sleep1','Sleep2');
    session_durations = [session_durations  diff(ti,1,2)*1e-6/60];
end
session_durations = session_durations';
bats = T.bat_num;
x = session_durations(:);
g = [bats; bats];
violinplot(x,g);
xlabel('bat');
ylabel('Sleep session duration (s)');
min(session_durations(:))
max(session_durations(:))
mean(session_durations(:))
median(session_durations(:))
prctile(session_durations(:),[25 75])

%% 19/03/2022 - how many TTs were in CA1 or CA3 per session
[exp_list, T] = decoding_get_inclusion_list();
T(~T.included_auto,:)=[];
for ii_exp = 1:height(T)
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID, 'details');
    T.num_tt_ca1(ii_exp) = sum(contains(exp.details.TT_loc,'CA1'));
    T.num_tt_ca3(ii_exp) = sum(contains(exp.details.TT_loc,'CA3'));
end

%% rest duration hist per bat
[exp_list, T] = decoding_get_inclusion_list();
T(~T.included_auto,:)=[];
exps = cellfun(@(exp_ID)(exp_load_data(exp_ID,'details','rest')),T.exp_ID);
rest = [exps.rest];
details = [exps.details];
for ii = 1:length(rest)
    [rest(ii).events.bat_num] = deal(details(ii).batNum);
end
events = [rest.events];
[~,sorted_IX] = sort([events.duration],'ascend');
events = events(sorted_IX);

%%
figure
subplot(211)
hold on
n = 200;
edges = linspace(0,120,n);
durations = [events.duration];
histogram(durations,edges,'Normalization','count');
title('all bats')
xlabel('Rest on balls duration (s)')
ylabel('Counts')
subplot(212)
hold on
[G,ID] = findgroups([events.bat_num]);
for ii_bat=1:length(ID)
    IX = G==ii_bat;
    histogram(durations(IX),edges,'Normalization','count','DisplayStyle','stairs','LineWidth',2)
end
legend("bat "+ID);
title('Per bat')
xlabel('Rest on balls duration (s)')
ylabel('Counts')
linkaxes(findobj(gcf, 'type', 'axes'), 'x')

%%
figure
hold on
duration_limits = [1 10];
% duration_limits = [1.2 1.8];
nClrPoints = 256;
duration_clr_points = linspace(duration_limits(1),duration_limits(2),nClrPoints);
cmap = parula(nClrPoints);
% duration_clr_points = [duration_limits(1) mean(duration_limits) duration_limits(2)];
% cmap = [1 0 0;0 1 0; 0 0 1];
% cmap = parula;
TF = [events.duration]>duration_limits(1) & [events.duration]<duration_limits(2);
% arrayfun(@(x,clr)plot(x.pos_original-median(x.pos_original), 'Color',interp1([0 .5 1],[1 0 0;0 1 0; 0 0 1],clr)), events(TF), linspace(0,1,sum(TF)));
arrayfun(@(x)plot(x.pos_original-median(x.pos_original), 'Color',...
    interp1(duration_clr_points,cmap,x.duration)), events(TF));
rescale_plot_data('x',[1/100 0]);
xlabel('Time (s)')
ylabel('{\Delta}Position from median position (m)')
colormap( interp1(duration_clr_points,cmap,duration_clr_points))
hcb=colorbar ;
hcb.Ticks=[0 1];
hcb.TickLabels = duration_limits;
hcb.Label.String = 'Duration (s)'

%%
bubblechart([events.mean_prob]-[seqs.score],[seqs.duration],[seqs.score],[seqs.score])
bubblesize([1 7])

%% 07/04/2022 - test rest detection algorithm (mainly the median position filtering)
clear
clc
fs = 100;
dt = 1/fs;
T = 60;
t = 0:dt:T;
for rest_win_s = [0.1:0.1:2 2.5 3 3.5 4]

    x = [linspace(0,200,fs*20) 200*ones(1,round(fs*rest_win_s)) linspace(200,0,fs*20)];
    x(length(t))=0;
    noise_std = .1;
    rng(0);
    x = x+randn(size(x)).*noise_std;
    
    win_s = 1;
    win = round(1 * fs);
%     win2 = round(.1 * fs);
    xfilt = smoothdata(x,2, "movmedian",win);
%     csaps_p = 1e-5;
%     pos_csaps_pp = csaps(1:length(t), xfilt', csaps_p);
%     xfilt = fnval(pos_csaps_pp, 1:length(t) );
%     v = fnval( fnder(pos_csaps_pp), 1:length(t)) .* fs;
    v = [0 diff(xfilt)].*fs;
%     v = smoothdata(v,"movmedian",win2);
    xdev = x-xfilt;

    speed_thr_1 = 1;
    speed_thr_2 = .1;
    [events, g, opts] = detect_events(100-abs(v),"thr",100-speed_thr_1,"thr_edges",100-speed_thr_1,"min_width",.25*win);
    IX = events(1).start_IX:events(1).end_IX;
    IX1 = find(abs(v(IX))<speed_thr_2,1,'first');
    IX2 = find(abs(v(IX))<speed_thr_2,1,'last');
    IX = IX(IX1:IX2);
    
    close all
    fig=figure;
    fig.WindowState = 'maximized';
    tiledlayout(3,1,'TileSpacing','compact')
    ax=[];
    
    ax(end+1) = nexttile;
    hold on
    plot(t,x,'x')
    plot(t,xfilt,'-g','LineWidth',2)
    plot(t(IX),x(IX),'.r')
    legend("original","Filtered","detected")
    ylabel('Position (m)');
    
    ax(end+1) = nexttile;
    hold on
    plot(t,xdev,'.')
    plot(t(IX),xdev(IX),'.r')
    ylabel('Original position - Filtered position (m)');
    
    ax(end+1) = nexttile;
    hold on
    plot(t,abs(v),'.-')
    plot(t(IX),abs(v(IX)),'.r')
    xlabel('time (s)');
    ylabel('Speed (m/s)');
    
    linkaxes(ax,'x')
    xlim([19 25])
    sgtitle({"real rest duration: " + rest_win_s + "s";...
            sprintf('detected rest duration: %.2f',range(IX)/fs);...
            sprintf('position median filter size: %.1fs',win_s)})
    dir_out = 'L:\Analysis\Results\exp\rest\tests';
    mkdir(dir_out);
    filename = fullfile(dir_out, sprintf('rest_detection_median_smothing_duration_%.2fs', rest_win_s) );
    filename = strrep(filename,'.','_');
    saveas(fig,filename,'fig');
    saveas(fig,filename,'jpg');

end

%% 24/04/2022 - time in rest histogram (forward vs reverse) - like in Diba 2007
for time_bin_width = 1%0.1:0.1:1
t1=[events.time_from_epoch_start];
t2=[events.time_to_epoch_end];
exps = cellfun(@(exp_ID)(exp_load_data(exp_ID,'details','rest')),T.exp_ID);
rest = [exps.rest];
rest_durations = range(cat(1,rest.ti),2)*1e-6;
[N,EDGES] = histcounts(rest_durations,'BinWidth',time_bin_width,'BinLimits',[0 max(rest_durations)]);
N_rest = flip(cumsum(flip(N)));
N1_forward = histcounts(t1([seqs.forward]),EDGES);
N1_reverse = histcounts(t1(~[seqs.forward]),EDGES);
N2_forward = histcounts(t2([seqs.forward]),EDGES);
N2_reverse = histcounts(t2(~[seqs.forward]),EDGES);
CENTERS = edges2centers(EDGES);

max_time = 25;
t0_lw = 2;
close all
fig=figure;
fig.WindowState = 'maximized';
subplot(221)
hold on
stairs(CENTERS, N1_forward./(N_rest.*time_bin_width),'LineWidth',2);
stairs(CENTERS, N1_reverse./(N_rest.*time_bin_width),'LineWidth',2);
xline(0,'LineWidth',t0_lw,'Color','k')
hax=gca;
hax.XLim = [0 max_time];
hax.YLim(1) = 0;
xlabel('Time from rest start (s)')
ylabel('Events rate (events/s)');
legend("forward","Reverse",'Location','northeast')
subplot(222)
hold on
stairs(-CENTERS, N2_forward./(N_rest.*time_bin_width),'LineWidth',2);
stairs(-CENTERS, N2_reverse./(N_rest.*time_bin_width),'LineWidth',2);
xline(0,'LineWidth',t0_lw,'Color','k')
hax=gca;
hax.XLim = [-max_time 0];
hax.YLim(1) = 0;
xlabel('Time from rest end (s)')
ylabel('Events rate (events/s)');
legend("forward","Reverse",'Location','northwest')
subplot(223)
hold on
stairs(CENTERS, N1_forward,'LineWidth',2);
stairs(CENTERS, N1_reverse,'LineWidth',2);
xline(0,'LineWidth',t0_lw,'Color','k')
hax=gca;
hax.XLim = [0 max_time];
hax.YLim(1) = 0;
xlabel('Time from rest start (s)')
ylabel('Events count');
legend("forward","Reverse",'Location','northeast')
subplot(224)
hold on
stairs(-CENTERS, N2_forward,'LineWidth',2);
stairs(-CENTERS, N2_reverse,'LineWidth',2);
xline(0,'LineWidth',t0_lw,'Color','k')
hax=gca;
hax.XLim = [-max_time 0];
hax.YLim(1) = 0;
xlabel('Time from rest end (s)')
ylabel('Events count');
legend("forward","Reverse",'Location','northwest')

sgtitle({'Reverse replay immediately after landing but no forward replay before takeoff';...
         sprintf('bin width = %.2f',time_bin_width)});
filename = fullfile('F:\sequences\figures\',sprintf('timing_of_reverse_vs_forward_replay_in_rest_bin_size_%d_ms.jpg',round(time_bin_width*1e3)));
saveas(fig,filename,'jpg');
end

%%
IX = true(size(seqs));
IX = IX & ~[seqs.forward];
IX = IX & [events.time_from_epoch_start]<25;
figure
tiledlayout('flow')
nexttile
plot([events(IX).time_from_epoch_start], [seqs(IX).compression],'.')
xlabel('Time from rest start (s)'); ylabel('Compression')
nexttile
plot([events(IX).time_from_epoch_start], [seqs(IX).distance],'.')
xlabel('Time from rest start (s)'); ylabel('distance (m)')
nexttile
plot([events(IX).time_from_epoch_start], [seqs(IX).duration],'.')
xlabel('Time from rest start (s)'); ylabel('duration (s)')
nexttile
plot([events(IX).time_from_epoch_start], [seqs(IX).score],'.')
xlabel('Time from rest start (s)'); ylabel('score')

%% classify events
gEarlyLate = categorical([events.time_from_epoch_start]<1,[false,true],["Late","Early"])';
gForRev = categorical([seqs.forward],[true false],{'Forward','Reverse'})';
gTakeLand = classify_replay_landing_takeoff_other(seqs);
gPastFuture = categorical( ([events.rest_ball_num] == 1 & [seqs.state_direction] == 1) | ...
                           ([events.rest_ball_num] == 2 & [seqs.state_direction] == -1), ...
    [true false],{'Future','Past'})';

%%
fig=figure; fig.WindowState = 'maximized';
tiledlayout('flow');
nexttile; histogram(gEarlyLate .* gForRev .* gTakeLand .* gPastFuture, 'DisplayOrder','descend');
nexttile; histogram(gEarlyLate .* gForRev .* gTakeLand,     'DisplayOrder','descend');
nexttile; histogram(gEarlyLate .* gForRev .* gPastFuture,   'DisplayOrder','descend');
nexttile; histogram(gEarlyLate .* gTakeLand .* gPastFuture, 'DisplayOrder','descend');
nexttile; histogram(gForRev .* gTakeLand .* gPastFuture,    'DisplayOrder','descend');
nexttile; histogram(gEarlyLate .* gForRev,     'DisplayOrder','descend');
nexttile; histogram(gEarlyLate .* gTakeLand,   'DisplayOrder','descend');
nexttile; histogram(gEarlyLate .* gPastFuture, 'DisplayOrder','descend');
nexttile; histogram(gForRev .* gTakeLand,      'DisplayOrder','descend');
nexttile; histogram(gForRev .* gPastFuture,    'DisplayOrder','descend');
nexttile; histogram(gTakeLand .* gPastFuture,  'DisplayOrder','descend');
nexttile; histogram(gEarlyLate , 'DisplayOrder','descend');
nexttile; histogram(gForRev,     'DisplayOrder','descend');
nexttile; histogram(gTakeLand ,  'DisplayOrder','descend');
nexttile; histogram(gPastFuture, 'DisplayOrder','descend');
sgtitle('Rest replay classifications')
filename = 'F:\sequences\figures\Rest_replay_classifications';
saveas(fig,filename,'fig');
saveas(fig,filename,'jpg');

%%
X = [[seqs.compression]; 
     [seqs.distance];
     [seqs.duration];
     [seqs.score];
     [seqs.middle_pos_norm];]';
features_names = ["compression","distance","duration","score","middle_pos_norm",];
% g = gEarlyLate .* gForRev .* gTakeLand .* gPastFuture;
g = gEarlyLate .* gTakeLand .* gPastFuture;
g = categorical(g=='Early Landing Past',[false true],["other","Early Landing Past"]);
% g = gEarlyLate .* gPastFuture;
% g = gEarlyLate;
nCat = length(categories(g));
clrs = brewermap(nCat, 'Set1');
clrs=[];
figure
gplotmatrix(X,[],g,clrs,'o+d*',5,[],[],features_names);

%%
figure
plot([events.rest_ball_num],[seqs.middle_pos_norm],'.')


    
%%