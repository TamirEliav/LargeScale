%%
clear
clc
close all

%% params
epoch_type = 'rest';
params_opt = 11;
event_type = 'posterior';
features_fn = {
    'score',
    'compression',
    'duration',
    'distance',
    'forward',
    'g_exp_ID',
    'prev_co_click_rate',
    'prev_co_time_diff',
    'prev_co_same_map',
    };

%% choose exp days for 2 bats with good behavioral decoding
exp_list = {
    'b2299_d191202',
    'b2299_d191203',
    'b2299_d191204',
    'b2299_d191205',
    'b2299_d191206',
    'b2299_d191207', % no data
    'b2299_d191208', % no data
    'b2299_d191209',
    'b2299_d191210',
%     'b2299_d191211', % 
%     'b2299_d191212',
    'b2299_d191213',
    }


%% load data
events_all_sessions = {};
for ii_exp = 1:length(exp_list)
    %%
    try 
    exp_ID = exp_list{ii_exp};
    exp = exp_load_data(exp_ID, 'details', 'pos');
    [events, params]= decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
    [seqs, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
    events(~TF)=[];
    co = exp.pos.proc_1D.co;
    prev_co_IX = interp1(co.ts, 1:length(co.ts), [events.peak_ts], "previous","extrap");
    prev_co_invalid = isnan(prev_co_IX);
    prev_co_IX(prev_co_invalid)=1; % workaround because matlab don't allow using nan as index, we will remove those invalid events at the end
    prev_co_ts = co.ts(prev_co_IX);
    prev_co_direction = co.direction(prev_co_IX);
    prev_co_direction(prev_co_direction==2)=-1;
    [events.prev_co_time_diff] = disperse(([events.peak_ts] - prev_co_ts').*1e-6);
    [events.prev_co_same_map] = disperse(prev_co_direction'==[seqs.state_direction]);
    [events.prev_co_pos] = disperse(co.pos(prev_co_IX));
    [events.prev_co_direction] = disperse(co.direction(prev_co_IX));
    [events.prev_co_click_rate] = disperse(co.click_rate(prev_co_IX));
    [events.prev_co_excluded_from_audio_analysis] = disperse(co.excluded_from_audio_analysis(prev_co_IX));
    next_co_IX = interp1(co.ts, 1:length(co.ts), [events.peak_ts], "next","extrap");
    next_co_invalid = isnan(next_co_IX);
    next_co_IX(next_co_invalid)=1; % workaround because matlab don't allow using nan as index, we will remove those invalid events at the end
    next_co_ts = co.ts(next_co_IX);
    [events.next_co_pos] = disperse(co.pos(next_co_IX));
    [events.exp_ID] = deal(exp_ID);
    events(prev_co_invalid)=[];
    events_all_sessions{ii_exp} = events;
    catch err
        fprintf('something happened (%s):',exp_ID);
        getReport(err)
    end
end
events_all = [events_all_sessions{:}];
seqs_all = [events_all.seq_model];
[seqs_all.g_exp_ID] = disperse(findgroups({events_all.exp_ID}));
[seqs_all.prev_co_click_rate] = disperse([events_all.prev_co_click_rate]);
[seqs_all.prev_co_time_diff] = disperse([events_all.prev_co_time_diff]);
[seqs_all.prev_co_same_map] = disperse([events_all.prev_co_same_map]);
[seqs_all.prev_co_pos] = disperse([events_all.prev_co_pos]);
[seqs_all.next_co_pos] = disperse([events_all.next_co_pos]);

%% prepare data info 
data_info = struct();
data_info.exp_ID = string({events_all.exp_ID});
data_info.ts = [events_all.peak_ts];
data_info.evnet_num = [events_all.num];

%%
replay_pos_limits = [25 115];
fig=figure;
fig.WindowState = 'maximized';
tiledlayout('flow','Padding','tight');
X = [events_all.prev_co_pos];
Y = [seqs_all.middle_pos];
for ii_fn = 1:length(features_fn)
    nexttile
    hold on
    fn = features_fn{ii_fn};
    c = [seqs_all.(fn)];
    scatter(X,Y,20,c,'filled');
    refline(1,0)
    xlabel('Last crossover position (m)')
    ylabel('Replay position (m)')
    hb=colorbar;
    hb.Label.String = fn;
    hb.Label.FontSize = 12;
    hb.Label.Interpreter = 'none';
    colormap jet
%     colormap parula
    ylim(replay_pos_limits);
    set(gca,'Color','k')
end
sgtitle(sprintf('Replay position vs. crossover positions (2 bats) - replay pos lim %dm-%dm',replay_pos_limits))
filename = fullfile('F:\sequences\figures\replay_vs_crossover\',sprintf('replay_vs_crossover_opt_%d_pos_limits_%d_%d',params_opt,replay_pos_limits));
saveas(gcf,filename,'fig');
exportgraphics(fig,[filename '.pdf']);

%% calc correlations 
clc
% close all
replay_pos_limits = [25 115];
% replay_pos_limits = [30 110];

TF = Y>replay_pos_limits(1) & Y<replay_pos_limits(2);
msg_str = sprintf('replay position between %d-%d', replay_pos_limits);
run_corr(X,Y,TF,msg_str,params_opt,data_info,seqs_all);

replay_dist_thr = 5;
TF = Y>replay_pos_limits(1) & Y<replay_pos_limits(2) & [seqs_all.distance]>replay_dist_thr;
msg_str = sprintf('replay position between %d-%d & replay distance>%d', replay_pos_limits, replay_dist_thr);
run_corr(X,Y,TF,msg_str,params_opt,data_info,seqs_all);

replay_dist_thr = 10;
TF = Y>replay_pos_limits(1) & Y<replay_pos_limits(2) & [seqs_all.distance]>replay_dist_thr;
msg_str = sprintf('replay position between %d-%d & replay distance>%d', replay_pos_limits, replay_dist_thr);
run_corr(X,Y,TF,msg_str,params_opt,data_info,seqs_all);

TF = Y>replay_pos_limits(1) & Y<replay_pos_limits(2) & [seqs_all.forward];
msg_str = sprintf('replay position between %d-%d & forward', replay_pos_limits);
run_corr(X,Y,TF,msg_str,params_opt,data_info,seqs_all);

TF = Y>replay_pos_limits(1) & Y<replay_pos_limits(2) & ~[seqs_all.forward];
msg_str = sprintf('replay position between %d-%d & reverse', replay_pos_limits);
run_corr(X,Y,TF,msg_str,params_opt,data_info,seqs_all);

TF = Y>replay_pos_limits(1) & Y<replay_pos_limits(2) & [seqs_all.prev_co_same_map];
msg_str = sprintf('replay position between %d-%d & same map', replay_pos_limits);
run_corr(X,Y,TF,msg_str,params_opt,data_info,seqs_all);

TF = Y>replay_pos_limits(1) & Y<replay_pos_limits(2) & [seqs_all.prev_co_same_map] & [seqs_all.forward];
msg_str = sprintf('replay position between %d-%d & same map & forward', replay_pos_limits);
run_corr(X,Y,TF,msg_str,params_opt,data_info,seqs_all);

TF = Y>replay_pos_limits(1) & Y<replay_pos_limits(2) & [seqs_all.prev_co_same_map] & ~[seqs_all.forward];
msg_str = sprintf('replay position between %d-%d & same map & reverse', replay_pos_limits);
run_corr(X,Y,TF,msg_str,params_opt,data_info,seqs_all);

%%
replay_dist_thr = 10;
TF = Y>replay_pos_limits(1) & Y<replay_pos_limits(2) & [seqs_all.distance]>replay_dist_thr;
msg_str = sprintf('replay position between %d-%d & replay distance>%d', replay_pos_limits, replay_dist_thr);
xxx = [seqs_all.prev_co_time_diff];
yyy = abs(X-Y);
ccc = [seqs_all.prev_co_same_map];
figure
scatter(xxx(TF),yyy(TF),10,ccc(TF),'filled')
xlabel('Time diff')
ylabel('Position diff')
xlim([0 100])
h=colorbar
h.Label.String = 'same map?'
cmap = colormap;
cmap = flip(cmap);
cmap = [1 0 0; 0 0 1];
colormap(cmap)

%% weighted linear regression
% clc
% TF = Y>replay_pos_limits(1) & Y<replay_pos_limits(2) & ~[seqs_all.forward];
TF = Y>replay_pos_limits(1) & Y<replay_pos_limits(2);
W = [seqs_all.score];
% W = [seqs_all.distance];
% W = ones(size(X));
% W = [seqs_all.compression];
x = X(TF);
y = Y(TF);
w = W(TF);
lm = fitlm(x,y,Weights=w)
figure
hold on
lm.plot();

%%
function stats = calc_corr(x,y,TF)
arguments
    x
    y
    TF=[]
end
if isempty(TF)
    TF = true(size(x));
end
x = x(TF);
y = y(TF);
x=x(:);
y=y(:);
tail = 'right';
[stats.Pearson.r, stats.Pearson.p] = corr(x,y,'type','Pearson',rows='pairwise',tail=tail);
[stats.Spearman.r, stats.Spearman.p] = corr(x,y,'type','Spearman',rows='pairwise',tail=tail);
stats.tail = tail;
stats.n = sum(~isnan(x) & ~isnan(y));
end

%%
function report_corr_res(stats,msg_str)
fprintf('-------------------------------------------------------------\n');
fprintf('Filters: %s\n',msg_str);
fprintf('Pearson: r=%.2g, p=%.2g\n',stats.Pearson.r, stats.Pearson.p);
fprintf('Spearman: r=%.2g, p=%.2g\n',stats.Spearman.r, stats.Spearman.p);
fprintf('n=%d\n',stats.n);
fprintf('tail=%s\n',stats.tail);
% fprintf('\n');
end

%%
function print_scatter(x,y,TF,msg_str,stats,params_opt,data_info,seqs_all)
dir_out = 'F:\sequences\figures\replay_vs_crossover';
filename = fullfile(dir_out, sprintf('replay_vs_crossover_opt_%d_%s',params_opt,msg_str) );
filename = strrep(filename, '>','_gt_');
fig=figure;
hold on
axis equal
hs = scatter(x(TF),y(TF),20,'k','filled');
row = dataTipTextRow('exp_ID',data_info.exp_ID(TF)); hs.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('ts',data_info.ts(TF),'%d'); hs.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('event_num',data_info.evnet_num(TF),'%d'); hs.DataTipTemplate.DataTipRows(end+1) = row;
xlim([10 125])
ylim([10 125])
refline(1,0);
xlabel('Last crossover position (m)','FontSize',16);
ylabel('Replay position (m)','FontSize',16);
title(msg_str,'FontSize',16);
text(0.88,0.30, sprintf('n=%d\n',stats.n), 'Units','normalized','FontSize',9);
text(0.88,0.20, sprintf('Pearson: \n r=%.2g, p=%.2g\n',stats.Pearson.r, stats.Pearson.p), 'Units','normalized','FontSize',9);
text(0.88,0.05, sprintf('Spearman: \n r=%.2g, p=%.2g\n',stats.Spearman.r, stats.Spearman.p), 'Units','normalized','FontSize',9);
mkdir(dir_out);
exportgraphics(fig,[filename '.pdf']);
save(filename,'x','y','TF','stats','msg_str','params_opt','seqs_all','data_info');
end


%%
function run_corr(x,y,TF,msg_str,params_opt,data_info,seqs_all)
stats = calc_corr(x,y,TF);
report_corr_res(stats,msg_str);
print_scatter(x,y,TF,msg_str,stats,params_opt,data_info,seqs_all);

end













%%
