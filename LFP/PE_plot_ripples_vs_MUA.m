function PE_plot_ripples_vs_MUA(exp_ID)

%% get exp info
exp = exp_load_data(exp_ID,'details','path','ripples','MUA');
prm = PARAMS_GetAll();
win_sec = 0.5;

%% trigger the MUA by ripple times
rpl = exp.ripples.events;
rpl_ts = [rpl.peak_ts];
MUA_ts = [exp.MUA.events.peak_ts];

% trig MUA FR by ripple evnets
win_n_sample = round(win_sec * exp.MUA.fs);
rpl_t_IX = find_nearest_point(rpl_ts, exp.MUA.t);
trig_MUA.FR = trigger_signal_by_IX(exp.MUA.FR, rpl_t_IX, win_n_sample, 0);
trig_MUA.zFR = trigger_signal_by_IX(exp.MUA.zFR, rpl_t_IX, win_n_sample, 0);
trig_MUA.t = linspace(-win_sec,win_sec,size(trig_MUA.zFR,2));

% trig ripple power (z) by MUA evnets
win_n_sample = round(win_sec * exp.ripples.fs);
MUA_t_IX = find_nearest_point(MUA_ts, exp.ripples.t);
trig_ripple.zpripple = trigger_signal_by_IX(exp.ripples.zpripple, MUA_t_IX, win_n_sample, 0);
trig_ripple.t = linspace(-win_sec,win_sec,size(trig_ripple.zpripple,2));

% why we do the following 3 lines??...
trig_MUA.FR(isnan(trig_MUA.FR)) = 0;
trig_MUA.zFR(isnan(trig_MUA.zFR)) = 0;
trig_ripple.zpripple(isnan(trig_ripple.zpripple)) = 0;

%% plot
% h={};
hf=figure;
hf.Units = 'centimeters';
hf.Position = [5 5 40 20];
% tiledlayout('flow','TileSpacing', 'compact')
pnl = panel();
pnl.pack('h',[.2 .3 .3 .2]);
pnl(1).pack('v',3);
pnl(2).pack('v',3);
pnl(3).pack('v',3);
pnl(4).pack('v',3);
pnl.margin = [20 20 30 20];
% pnl.de.margin = [30 15 30 15];
pnl(1).margin = [20 15 20 15];
pnl(2).margin = [20 15 30 15];
pnl(3).margin = [35 15 35 15];
pnl(4).margin = [20 15 20 15];
clrs = [1 0 0; 0 0 1];

% ------ MUA zFR trig by ripples events -----------

% cluster
X=[];
X(:,end+1) = [rpl.ripple_gamma_power_ratio_at_peak];
X(:,end+1) = [rpl.ripple_gamma_power_ratio_mean];
X(:,end+1) = exp.MUA.zFR(rpl_t_IX);
X = min(X,prctile(X,[99])); % trim outliers
dist_metric = 'sqeuclidean';
rng(0)
g = kmeans(X,2,'Distance',dist_metric,'Replicates',100);
g(isnan(g)) = 0;
g(g==0) = 3;
[~,g_sort_IX] = sort(splitapply(@(x)(nanmean(x,'all')),X,g),'descend'); % make sure cluster 1 contains high values
g2 = zeros(size(g));
for ii = 1:length(g_sort_IX)
    g2(g==ii) = g_sort_IX(ii);
end
g=g2;
clear g2;
g(g==3) = 2; % workaround for missing data points (nan)
[~,sort_IX]=sort(g);
grp_rpl = g;

% trig mean
pnl(2,1).select(); hold on;
clear h
h(1) = shadedErrorBar(trig_MUA.t, trig_MUA.zFR, {@nanmean,@nansem});
h(2) = shadedErrorBar(trig_MUA.t, trig_MUA.zFR(g==1,:), {@nanmean,@nansem},'lineprops','r');
h(3) = shadedErrorBar(trig_MUA.t, trig_MUA.zFR(g==2,:), {@nanmean,@nansem},'lineprops','b');
legend([h.mainLine], 'All','Strong cluster','Weak cluster');
xlabel('Time relative to ripple peak (s)')
ylabel('MUA firing rate (z)')

% trig image 
pnl(2,2).select(); hold on
imagesc(trig_MUA.t, 1:size(trig_MUA.zFR,1), trig_MUA.zFR);
y_ticks = [1 size(trig_MUA.zFR,1)];
xlim(trig_MUA.t([1 end]));
ylim([0 size(trig_MUA.zFR,1)+1]);
yticks(y_ticks);
xlabel('Time relative to ripple peak (s)')
ylabel({'Ripple events';'ordered by time'},'Units','Normalized','Position',[-0.03 0.5]);
h = colorbar;
h.Location = 'manual';
h.Units = 'normalized';
h.TickLength = 0.03;
h.Position = [0.44 0.45 0.005 0.1];
h.Label.String = 'MUA firing rate (z)';
h.Label.Units = 'normalized';
h.Label.Position(1) = -h.Label.Position(1);
hax = gca;
hax.CLim(2) = 5;

% trig image (reordered by clusters)
pnl(2,3).select(); hold on
imagesc(trig_MUA.t, 1:size(trig_MUA.zFR,1), trig_MUA.zFR(sort_IX,:));
y_ticks = [1 sum(g==1) size(trig_MUA.zFR,1)];
plot(-win_sec*[1 1], y_ticks([1 2]),'r','Clipping','off','LineWidth',2);
plot(-win_sec*[1 1], y_ticks([2 3]),'b','Clipping','off','LineWidth',2);
xlim(trig_MUA.t([1 end]));
ylim([0 size(trig_MUA.zFR,1)+1]);
xlabel('Time relative to ripple peak (s)')
ylabel({'Ripple events';'ordered by clusters'},'Units','Normalized','Position',[-0.1 0.5]);
yticks(y_ticks)
h = colorbar;
h.Location = 'manual';
h.Units = 'normalized';
h.TickLength = 0.03;
h.Position = [0.44 0.16 0.005 0.1];
h.Label.String = 'MUA firing rate (z)';
h.Label.Units = 'normalized';
h.Label.Position(1) = -h.Label.Position(1);
hax = gca;
hax.CLim(2) = 5;
hax.TickDir = 'out';

% ---------------------------------------------
pnl(1,1).select()
x = [rpl.peak_zpripple];
y = exp.MUA.zFR(rpl_t_IX);
% IX = kmeans([x;y]',2,'Replicates',1000,'Distance','cosine');
scatter(x,y,20,clrs(g,:),'.')
xlabel('Ripple peak amplitude (z)')
ylabel('MUA at ripple peak (z)')

pnl(1,2).select()
x = [rpl.peak_zpripple];
y = [rpl.ripple_gamma_power_ratio_at_peak];
scatter(x,y,20,clrs(g,:),'.')
xlabel('Ripple peak amplitude (z)')
ylabel('Ripple-gamma ratio (@peak)')
hax=gca;
hax.YScale = 'log';

pnl(1,3).select()
x = [rpl.peak_zpripple];
y = [rpl.ripple_gamma_power_ratio_mean];
scatter(x,y,20,clrs(g,:),'.')
xlabel('Ripple peak amplitude (z)')
ylabel('Ripple-gamma ratio (mean)')
hax=gca;
hax.YScale = 'log';



% ------ ripples trig by MUA events -----------

% cluster
X=[];
X(:,end+1) = [exp.MUA.events.peak_zFR];
X(:,end+1) = exp.ripples.zpripple(MUA_t_IX);
X = min(X,prctile(X,[99])); % trim outliers
dist_metric = 'sqeuclidean';
% dist_metric = 'correlation';
% dist_metric = 'cosine';
rng(0)
g = kmeans(X,2,'Distance',dist_metric,'Replicates',1000);
g(isnan(g)) = 0;
g(g==0) = 3;
[~,g_sort_IX] = sort(splitapply(@(x)(nanmean(x,'all')),X,g),'descend'); % make sure cluster 1 contains high values
g2 = zeros(size(g));
for ii = 1:length(g_sort_IX)
    g2(g==ii) = g_sort_IX(ii);
end
g=g2;
clear g2;
g(g==3) = 2; % workaround for missing data points (nan)
[~,sort_IX]=sort(g);
grp_mua = g;

pnl(3,1).select(); hold on
clear h
h(1) = shadedErrorBar(trig_ripple.t, trig_ripple.zpripple, {@nanmean,@nansem});
h(2) = shadedErrorBar(trig_ripple.t, trig_ripple.zpripple(g==1,:), {@nanmean,@nansem},'lineprops','r');
h(3) = shadedErrorBar(trig_ripple.t, trig_ripple.zpripple(g==2,:), {@nanmean,@nansem},'lineprops','b');
legend([h.mainLine], 'All','Strong cluster','Weak cluster');
xlabel('Time relative to MUA peak (s)')
ylabel('Ripples power (z)')

pnl(3,2).select()
imagesc(trig_ripple.t, 1:size(trig_ripple.zpripple,1), trig_ripple.zpripple);
xlim(trig_ripple.t([1 end]));
ylim([0 size(trig_ripple.zpripple,1)+1]);
xlabel('Time relative to MUA peak (s)')
ylabel({'MUA events';'ordered by time'},'Units','Normalized','Position',[-0.03 0.5]);
yticks([1 size(trig_ripple.zpripple,1)])
hax = gca;
hax.CLim(2) = 5;
h = colorbar;
h.Location = 'manual';
h.Units = 'normalized';
h.TickLength = 0.03;
h.Position = [0.73 0.45 0.005 0.1];
h.Label.String = 'Ripples power (z)';
h.Label.Units = 'normalized';
h.Label.Position(1) = -h.Label.Position(1);

% trig image (ordered)
pnl(3,3).select(); hold on
imagesc(trig_ripple.t, 1:size(trig_ripple.zpripple,1), trig_ripple.zpripple(sort_IX,:));
y_ticks = [1 sum(g==1) size(trig_ripple.zpripple,1)];
plot(-win_sec*[1 1], y_ticks([1 2]),'r','Clipping','off','LineWidth',2);
plot(-win_sec*[1 1], y_ticks([2 3]),'b','Clipping','off','LineWidth',2);
xlim(trig_ripple.t([1 end]));
ylim([0 size(trig_ripple.zpripple,1)+1]);
xlabel('Time relative to MUA peak (s)')
% ylabel({'Ripple events';'ordered by time'},'Units','Normalized','Position',[-0.1 0.5]);
ylabel({'MUA events';'ordered by clusters'});
yticks(y_ticks)
h = colorbar;
h.Location = 'manual';
h.Units = 'normalized';
h.TickLength = 0.03;
h.Position = [0.73 0.16 0.005 0.1];
h.Label.String = 'Ripples power (z)';
h.Label.Units = 'normalized';
h.Label.Position(1) = -h.Label.Position(1);
hax = gca;
hax.CLim(2) = 5;
hax.TickDir = 'out';

% ---------------------------------------------
pnl(4,1).select()
x = [exp.MUA.events.peak_zFR];
y = exp.ripples.zpripple(MUA_t_IX);
scatter(x,y,20,clrs(g,:),'.')
ylabel('Ripple power at MUA peak (z)')
xlabel('MUA peak FR (z)')

% ----------------- MUA-ripples events xcorr ---------
pnl(4,2).select(); hold on
tdiff = MUA_ts - rpl_ts';
bins = linspace(-1000,1000,40);
h=histogram(tdiff(:).*1e-3, bins);
h.FaceColor = 'k';
xlabel('Time relative to ripple peak (ms)')
ylabel('Counts')
title('Ripple vs MUA events timing')
text(0.55,0.9,"#ripples = "+size(tdiff,1),'Units','normalized')
text(0.55,0.8,"#MUA = "+size(tdiff,2),'Units','normalized')
text(0.55,0.7,"BinWidth = "+round(h.BinWidth)+"ms",'Units','normalized')

pnl(4,3).select(); hold on
tdiff1 = MUA_ts(grp_mua==1) - rpl_ts((grp_rpl==1))';
tdiff2 = MUA_ts(grp_mua==2) - rpl_ts((grp_rpl==2))';
% tdiff2 = MUA_ts(grp_mua==1) - rpl_ts';
% tdiff2 = MUA_ts - rpl_ts((grp_rpl==1))';
bins = linspace(-1000,1000,40);
h2=histogram(tdiff2(:).*1e-3, bins);
h2.FaceColor = 'b';
h1=histogram(tdiff1(:).*1e-3, bins);
h1.FaceColor = 'r';
xlabel('Time relative to ripple peak (ms)')
ylabel('Counts')
title('Ripple vs MUA events timing')
text(0.05,0.9,"#ripples = "+size(tdiff1,1),'Units','normalized','FontSize',8,'Color','r');
text(0.05,0.8,"#MUA = "+size(tdiff1,2),'Units','normalized','FontSize',8,'Color','r');
text(0.55,0.9,"#ripples = "+size(tdiff2,1),'Units','normalized','FontSize',8,'Color','b');
text(0.55,0.8,"#MUA = "+size(tdiff2,2),'Units','normalized','FontSize',8,'Color','b');
text(0.55,0.7,"BinWidth = "+round(h1.BinWidth)+"ms",'Units','normalized','FontSize',8);

% ----------------- figure title ---------
sgtitle({exp.details.exp_ID,'Ripple vs. MUA events'},'interpreter','None')

% ----------------- save figure -------------
fig_filename = fullfile('L:\Analysis\Results\exp\PE', [exp_ID '_PE_ripples_vs_MUA']);
% saveas(gcf,fig_filename,'fig')
saveas(gcf,fig_filename,'jpeg')

%% merge ripples and MUA to PE
rpl = exp.ripples.events;
mua = exp.MUA.events;
[rpl.peak_zFR] = disperse(interp1(exp.MUA.t,exp.MUA.zFR, [rpl.peak_ts]));
[mua.peak_zpripple] = disperse(interp1(exp.ripples.t,exp.ripples.zpripple, [mua.peak_ts]));
t = exp.MUA.t; % use MUA timestamps
PE_all = merge_ripples_and_MUA(rpl,mua,t);
PE_strong = merge_ripples_and_MUA(rpl(grp_rpl==1),mua(grp_mua==1),t);
PE_thr = merge_ripples_and_MUA(...
    rpl([rpl.peak_zFR]      > prm.MUA.high_thr_std), ...
    mua([mua.peak_zpripple] > prm.ripples.high_thr_std), t);

%% save ripples detection results
PE = struct();
PE.all = PE_all;
PE.strong = PE_strong;
PE.thr = PE_thr;
file_name = fullfile('L:\Analysis\Results\exp\PE',[exp_ID '_exp_PE']);
save(file_name,'PE');


%%

end



%%
function PE = merge_ripples_and_MUA(rpl,mua,t)

%%
mua_ti = [mua.start_ts; mua.end_ts]';
rpl_ti = [rpl.start_ts; rpl.end_ts]';
is_mua = any(t>=mua_ti(:,1)&t<=mua_ti(:,2),1);
is_rpl = any(t>=rpl_ti(:,1)&t<=rpl_ti(:,2),1);

%%
mua_lbl = bwlabel(is_mua);
% max(mua_lbl)
if max(mua_lbl) ~= length(mua)
    error('Something went wrong...');
end
valid_events = unique(mua_lbl .* is_rpl);
valid_events(valid_events==0)=[];
% length(valid_events) / length(mua.events)
PE = mua(valid_events);
PE = rmfield(PE,'start_IX');
PE = rmfield(PE,'end_IX');
PE = rmfield(PE,'peak_IX');
% TODO: add fields for ripples quantifications (z/ripple-gamma-ratio/...)


end


%%





