function ripples_trigger_MUA(exp_ID)

%% get exp info
exp = exp_load_data(exp_ID,'details','path','ripples','MUA');
prm = PARAMS_GetAll();

%% trigger the MUA by ripple times
rpl = exp.ripples.all;
rpl_ts = [rpl.peak_ts];

win_sec = 0.5;
win_n_sample = round(win_sec * exp.MUA.fs);
rpl_t_IX = find_nearest_point(rpl_ts, exp.MUA.t);
trig_FR = trigger_signal_by_IX(exp.MUA.FR, rpl_t_IX, win_n_sample);
trig_zFR = trigger_signal_by_IX(exp.MUA.zFR, rpl_t_IX, win_n_sample);
trig_t = linspace(-win_sec,win_sec,size(trig_zFR,2));

tdiff = [exp.MUA.events.peak_ts] - rpl_ts';

%%
traces = trig_zFR;
% L=size(traces,2);
% traces = traces(:,0.4*L:0.6*L);
X = traces;
X=[];
X(:,end+1) = [rpl.ripple_gamma_power_ratio_at_peak];
X(:,end+1) = [rpl.ripple_gamma_power_ratio_mean];
X(:,end+1) = exp.MUA.zFR(rpl_t_IX);
% dist_metric = 'sqeuclidean';
dist_metric = 'correlation';
% dist_metric = 'cosine';
g=kmeans(X,2,'Distance',dist_metric,'Replicates',10);
g(isnan(g))=3;
[~,IX]=sort(g);

figure
tiledlayout('flow')
nexttile
imagesc(traces(IX,:))
nexttile
histogram(g)
nexttile
x = [rpl.peak_zpripple];
y = exp.MUA.zFR(rpl_t_IX);
clrs = [1 0 0; 0 0 1;0.5 0.5 0.5];
scatter(x,y,20,clrs(g,:),'.')
sgtitle(dist_metric)
nexttile


%% plot
% h={};
hf=figure;
hf.Units = 'centimeters';
hf.Position = [5 5 25 20];
% tiledlayout('flow','TileSpacing', 'compact')
pnl = panel();
pnl.pack('h',2);
pnl(1).pack('v',3);
pnl(2).pack('v',3);
pnl.margin = [20 20 30 20];
pnl.de.margin = 15;
clrs = [1 0 0; 0 0 1];

% nexttile
pnl(1,1).select()
x = [rpl.peak_zpripple];
y = exp.MUA.zFR(rpl_t_IX);
IX = kmeans([x;y]',2,'Replicates',1000,'Distance','cosine');
scatter(x,y,20,clrs(IX,:),'.')
% h{end+1}=plot(x, y,'.');
xlabel('Ripple peak amplitude (z)')
ylabel('MUA at ripple peak (z)')
title('clustering based on MUA vs. ripple (this scatter)')

% nexttile
pnl(1,2).select()
x = [rpl.peak_zpripple];
y = [rpl.ripple_gamma_power_ratio_at_peak];
scatter(x,y,20,clrs(IX,:),'.')
% h{end+1}=plot(x, y,'.');
xlabel('Ripple peak amplitude (z)')
ylabel('Ripple-gamma ratio (@peak)')
hax=gca;
hax.YScale = 'log';

% nexttile
pnl(1,3).select()
x = [rpl.peak_zpripple];
y = [rpl.ripple_gamma_power_ratio_mean];
scatter(x,y,20,clrs(IX,:),'.')
% h{end+1}=plot(x, y,'.');
xlabel('Ripple peak amplitude (z)')
ylabel('Ripple-gamma ratio (mean)')
hax=gca;
hax.YScale = 'log';

% linkprop([h{:}],'BrushData');

% nexttile
pnl(2,1).select()
h=histogram(tdiff(:).*1e-3,linspace(-1000,1000,40));
xlabel('Time relative to ripple peak (s)')
ylabel('MUA events (counts)')
title('Ripple vs MUA events timing')
text(0.55,0.9,"#ripples = "+length(rpl_ts),'Units','normalized')
text(0.55,0.8,"#MUA = "+length(exp.MUA.events),'Units','normalized')
text(0.55,0.7,"BinWidth = "+round(h.BinWidth)+"ms",'Units','normalized')

pnl(2,2).select()
shadedErrorBar(trig_t,trig_zFR,{@nanmean,@nansem});
xlabel('Time relative to ripple peak (s)')
ylabel('MUA firing rate (z)')

pnl(2,3).select()
imagesc(trig_t, 1:size(trig_zFR,1), trig_zFR);
xlim(trig_t([1 end]));
ylim([0 size(trig_zFR,1)+1]);
h = colorbar;
% h.Location = 'northoutside';
% h.Location = 'southoutside';
h.Location = 'manual';
h.Units = 'normalized';
h.TickLength = 0.03;
h.Position = [0.9 0.15 0.015 0.1];
h.Label.String = 'MUA firing rate (z)';
xlabel('Time relative to ripple peak (s)')
ylabel('MUA Event');


sgtitle({exp.details.exp_ID,'Ripple vs. MUA events'},'interpreter','None')

%% save figure
fig_filename = fullfile('L:\Analysis\Results\exp\ripples', [exp_ID '_ripples_triggered_MUA']);
saveas(gcf,fig_filename,'fig')
saveas(gcf,fig_filename,'jpeg')

end



