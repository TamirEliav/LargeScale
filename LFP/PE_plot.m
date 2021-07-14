function PE_plot(exp_ID)

%% get exp info
exp = exp_load_data(exp_ID,'details','ripples','MUA','PE');

%%
PE = exp.PE;
[~,sort_IX] = sort([PE.peak_zFR],'descend');
PE = PE(sort_IX);

%% trigger MUA by PE peak/start/end
win_s = 0.5;
fs = exp.MUA.fs;
win_n_samples = round(win_s * fs);
trig_IX_by_peak = find_nearest_point([PE.peak_ts], exp.MUA.t);
trig_IX_by_start = find_nearest_point([PE.start_ts], exp.MUA.t);
trig_IX_by_end = find_nearest_point([PE.end_ts], exp.MUA.t);
trig_MUA.traces.by_peak = trigger_signal_by_IX(exp.MUA.zFR, trig_IX_by_peak, win_n_samples);
trig_MUA.traces.by_start = trigger_signal_by_IX(exp.MUA.zFR, trig_IX_by_start, win_n_samples);
trig_MUA.traces.by_end = trigger_signal_by_IX(exp.MUA.zFR, trig_IX_by_end, win_n_samples);
trig_MUA.t = linspace(-win_s,win_s,size(trig_MUA.traces.by_peak,2));

%% trigger zpripple by PE peak/start/end
win_s = 0.5;
fs = 1e6/median(diff(exp.ripples.t));
win_n_samples = round(win_s * fs);
trig_IX_by_peak = find_nearest_point([PE.peak_ts], exp.ripples.t);
trig_IX_by_start = find_nearest_point([PE.start_ts], exp.ripples.t);
trig_IX_by_end = find_nearest_point([PE.end_ts], exp.ripples.t);
trig_ripples.traces.by_peak = trigger_signal_by_IX(exp.ripples.zpripple_all, trig_IX_by_peak, win_n_samples);
trig_ripples.traces.by_start = trigger_signal_by_IX(exp.ripples.zpripple_all, trig_IX_by_start, win_n_samples);
trig_ripples.traces.by_end = trigger_signal_by_IX(exp.ripples.zpripple_all, trig_IX_by_end, win_n_samples);
trig_ripples.t = linspace(-win_s,win_s,size(trig_ripples.traces.by_peak,2));

%% plot MUA triggered FR
hf=figure;
hf.Units = 'centimeters';
hf.Position = [2 2 25 24];
pnl = panel();
pnl.pack(4,2);
pnl.margintop = 20;
pnl.de.margintop = 10;
% plot MUA firing rate
trig = trig_MUA; 
t = trig_MUA.t;
pnl(1,1).select();
traces = trig.traces.by_peak;
imagesc(t, 1:size(traces), traces);
xlabel('Time from peak (ms)')
ylabel('Event')
axis tight
h = colorbar;
h.Location = 'northoutside';
h.Location = 'manual';
h.Units = 'normalized';
h.TickLength = 0.03;
h.TickDirection = 'out';
h.Position = [0.07 0.93 0.1 0.015];
h.Label.String = 'MUA firing rate (z)';
h.Label.Rotation = 0;
pnl(2,1).select();
traces = trig.traces.by_start;
imagesc(t, 1:size(traces), traces);
xlabel('Time from start (ms)')
ylabel('Event')
axis tight
pnl(3,1).select();
traces = trig.traces.by_end;
imagesc(t, 1:size(traces), traces);
xlabel('Time from end (ms)')
ylabel('Event')
axis tight
pnl(4,1).select();
hold on
clear h
h(1)=shadedErrorBar(t,trig.traces.by_peak,{@nanmean,@nansem},'lineprops','k');
h(2)=shadedErrorBar(t,trig.traces.by_start,{@nanmean,@nansem},'lineprops','b');
h(3)=shadedErrorBar(t,trig.traces.by_end,{@nanmean,@nansem},'lineprops','r');
xlabel('Time (ms)')
ylabel('Firing rate (z)')
legend([h.mainLine],{'peak','start','end'})

% plot ripple power
trig = trig_ripples;
t = trig_ripples.t;
pnl(1,2).select();
traces = trig.traces.by_peak;
imagesc(t, 1:size(traces), traces);
hax=gca;
hax.CLim(2) = 7;
xlabel('Time from peak (ms)')
ylabel('Event')
axis tight
h = colorbar;
h.Location = 'northoutside';
h.Location = 'manual';
h.Units = 'normalized';
h.TickLength = 0.03;
h.TickDirection = 'out';
h.Position = [0.85 0.93 0.1 0.015];
h.Label.String = 'Ripple power (z)';
h.Label.Rotation = 0;
pnl(2,2).select();
traces = trig.traces.by_start;
imagesc(t, 1:size(traces), traces);
hax=gca;
hax.CLim(2) = 7;
xlabel('Time from start (ms)')
ylabel('Event')
axis tight
pnl(3,2).select();
traces = trig.traces.by_end;
imagesc(t, 1:size(traces), traces);
hax=gca;
hax.CLim(2) = 7;
xlabel('Time from end (ms)')
ylabel('Event')
axis tight
pnl(4,2).select();
hold on
clear h
h(1)=shadedErrorBar(t,trig.traces.by_peak,{@nanmean,@nansem},'lineprops','k');
h(2)=shadedErrorBar(t,trig.traces.by_start,{@nanmean,@nansem},'lineprops','b');
h(3)=shadedErrorBar(t,trig.traces.by_end,{@nanmean,@nansem},'lineprops','r');
xlabel('Time (ms)')
ylabel('Ripple power (z)')
legend([h.mainLine],{'peak','start','end'})

sgtitle({exp.details.exp_ID;'Population event triggered MUA Firing Rate and ripple power'},'Interpreter','None');

fig_filename = fullfile('L:\Analysis\Results\exp\PE', [exp_ID '_PE_triggered_MUA_FR_zpripple']);
saveas(gcf,fig_filename,'jpeg')



end