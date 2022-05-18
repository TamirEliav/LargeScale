function decoding_plot_MAP(decode,flight_decoding_param_opt)
arguments
    decode
    flight_decoding_param_opt = 4
end
%%
exp_ID = decode.exp_ID;
epoch_type = decode.epoch_type;
params_opt = decode.params_opt;

%% load data
exp = exp_load_data(exp_ID, 'ripples','MUA','PE','LM','pos');
if strcmp(epoch_type, 'flight')
    decode_flight = decode;
else
    decode_flight = decoding_load_data(exp_ID, 'flight', flight_decoding_param_opt);
end
dir_OUT = 'F:\sequences\decoded_figs\MAP_entire_session';
mkdir(dir_OUT);

%%
bin_centers = decode.pos;
bin_edges = centers2edges(bin_centers);

%% plot figure
fig=figure;
fig.WindowState = 'maximized';
pnl = panel();
pnl.pack('v',[.15 0.15 .7])
pnl(1).pack('h',[.85 .05 .05]);
pnl(2).pack('h',[.85 .05 .05]);
pnl(3).pack('h',[.85 .05 .05]);
pnl.de.margin = 10;
pnl.margintop = 15;
% pnl.margin = 25;
h = pnl.title({exp_ID;sprintf('MAP estimate - %s',epoch_type)});
h.Interpreter = 'none';
h.FontSize = 12;

% MAP positio/state
pnl(3,1).select();
hax=gca;
hold on
if strcmp(epoch_type,'flight')
    pos_ts_IX = interp1(decode.time, 1:length(decode.time), exp.pos.proc_1D.ts, 'previous');
    plot(pos_ts_IX,exp.pos.proc_1D.pos,'r-','LineWidth',0.01);
end
clear h
hax.ColorOrderIndex = 1;
h = arrayfun(@(x)(plot(x,nan,'.','MarkerSize',20)),[1:length(decode.state)]); % dummy points for legend
x = 1:length(decode.time);
y = decode.MAP_pos;
g = decode.MAP_state_IX;
nStates = length(decode.state);
g = [g 1:nStates];
x = [x nan(1,nStates)];
y = [y nan(1,nStates)];
hax.ColorOrderIndex = 1;
splitapply(@(IX,x)(plot(IX,x,'.')),x, y, g);
gaps_IX = find(diff(decode.time) > median(diff(decode.time)));
arrayfun(@(x)(xline(x,'g','time gap','LineWidth',0.5)), gaps_IX);
rescale_plot_data('x',[1/decode.Fs 0]);
arrayfun(@(y,str)(yline(y,'-',str,'Color',0.5.*[1 1 1],'LineWidth',0.5,'LabelVerticalAlignment','middle')),[exp.LM.pos_proj], string({exp.LM.name}))
hl = legend(h,decode.state,'Interpreter','none','Location','northoutside');
hl.Position = [0.84 0.87 0.1 0.1];
xlabel('Time (s)')
ylabel('Position (m)')
ylimits = [min(decode.pos) max(decode.pos)];
if range(ylimits)<10
    ylimits = ylimits + 0.2.*range(ylimits).*[-1 1];
end
ylim(ylimits);

% MAP pos hist
pnl(3,2).select();
hold on
for ii_state = 1:length(decode.state)
    if contains(decode.state(ii_state),'empirical')
        IX = g == ii_state;
        histogram(y(IX), 'BinEdges',bin_edges, 'Orientation','horizontal','Normalization','probability','DisplayStyle','stairs','LineWidth',1.5);
    else
        plot(nan,nan); % to keep the color progression...
    end
end
arrayfun(@(y,str)(yline(y,'Color',0.5.*[1 1 1],'LineWidth',0.5)),[exp.LM.pos_proj], string({exp.LM.name}))
title({decode.epoch_type;'over representation'})
ylim(ylimits);

% flight error hist
pnl(3,3).select();
% TODO: take from params function
err_thr_prc = 5;
err_thr_m = range(bin_edges)*err_thr_prc/100;
pos_real = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos, decode_flight.time);
pos_predict = decode_flight.MAP_pos;
TF = abs(pos_real-pos_predict) < err_thr_m;
% histogram(pos_predict(~TF), 'BinEdges',bin_edges, 'Orientation','horizontal','Normalization','probability');
N_error = histcounts(pos_predict(~TF), 'BinEdges',bin_edges, 'Normalization','count');
% N_total = histcounts(pos_predict,      'BinEdges',bin_edges, 'Normalization','count');
% predict_err_prob = N_error./N_total;
predict_err_prob = N_error./length(pos_predict);
barh(bin_centers,predict_err_prob,'r');
arrayfun(@(y,str)(yline(y,'Color',0.5.*[1 1 1],'LineWidth',0.5)),[exp.LM.pos_proj], string({exp.LM.name}));
title({decode_flight.epoch_type;'Predict errors'})
ylim(ylimits);

linkaxes(pnl(3).de.axis,'y')

% -------- ripple power and MUA firing rate
% arrange data
MUA_ts_IX = interp1(decode.time, 1:length(decode.time), exp.MUA.t, 'previous');
ripples_ts_IX = interp1(decode.time, 1:length(decode.time), exp.ripples.t, 'previous');
PE_events_ts_IX = interp1(decode.time, 1:length(decode.time), [exp.PE.thr.peak_ts], 'previous');
MUA_ts_IX(ismember(MUA_ts_IX, gaps_IX)) = nan; % remove points out of the current decoding epoch time
PE_events_ts_IX(ismember(PE_events_ts_IX, gaps_IX)) = nan; % remove points out of the current decoding epoch time
% ripple
pnl(1,1).select();
hold on
plot(ripples_ts_IX , exp.ripples.zpripple);
plot(PE_events_ts_IX, [exp.PE.thr.peak_zpripple],'*r');
arrayfun(@(x)(xline(x,'g','LineWidth',0.5)), gaps_IX);
legend("Ripple","Pop events")
xlim([0 1.05*length(decode.time)])
ylim([0 10])
rescale_plot_data('x',[1/decode.Fs 0]);
xlabel('Time (s)')
ylabel('Ripple power (z)')
% MUA 
pnl(2,1).select();
hold on
plot(MUA_ts_IX, exp.MUA.zFR);
plot(PE_events_ts_IX, [exp.PE.thr.peak_zFR],'*r');
arrayfun(@(x)(xline(x,'g','LineWidth',0.5)), gaps_IX);
legend("MUA","Pop events")
xlim([0 1.05*length(decode.time)])
ylim([0 10])
rescale_plot_data('x',[1/decode.Fs 0]);
xlabel('Time (s)')
ylabel('Firing rate (z)')
% link time
linkaxes([pnl(1,1).axis pnl(2,1).axis pnl(3,1).axis],'x')

params_str = {
    sprintf('params opt: %d',params_opt),
    sprintf('bin size: %.2gm',decode.params.pos_bin_size),
    sprintf('replay speed: x%d',decode.params.replay_speed),
    sprintf('state decay timescale: %.3g s',decode.params.state_decay_timescale),
    };
annotation('textbox', [0.84 0.8 0.2 0.05], 'String',params_str,'LineStyle','None')

%%
fig_name = sprintf('%s_MAP_%s_opt_%d',exp_ID, epoch_type, params_opt);
filename = fullfile(dir_OUT, fig_name);
saveas(fig, filename , 'jpg');


%%

