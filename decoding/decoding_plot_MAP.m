function decoding_plot_MAP(exp_ID, epoch_type, params_opt)
arguments
    %% 
    exp_ID = 'b9861_d180526'
    epoch_type {mustBeMember(epoch_type,{'sleep','rest','flight'})} = 'sleep'
    params_opt = 11;
end

%% load data
exp = exp_load_data(exp_ID, 'ripples','MUA','PE','LM');
decode = decoding_load_data(exp_ID, epoch_type, params_opt);
if strcmp(epoch_type, 'flight')
    decode_flight = decode;
else
    flight_decoding_param_opt = 4;
    decode_flight = decoding_load_data(exp_ID, 'flight', flight_decoding_param_opt);
end
dir_OUT = 'F:\sequences\decoded_figs\MAP_entire_session';
mkdir(dir_OUT);

%% plot figure
fig=figure;
fig.WindowState = 'maximized';
pnl = panel();
pnl.pack('v',[.3 .7])
pnl(1).pack('h',[.85 .05 .05]);
pnl(2).pack('h',[.85 .05 .05]);
pnl.de.margin = 10;
pnl.margintop = 15;
% pnl.margin = 25;
h = pnl.title({exp_ID;sprintf('MAP estimate - %s',epoch_type)});
h.Interpreter = 'none';
h.FontSize = 12;
pnl(2,1).select();
hold on
clear h
h = arrayfun(@(x)(plot(x,nan,'.','MarkerSize',20)),[1:length(decode.state)]); % dummy points for legend
hax=gca;
hax.ColorOrderIndex = 1;
splitapply(@(IX,x)(plot(IX,x,'.')),1:length(decode.time), decode.MAP_pos, decode.MAP_state_IX)
gaps_IX = find(diff(decode.time) > median(diff(decode.time)));
arrayfun(@(x)(xline(x,'g','time gap','LineWidth',0.5)), gaps_IX);
rescale_plot_data('x',[1/decode.Fs 0]);
arrayfun(@(y,str)(yline(y,'-',str,'Color',0.5.*[1 1 1],'LineWidth',0.5,'LabelVerticalAlignment','middle')),[exp.LM.pos_proj], string({exp.LM.name}))
hl = legend(h,decode.state,'Interpreter','none','Location','northoutside');
hl.Position = [0.84 0.87 0.1 0.1];
xlabel('Time (s)')
ylabel('Position (m)')
pnl(2,2).select(); hold on
histogram(decode.MAP_pos,       'BinWidth',decode.params.pos_bin_size, 'Orientation','horizontal','Normalization','probability');
arrayfun(@(y,str)(yline(y,'Color',0.5.*[1 1 1],'LineWidth',0.5)),[exp.LM.pos_proj], string({exp.LM.name}))
title({decode.epoch_type;'over representation'})
pnl(2,3).select();
histogram(decode_flight.MAP_pos,'BinWidth',decode.params.pos_bin_size, 'Orientation','horizontal','Normalization','probability');
arrayfun(@(y,str)(yline(y,'Color',0.5.*[1 1 1],'LineWidth',0.5)),[exp.LM.pos_proj], string({exp.LM.name}))
title({decode_flight.epoch_type;'over representation'})
linkaxes(pnl(2).de.axis,'y')
pnl(1,1).select();
hold on
MUA_ts_IX = interp1(decode.time, 1:length(decode.time), exp.MUA.t, 'nearest');
% ripples_ts_IX = interp1(decode.time, 1:length(decode.time), exp.ripples.t, 'nearest');
PE_events_ts_IX = interp1(decode.time, 1:length(decode.time), [exp.PE.thr.peak_ts], 'nearest');
plot(MUA_ts_IX, exp.MUA.zFR);
% plot(ripples_ts_IX , exp.ripples.zpripple_all);
plot(PE_events_ts_IX, [exp.PE.thr.peak_zFR],'*r');
legend("Firing rate","Pop events")
xlim([0 1.05*length(decode.time)])
ylim([0 10])
rescale_plot_data('x',[1/decode.Fs 0]);
xlabel('Time (s)')
ylabel('Firing rate (z)')
linkaxes([pnl(1,1).axis pnl(2,1).axis],'x')

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
% saveas(fig, filename , 'tif');
saveas(fig, filename , 'jpg');
% for dpi = 100:50:300
%     exportgraphics(fig, sprintf('%s_%d_dpi.jpg',filename,dpi),'Resolution',dpi);
% end


%%

