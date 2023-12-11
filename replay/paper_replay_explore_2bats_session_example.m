%%
clear
clc
close all

%%
win_s = 60*5; % 5 min +- window
remove_near_balls = true;

%%
data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115.mat';
data = load(data_filename);
exp_ID_list = data.data_info.exp_ID(data.TF);
ts_list = data.data_info.ts(data.TF);
events_num = data.data_info.evnet_num(data.TF);

%%
fig=figure;
fig.WindowState = 'maximized';
tiledlayout(8,6,'TileSpacing','loose');

%%
for ii_ex = 1:length(exp_ID_list)
    %%
    nexttile
    cla
    hold on

    %% load data
    exp_ID = exp_ID_list{ii_ex};
    exp = exp_load_data(exp_ID, 'details', 'pos');
    epoch_type = 'rest';
    params_opt = 11;
    event_type = 'posterior';
    [events, params]= decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
    [seqs, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
    events(~TF)=[];
    session_ti = exp_get_sessions_ti(exp_ID,'Behave');
    ti = ts_list(ii_ex)+[-1 1].*win_s*1e6;
    lw = 2;

    %% plot
    if remove_near_balls
        replay_pos_limits = [25 115];
        TF = [seqs.middle_pos] > replay_pos_limits(1) & [seqs.middle_pos] < replay_pos_limits(2);
        seqs = seqs(TF);
    end
    plot(exp.pos.proc_1D.ts, interp_nans(exp.pos.proc_1D.pos),'LineWidth',lw);
    plot(exp.pos.proc_1D.other.ts, interp_nans(exp.pos.proc_1D.other.pos(1,:)),'LineWidth',1);
    plot(exp.pos.proc_1D.co.ts,exp.pos.proc_1D.co.pos,'xk','MarkerSize',8)
    plot([seqs.start_ts;seqs.end_ts],[seqs.start_pos; seqs.end_pos],'-m','LineWidth',1.3);
    % plot([seqs.start_ts],[seqs.start_pos],'.m','MarkerSize',10)
    plot([seqs([seqs.direction]== 1).end_ts],[seqs([seqs.direction]== 1).end_pos],'^m','MarkerSize',2,'MarkerFaceColor','m');
    plot([seqs([seqs.direction]==-1).end_ts],[seqs([seqs.direction]==-1).end_pos],'vm','MarkerSize',2,'MarkerFaceColor','m');
    xline(35969648929.0)
    xlim(ti)
    rescale_plot_data('x',[1e-6/60 ti(1)]);
    ylim([0 135])
    % xticks(linspace(0,135,4))
    yticks(linspace(0,135,4))
    xlabel('Time (min)', 'Units','normalized', 'Position',[0.5 -0.11]);
    ylabel('Position (m)', 'Units','normalized', 'Position',[-0.13 .5]);
    hax=gca;
    hax.XRuler.TickLabelGapOffset = -1.8;
    hax.YRuler.TickLabelGapOffset = 1;
    str = sprintf('%s_%d',exp_ID,events_num(ii_ex));
    title(str,'Interpreter','none');

end

%%
filename = 'fig_3b_options';
switch remove_near_balls
    case true
        filename = [filename '_exc_near_balls'];
    case false
        filename = [filename '_inc_near_balls'];
end
exportgraphics(fig,[filename '.pdf'])





