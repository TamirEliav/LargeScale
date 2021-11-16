function decoding_plot_flight_posterior(exp_ID, params_opt)

%% temp for development...
% exp_ID = 'b0184_d191201';
% exp_ID = 'b9861_d180527';
% params_opt = 4;
epoch_type = 'flight';

%% load data
exp = exp_load_data(exp_ID, 'details','path','rest','ripples','MUA','PE','pos','flight');
dir_IN = 'F:\sequences\decoded';
dir_OUT = 'F:\sequences\decoded_figs';
figs_dir = fullfile(dir_OUT, epoch_type, exp_ID, "opt_"+params_opt);
decode_filename = fullfile(dir_IN, epoch_type, exp_ID, sprintf('%s_%s_opt_%d.nc',exp_ID,epoch_type,params_opt));
decode = decoding_read_decoded_file(decode_filename);
mkdir(figs_dir);

%% FE to use
FE_to_use = exp.flight.FE;
FE_to_use([FE_to_use.distance]<100) = [];
[FE_to_use.num] = disperse(1:length(FE_to_use));

%%
nRows = 5;
nCols = 5;
nPanels = nRows * nCols;
nFigs = ceil(length(FE_to_use) / nPanels);

for ii_fig = 1:nFigs
    
    %%
    FE_IX_start = (ii_fig-1)*nPanels +1;
    FE_IX_end = min(length(FE_to_use), FE_IX_start+nPanels-1);
    FE_IX = FE_IX_start : FE_IX_end;
    FE = FE_to_use(FE_IX);
    hf = figure;
    hf.Units = 'centimeters';
    hf.WindowState = 'maximized';
    pnl = panel();
    pnl.pack(nRows,nCols);
    pnl.margin = [25 25 10 20];
    pnl.de.margin = [10 10 10 10];
    for r = 1:nRows
        for c = 1:nCols
            pnl(r,c).pack('v',[.2 .2 .6]);
            pnl(r,c).de.margin = 1;
        end
    end
    for ii_FE = 1:length(FE)
        r = ceil(ii_FE/(nRows));
        c = mod(ii_FE-1, nCols)+1;

        % time window to display
        t0 = FE(ii_FE).start_ts;
        ti = [FE(ii_FE).start_ts FE(ii_FE).end_ts];

        % get ripples/MUA
%         IX = get_data_in_ti(exp.ripples.t, ti);
%         zpripple = exp.ripples.zpripple_all(IX);
%         zpripple_t = exp.ripples.t(IX);
        IX = get_data_in_ti(exp.MUA.t, ti);
        zFR = exp.MUA.zFR(IX);
        FR = exp.MUA.FR(IX);
        zFR_t = exp.MUA.t(IX);

        % get states/position prob
        IX = get_data_in_ti(decode.time, ti);
        prob_state = decode.posterior_state(:,IX);
        prob_pos = decode.posterior_pos(:,IX);
        prob_t = decode.time(IX);
        real_pos = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos, prob_t);
        if isempty(IX)
            continue;
        end

        % change time to ms aligned to flight peak
%         zpripple_t = (zpripple_t-t0) * 1e-6;
        zFR_t = (zFR_t-t0) * 1e-6;
        prob_t = (prob_t-t0) * 1e-6;

        % plot ripple power / MUA
        pnl(r,c,1).select();
        title("flight #" + FE(ii_FE).num)
%         yyaxis left
        plot(zFR_t, FR,'LineWidth',1.5);
%         yticks([0 ceil(max(FR))])
        set(gca,'tickdir','out')
%         yyaxis right
%         plot(zFR_t, zFR,'LineWidth',1.5);
%         yticks([0 ceil(max(zFR))])
%         set(gca,'tickdir','out')
        xlim(prob_t([1 end]))
        box on
        hax=gca;
        hax.XTickLabel = [];

        % plot state prob
        pnl(r,c,2).select();
        plot(prob_t, prob_state','LineWidth',1.5);
        box on
        hax=gca;
        hax.XLim = prob_t([1 end]);
        hax.YLim = [0 1];
        hax.YTick = [0 1];
        hax.TickDir = 'out';
        hax.XTickLabel = [];

        % plot position prob
        pnl(r,c,3).select();
        hold on
        imagesc(prob_t, decode.pos, prob_pos);
        axis tight
        plot(prob_t, real_pos, 'r', 'LineWidth',0.01)
        hax = gca;
        hax.TickDir = 'out';
        hax.TickLength(1) = 0.008;
        climits = quantile(prob_pos(:),[0.01 0.99]);
%         climits  = [0 10/length(decode.pos)];
        if all(~isnan(climits))
            hax.CLim = climits;
        end
        hax.XLim = prob_t([1 end]);
    %     hax.YLim = pos([1 end])  + [-1;1].*10;
        hax.YLim = [min(decode.pos) max(decode.pos)] + [-1 1].*median(diff(decode.pos))*3;
        box on
        hax.XRuler.TickLabelGapOffset = -4;

    end
    cmap = bone;
    cmap = flipud(cmap);
    colormap(cmap);

    %% add title/labels/legend/text
    pnl.margin = [15 22 12 12];
    
    fig_name = sprintf('%s_flights_%d-%d', exp_ID, FE_IX_start, FE_IX_end );
    h = pnl.title(fig_name);
    h.Interpreter = 'none';
    h.Position(2) = 1.02;
    h.FontSize = 14;

    h=pnl.xlabel('Time (ms)');
    h.FontSize = 12;
%     pnl(1,1,1).select();
%     yyaxis left
    h=pnl(1,1,1).ylabel('FR (Hz)');
    h.FontSize = 8;
    h=pnl(1,1,2).ylabel({'state';'prob.'});
    h.FontSize = 8;
    h=pnl(1,1,3).ylabel('Position (m)');
    h.FontSize = 8;
    
    pnl(1,1,1).select();
%     h=legend({'MUA firing rate (Hz)','MUA firing rate (z)'},'NumColumns',1,'Location','southoutside');
    h=legend({'MUA firing rate (Hz)'},'NumColumns',1,'Location','southoutside');
    h.Position([1 2]) = [0.25 0.01];
    h.Interpreter = 'none';
    h.Box = 'on';
    
    pnl(1,1,2).select();
    h=legend(decode.state,'NumColumns',round(length(decode.state)/2),'Location','southoutside');
    h.Position([1 2]) = [0.025 0.01];
    h.Interpreter = 'none';
    h.Box = 'on';
    
    params_str = {
        sprintf('params opt: %d',params_opt),
        sprintf('bin size: %.2gm',decode.params.pos_bin_size),
        sprintf('replay speed: x%d',decode.params.replay_speed),
        sprintf('state_decay_timescale: %.3g s',decode.params.state_decay_timescale),
        };
    annotation('textbox', [0.65 0.01 0.2 0.05], 'String',params_str,'LineStyle','None')
    
    %% save fig
    filename = fullfile(figs_dir, fig_name);
    saveas(gcf, filename , 'jpg');
    close(gcf)
    
end



%%





end