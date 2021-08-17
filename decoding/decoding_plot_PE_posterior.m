function decoding_plot_PE_posterior(exp_ID, win_s, replay_speed, pos_BW, pos_binsize)

%%
% win_s = 0.5;
% exp_ID = 'b9861_d180527';
% replay_speed = 200;
% pos_BW = 2;
% pos_binsize = 1;
dir_IN = 'D:\sequences\seq_uri_eden\decoded';
dir_OUT = 'D:\sequences\seq_uri_eden\decoded_figs';
h5_filename = fullfile(dir_IN, sprintf('%s_decoding_sleep_replay_speed_%d_pos_BW_%g_binsize_%gm.h5',exp_ID,replay_speed,pos_BW,pos_binsize));
% filename = "D:\sequences\seq_uri_eden\decoded\b9861_d180526_sleep_decoding.h5";
exp = exp_load_data(exp_ID, 'details','path','PE','MUA','ripples');

%% arrange PE to use
PE_to_use = exp.PE.strong;
[PE_to_use.num] = disperse(1:length(PE_to_use));

%% read
pos = h5read(h5_filename,'/position');
posterior = h5read(h5_filename,'/acausal_posterior');
state = h5read(h5_filename,'/state');
time = h5read(h5_filename,'/time');
nPos = length(pos);
nStates = length(state);
nTimes = length(time);
dim_pos = find(size(posterior) == nPos);
dim_state = find(size(posterior) == nStates);
dim_time = find(size(posterior) == nTimes);

%% calc posterior marginals and MAP
posterior_state = squeeze(sum(posterior,dim_pos));
posterior_pos = squeeze(sum(posterior,dim_state));
[~,MAP_pos_IX] = max(posterior_pos,[],1);
[~,MAP_state_IX] = max(posterior_state,[],1);
MAP_pos = pos(MAP_pos_IX);
MAP_state = state(MAP_state_IX);

%%
nRows = 5;
nCols = 5;
nPanels = nRows * nCols;
nFigs = ceil(length(PE_to_use) / nPanels);
for ii_fig = 1:nFigs
    
    %%
    PE_IX_start = (ii_fig-1)*nPanels +1;
    PE_IX_end = min(length(PE_to_use), PE_IX_start+nPanels-1);
    PE_IX = PE_IX_start : PE_IX_end;
    PE = PE_to_use(PE_IX);
    hf = figure;
    hf.Units = 'centimeters';
    hf.WindowState = 'maximized';
    pnl = panel();
    pnl.pack(nRows,nCols);
    pnl.margin = [25 20 15 20];
    pnl.de.margin = [10 10 10 10];
    for r = 1:nRows
        for c = 1:nCols
            pnl(r,c).pack('v',[.2 .2 .6]);
            pnl(r,c).de.margin = 1;
        end
    end
    for ii_PE = 1:length(PE)
        r = ceil(ii_PE/(nRows));
        c = mod(ii_PE-1, nCols)+1;

        % time window to display
        t0 = PE(ii_PE).peak_ts;
        ti = t0 + [-1 1].*win_s*1e6;

        % get ripples/MUA
        IX = get_data_in_ti(exp.ripples.t, ti);
        zpripple = exp.ripples.zpripple_all(IX);
        zpripple_t = exp.ripples.t(IX);
        IX = get_data_in_ti(exp.MUA.t, ti);
        zFR = exp.MUA.zFR(IX);
        zFR_t = exp.MUA.t(IX);

        % get states/position prob
        IX = get_data_in_ti(time', ti);
        prob_state = posterior_state(:,IX);
        prob_pos = posterior_pos(:,IX);
        prob_t = time(IX);
        if isempty(IX)
            continue;
        end

        % change time to ms aligned to event peak
        zpripple_t = (zpripple_t-t0) * 1e-3;
        zFR_t = (zFR_t-t0) * 1e-3;
        prob_t = (prob_t-t0) * 1e-3;

        % plot ripple power / MUA
        pnl(r,c,1).select();
        title("event #" + PE(ii_PE).num)
        yyaxis left
        plot(zpripple_t, zpripple,'LineWidth',1.5);
        yticks([0 floor(max(zpripple))])
        set(gca,'tickdir','out')
        yyaxis right
        plot(zFR_t, zFR,'LineWidth',1.5);
        yticks([0 floor(max(zFR))])
        set(gca,'tickdir','out')
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
        imagesc(prob_t, pos, prob_pos);
        axis tight
        hax = gca;
        hax.TickDir = 'out';
        hax.TickLength(1) = 0.008;
    % 	hax.CLim = quantile(prob_pos(:),[0.05 0.95]);
        hax.CLim = [0 10/length(pos)];
        hax.XLim = prob_t([1 end]);
    %     hax.YLim = pos([1 end])  + [-1;1].*10;
        hax.YLim = [0 200];
        box on
        hax.XRuler.TickLabelGapOffset = -4;

    end
    cmap = bone;
    cmap = flipud(cmap);
    colormap(cmap);

    %% add title/labels/legend
    pnl.marginbottom = 12;
    pnl.margin = [15 12 12 18];
    
    fig_name = sprintf('%s_win_%dms_events_%d-%d', exp_ID, win_s*1e3, PE_IX_start, PE_IX_end );
    h = pnl.title(fig_name);
    h.Interpreter = 'none';
    h.Position(2) = 1.05;
    h.FontSize = 14;

    h=pnl.xlabel('Time (ms)');
    h.FontSize = 12;
    pnl(1,1,1).select();
    yyaxis left
    h=pnl(1,1,1).ylabel('(z)');
    h.FontSize = 8;
    h=pnl(1,1,2).ylabel({'state';'prob.'});
    h.FontSize = 8;
    h=pnl(1,1,3).ylabel('Position (m)');
    h.FontSize = 8;
    
    pnl(1,1,1).select();
    h=legend({'Ripple power (z)','MUA firing rate (z)'},'NumColumns',1,'Location','northeastoutside');
    h.Position([1 2]) = [0.25 0.96];
    h.Interpreter = 'none';
    h.Box = 'on';
    
    pnl(1,1,2).select();
    h=legend(state,'NumColumns',round(length(state)/2),'Location','northeastoutside');
    h.Position([1 2]) = [0.025 0.959];
    h.Interpreter = 'none';
    h.Box = 'on';
    
    %% save fig
    folder = fullfile(dir_OUT, exp_ID, sprintf('sleep_replay_speed=%d_pos_BW=%g_binsize=%g',replay_speed,pos_BW,pos_binsize));
    mkdir(folder);
    filename = fullfile(folder, fig_name );
    saveas(gcf, filename , 'jpg');
    close(gcf)
    
end



%%





end