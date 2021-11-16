function decoding_plot_rest_posterior(exp_ID, opt_params)

%%
exp_ID = 'b0184_d191201';
opt_params = 4;

%% load data
exp = exp_load_data(exp_ID, 'details','path','rest','ripples','MUA','PE','pos');
dir_IN = 'F:\sequences\decoded';
dir_OUT = 'F:\sequences\decoded_figs';
figs_dir = fullfile(dir_OUT,'rest',exp_ID,"opt_"+opt_params);
decode_filename = fullfile(dir_IN,'rest',exp_ID, sprintf('%s_rest_opt_%d.nc',exp_ID,opt_params));
decode = decoding_read_decoded_file(decode_filename);
mkdir(figs_dir);

%% arrange events ti to plot
PE_ti = [exp.PE.thr.start_ts;exp.PE.thr.end_ts]';
rest_ti = exp.rest.ti;
IX = any(PE_ti>shiftdim(rest_ti(:,1),-2) & PE_ti<shiftdim(rest_ti(:,2),-2), [2 3]);
events_ti = PE_ti(IX,:);

%%
for ii_event = 1:size(events_ti,1)
    
    %% create figure
    hf = figure;
    hf.WindowState = 'maximized';
    hf.Units = 'centimeters';
    hf.Position = [5 2 30 24];
    pnl = panel();
    pnl.pack('h',[0.8 0.2]);
    pnl(1).pack('v',[0.1 0.1 0.1 0.7]);
    
    fig_name = sprintf('%s_rest_%d', exp_ID, ii_event);
    h = pnl.title(fig_name);
    h.Interpreter = 'none';
    h.Position(2) = 0.95;
    h.FontSize = 14;
    
    % arrange data
    ti = events_ti(ii_event,:);
    t0 = ti(1);
    
    % get ripples/MUA
    IX = get_data_in_ti(exp.ripples.t, ti);
    zpripple = exp.ripples.zpripple_all(IX);
    zpripple_t = exp.ripples.t(IX);
    IX = get_data_in_ti(exp.MUA.t, ti);
    zFR = exp.MUA.zFR(IX);
    zFR_t = exp.MUA.t(IX);

    % get states/position prob
    IX = get_data_in_ti(decode.time', ti);
    prob_state = decode.posterior_state(:,IX);
    prob_pos = decode.posterior_pos(:,IX);
    prob_t = decode.time(IX);
%     if isempty(IX)
%         continue;
%     end
    
    % change time to ms aligned to rest start ts
    zpripple_t = (zpripple_t-t0) * 1e-6;
    zFR_t = (zFR_t-t0) * 1e-6;
    prob_t = (prob_t-t0) * 1e-6;

    % plot ripple power / MUA
%     pnl(2).select();
%     yyaxis left
%     plot(zpripple_t, zpripple,'LineWidth',1.5);
%     yticks([0 floor(max(zpripple))])
%     set(gca,'tickdir','out')
%     yyaxis right
%     plot(zFR_t, zFR,'LineWidth',1.5);
%     yticks([0 floor(max(zFR))])
%     set(gca,'tickdir','out')
%     xlim(prob_t([1 end]))
%     box on
%     hax=gca;
%     hax.XTickLabel = [];

    % plot state prob
    pnl(1,3).select();
    plot(prob_t, prob_state','LineWidth',1.5);
    box on
    hax=gca;
    hax.XLim = prob_t([1 end]);
    hax.YLim = [0 1];
    hax.YTick = [0 1];
    hax.TickDir = 'out';
    hax.XTickLabel = [];
    h=legend(decode.state,'NumColumns',round(length(decode.state)/2),'Location','northeastoutside');
    h.Position([1 2]) = [0.025 0.85];
    h.Interpreter = 'none';
    h.Box = 'on';
    ylabel('Probability')
    
    % plot position prob
    pnl(1,4).select();
    imagesc(prob_t, decode.pos, prob_pos);
    axis tight
    hax = gca;
    hax.TickDir = 'out';
    hax.TickLength(1) = 0.008;
% 	hax.CLim = quantile(prob_pos(:),[0.05 0.95]);
%     hax.CLim = [0 10/length(decode.pos)];
    hax.XLim = prob_t([1 end]);
%     hax.YLim = pos([1 end])  + [-1;1].*10;
    hax.YLim = [min(decode.pos) max(decode.pos)] + [-1 1];
    box on
    hax.XRuler.TickLabelGapOffset = -4;
    cmap = bone;
    cmap = flipud(cmap);
    colormap(cmap);
    xlabel('Time (s)')
    ylabel('Position (m)')
    
    % plot behavior and mark rest period
    pnl(2).select();
    hold on
    plot(exp.pos.proc_1D.pos, exp.pos.proc_1D.ts)
    hax = gca;
    xlimits = hax.XLim;
%     x1 = xlimits(1);
%     x2 = xlimits(2);
%     y1 = ti(1);
%     y2 = ti(2);
%     h = patch([x1 x1 x2 x2],[y1 y2 y2 y1],'r');
%     h.FaceAlpha = 0.3;
%     h.EdgeColor = 'none';
    clear hl
    h(1)=yline(ti(1));
    h(2)=yline(ti(2));
    [h.Color] = disperse(repelem('r',length(h)));
    [h.LineWidth] = disperse(repelem(0.5,length(h)));
    behave_ti = exp_get_sessions_ti(exp_ID,'Behave');
    ylim(behave_ti);
    rescale_plot_data('y',[1e-6/60 behave_ti(1)])
    xlabel('Position (m)')
    ylabel('Time (min)')
    
    linkaxes([pnl(1,3).axis pnl(1,4).axis],'x');
    
    %% save fig
    filename = fullfile(figs_dir, fig_name );
    saveas(gcf, filename , 'jpg');
    close(gcf)
end



%%



end