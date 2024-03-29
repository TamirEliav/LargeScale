function decoding_plot_event(exp_ID, epoch_type, params_opt, event_type, event_num, win_s,clim_prc)
arguments
    %% temp for development...
    exp_ID = 'b9861_d180526';
    epoch_type = 'sleep';
    params_opt = 11;
    event_type = 'posterior';
%     event_type = 'PE';
    event_num = 49;
%     event_num = 113;
    win_s = 0.5;
    clim_prc = 1;
end

%% params
if ~exist('win_s','var')
    win_s = [];
end
switch numel(win_s)
    case 1
        win_s = [-1 1].*win_s;
    case 0
        win_s = [-1 1].*0.5;
end

%% folders
dir_OUT = 'F:\sequences\decoded_figs\event_examples';
figs_dir = fullfile(dir_OUT);
mkdir(figs_dir);

%% load data
exp = exp_load_data(exp_ID, 'details','path','rest','ripples','MUA','PE','pos');
decode = decoding_load_data(exp_ID, epoch_type, params_opt);
events = decoding_load_events(exp_ID, epoch_type, params_opt, event_type);

%% arrange data
% event to use
event = events(event_num);
% time window to display
t0 = event.peak_ts;
ti = t0 + win_s.*1e6;
% ti = ti - win_s(1).*1e6;

% get ripples/MUA
IX = get_data_in_ti(exp.ripples.t, ti);
zpripple = exp.ripples.zpripple_all(IX);
zpripple_t = exp.ripples.t(IX);
IX = get_data_in_ti(exp.MUA.t, ti);
zFR = exp.MUA.zFR(IX);
zFR_t = exp.MUA.t(IX);

% get states/position prob
IX = get_data_in_ti(decode.time, ti);
prob_t = decode.time(IX);
prob_state = decode.posterior_state(:,IX);
prob_pos = decode.posterior_pos(:,IX);
directions_strs = ["Outbound","Inbound"];
map1_IX = contains( decode.state, directions_strs(1));
map2_IX = contains( decode.state, directions_strs(2));
prob_pos_by_map = zeros([2 size(prob_pos)]);
prob_pos_by_map(1,:,:) = squeeze(sum(decode.posterior(:,map1_IX,IX),2));
prob_pos_by_map(2,:,:) = squeeze(sum(decode.posterior(:,map2_IX,IX),2));
real_pos = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos, prob_t);
if isempty(IX)
    error('data is missing!')
end

% % change time to ms aligned to event peak
% zpripple_t = (zpripple_t-t0) * 1e-3;
% zFR_t = (zFR_t-t0) * 1e-3;
% prob_t = (prob_t-t0) * 1e-3;

%% create figure
hf = figure;
figsize = [22 25];
hf.Units = 'centimeters';
hf.Position = [1 1 figsize];
hf.PaperUnits = 'centimeters';
hf.PaperPosition = [1 1 figsize];
pnl = panel();
pnl.margin = [20 10 65 15];
pnl.de.margin = 5;
pnl.pack('v',[.08 .08 .27 .27 .27]);
pnl.de.margin = 5;
pnl(3).marginbottom = 15;
pnl(4).marginbottom = 15;
% pnl(5).marginbottom = 15;

% plot MUA
pnl(1).select();
plot(zFR_t, zFR,'LineWidth',1.5);
yticks([0 ceil(max(zFR))])
% ylabel({'Multi unit activity';'firing rate (z)'})
ylabel({'MUA';'(z)'},'fontsize',14)
set(gca,'tickdir','out')
box on
hax=gca;
hax.XTickLabel = [];
hax.YLim(2) = ceil(max(zFR));
rescale_plot_data('x',[1e-3 ti(1)]);
xlim((win_s-win_s(1)).*1e3)

% plot state prob
pnl(2).select();
plot(prob_t, prob_state','LineWidth',1.5);
box on
hax=gca;
hax.YLim = [0 1];
hax.YTick = [0 1];
hax.TickDir = 'out';
hax.XTickLabel = [];
ylabel({'State';'prob.'},'fontsize',14)
rescale_plot_data('x',[1e-3 ti(1)]);
xlim((win_s-win_s(1)).*1e3)
% h=legend(decode.state,'NumColumns',round(length(decode.state)/3),'Location','southoutside');
h=legend(decode.state,'NumColumns',1,'Location','southoutside');
h.Position([1 2]) = [0.74 0.78];
h.Interpreter = 'none';
h.Box = 'on';

% plot position prob
pnl(3).select();
hold on
plot_prob_map(prob_t, decode.pos, prob_pos, clim_prc);
rescale_plot_data('x',[1e-3 ti(1)]);
xlim((win_s-win_s(1)).*1e3)
% plot position prob (per direction)
for ii_dir = 1:2
    pnl(3+ii_dir).select();
    plot_prob_map(prob_t, decode.pos, squeeze(prob_pos_by_map(ii_dir,:,:)), clim_prc);
    rescale_plot_data('x',[1e-3 ti(1)]);
    h = title(directions_strs(ii_dir));
    h.Units = 'normalized';
    h.Position = [-0.06 0.98];
    xlim((win_s-win_s(1)).*1e3)
end

%%
fig_name = sprintf('%s_%s_opt_%d_%s_event_%d_win_%d-%dms', exp_ID, epoch_type, params_opt, event_type, event_num, 1e3*abs(win_s));
h = pnl.title(fig_name);
h.Interpreter = 'none';
h.Position(2) = 1.03;
h.FontSize = 12;

params_str = {
sprintf('params opt: %d',params_opt),
sprintf('bin size: %.2gm',decode.params.pos_bin_size),
sprintf('replay speed: x%d',decode.params.replay_speed),
sprintf('state_decay_timescale: %.3g s',decode.params.state_decay_timescale),
sprintf('clim prc: %d',clim_prc),
};
annotation('textbox', [0.8 0.2 0.2 0.05], 'String',params_str,'LineStyle','None','Interpreter','none');

linkaxes( [pnl.de.axis(:)] ,'x');
linkaxes( [pnl.de.axis([3 4 5])] ,'y');

% save fig
filename = fullfile(figs_dir, fig_name);
saveas(gcf, filename , 'jpg');
saveas(gcf, filename , 'pdf');
% close(gcf)



end


%%
function plot_prob_map(t, pos, prob, clim_prc)
    hold on    
    imagesc(t, pos, prob);
    axis tight
%     if strcmp(epoch_type,'rest')
%         plot(t, real_pos, 'r', 'LineWidth',0.01)
%     end
    hax = gca;
    hax.TickDir = 'out';
    hax.TickLength(1) = 0.008;
    % climits = quantile(prob_pos(:),[0.01 0.99]);
    climits = prctile(prob(:), [0 100]+[1 -1].*clim_prc);
    %         climits  = [0 10/length(decode.pos)];
    if all(~isnan(climits))
        hax.CLim = climits;
    else
        error('check what is the problem')
    end
    hax.XLim = t([1 end]);
    %     hax.YLim = pos([1 end])  + [-1;1].*10;
    % hax.YLim = [min(decode.pos) max(decode.pos)] + [-1 1].*median(diff(decode.pos))*3;
    box on
    hax.XRuler.TickLabelGapOffset = -4;
    xlabel('Time (ms)','fontsize',14)
    ylabel('Position (m)','fontsize',14)
    yticks(0:50:200)
%     switch exp.details.recordingArena
%         case '200m'
%             ylim([0 200])
%         case '120m'
%             ylim([0 150])
%         otherwise
%                 error('wrong recording arena')
%     end
    cmap = bone;
    cmap = flipud(cmap);
    colormap(cmap);
end

