%% Replay - Fig supp 2 - replay example deconstructed
%%
clear 
clc
close all

%% plotting options
% example_opt_IX = [5 5 12] % I chose option 5, which is example #386 in
% the new list
% example_options =[
% struct(exp_ID = 'b0184_d191130', epoch_type = 'sleep', params_opt = 11, event_num = 139)
% struct(exp_ID = 'b0184_d191201', epoch_type = 'sleep', params_opt = 11, event_num = 100)
% struct(exp_ID = 'b0184_d191202', epoch_type = 'sleep', params_opt = 11, event_num = 204)
% struct(exp_ID = 'b0184_d191202', epoch_type = 'sleep', params_opt = 11, event_num = 60)
% struct(exp_ID = 'b0184_d191203', epoch_type = 'sleep', params_opt = 11, event_num = 31)
% struct(exp_ID = 'b0184_d191203', epoch_type = 'sleep', params_opt = 11, event_num = 34)
% struct(exp_ID = 'b0184_d191204', epoch_type = 'sleep', params_opt = 11, event_num = 45)
% struct(exp_ID = 'b0184_d191205', epoch_type = 'sleep', params_opt = 11, event_num = 227)
% struct(exp_ID = 'b0184_d191212', epoch_type = 'sleep', params_opt = 11, event_num = 56)
% struct(exp_ID = 'b2382_d190712', epoch_type = 'sleep', params_opt = 11, event_num = 24)
% struct(exp_ID = 'b2382_d190712', epoch_type = 'sleep', params_opt = 11, event_num = 24)
% struct(exp_ID = 'b2299_d191213', epoch_type = 'rest',  params_opt = 11, event_num = 20) % bats crossover
% ];
% examples_list = example_options(example_opt_IX);

%% replay examples - new format
replay_examples_list_filename = "L:\Analysis\Code\inclusion_lists\replay_examples.xlsx";
replay_examples_list = table2struct(readtable(replay_examples_list_filename));
% replay_examples_options = [
%     536 247 1000;
%     
%     516 499 483 % main 2
%     354 536 552 % main 3   
%     193 197 131 % main 4
%     29 247 160  % main 5
% 
%     500 328 348 % supp 6
%     555 359 378 % supp 7
%     540 394 543 % supp 8
%     460 557 450 % supp 9
%     561 339 463 % supp 10
%     364 493 288 % supp 11
%     168  58  76 % supp 12
%     15   39 112 % supp 13
%     177 109 113 % supp 14
%     62  51  41  % supp 15
%     238 228 216 % supp 16
%     257 253  83 % supp 17
% ];
% replay_ex_opt = 17; 
replay_examples_options = [
536 552 493 450 41 76 253 1000;
450 493 552 536 41 76 253 1000;
493 450 552 536 41 76 253 1000;
];
replay_ex_opt = 3; 
[~,IX] = ismember(replay_examples_options(replay_ex_opt,:),[replay_examples_list.ex_num]);
replay_examples = replay_examples_list(IX);
disp('chosen examples:')
for ii_ex = 1:length(replay_examples)
    ex = replay_examples(ii_ex);
    fprintf('%d_%s_%s_%d\n',ex.ex_num,ex.epoch_type,ex.exp_ID,ex.event_num)
end

titles_str = {
    'Sleep replay',
    'Sleep replay',
    'Sleep replay',
    'Sleep replay',
    'Awake replay',
    'Awake replay',
    'Awake replay',
    'Awake replay, 2-bats',
    };

%% graphics params
win_s = 0.5;
cmap = bone;
cmap = flipud(cmap);
xlimits = [-1 1]*.3;
xtick = [-0.3:0.1:0.3];
states_clrs = [
133  63 0
160 225 157
235 223 99
0 98 190
230 136 180
0 251 246
]./255;

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Figure S2';
fig_caption_str = ' ';
log_name_str = [fig_name_str '_log_file' '.txt'];
log_name_str = strrep(log_name_str , ':', '-');
log_name_str = strrep(log_name_str , ' ', '_');
log_name_out = fullfile(res_dir, log_name_str);

%% open log file
diary off
diary(log_name_out)
diary on
disp('Log file');
disp(['created: ', datestr(clock)]);
disp('======================================================');
disp([fig_name_str ':' fig_caption_str]);   
disp('======================================================');
disp('');

%% create figure
% figure_size_cm = [21.0 29.7]; % ~A4
figure_size_cm = [21.6 27.9]; % ~US letter
fig = figure;
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',7);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'DefaultAxesUnits','centimeters');
set(gcf,'PaperType','usletter')
% set(gcf,'PaperType','<custom>');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 figure_size_cm]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]); % position on screen...
set(gcf, 'Renderer', 'painters');
% set(gcf, 'color', 'none');
set(groot, 'defaultAxesColor','None')
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

% create panels
offsets_x = linspace(2,17,4);
offsets_y = linspace(3,14,2)+0.2;
offsets_y = flip(offsets_y);
clear panels_ex panels
for ii=1:2
    for jj=1:4
        offset_x = offsets_x(jj);
        offset_y = offsets_y(ii);
        offset = [offset_x offset_y];
        W = 4;
        H1 = 3.5;
        H2 = 0.9;
        y = 0;
        H=0;
        ymarg = 0.2;
        y=y+H;       H=H1; panels_ex{jj,ii}(1) = axes('position', [offset+[0 y] W H]);
        y=y+H+ymarg; H=H2; panels_ex{jj,ii}(2) = axes('position', [offset+[0 y] W H]);
        y=y+H+ymarg; H=H2; panels_ex{jj,ii}(3) = axes('position', [offset+[0 y] W H]);
        y=y+H+ymarg; H=H2; panels_ex{jj,ii}(4) = axes('position', [offset+[0 y] W H]);
        y=y+H+ymarg; H=H2; panels_ex{jj,ii}(5) = axes('position', [offset+[0 y] W H]);
    end
end
panels_ex = {panels_ex{:}}';
% offset_y = 22;
% panels{1}(1) = axes('position', [2 offset_y 3 3]);
% panels{1}(2) = axes('position', [6 offset_y 3 3]);
% panels{1}(3) = axes('position', [10 offset_y 3 3]);
% panels{1}(4) = axes('position', [14 offset_y 3 3]);

%%
for ii_ex = 1:length(replay_examples)
    %% load data
    ex = replay_examples(ii_ex);
    addFieldsToWorkspace(ex);
    decode = decoding_load_data(exp_ID, epoch_type, params_opt );
    exp = exp_load_data(exp_ID,'details','path','MUA','ripples');
    events = decoding_load_events_quantification(exp_ID,epoch_type,params_opt,"posterior");
    event = events([events.num] ==event_num);
    seq = event.seq_model;
    seq_ti = [event.start_ts event.end_ts];
    t0 = mean(seq_ti);
    ti = t0 + [-1 1].*win_s*1e6;
    
    TT = exp.ripples.stats.best_TT;
    [LFP.signal, LFP.ts, LFP.fs, LFP.params] = LFP_load(exp_ID,TT,'limits_ts',ti);
    [LFP_ripples.signal, LFP_ripples.ts, LFP_ripples.fs, LFP_ripples.params] = LFP_load(exp_ID,TT,'band','ripple','limits_ts',ti);
    LFP.avg_signal = nanmean(LFP.signal,[2 3]);
    LFP_ripples.avg_signal = nanmean(LFP_ripples.signal,[2 3]);
    
    %% rename decoder states
    states = decode.state';
    states(states=="Inbound-empirical_movement") = "Movement state dir 2";
    states(states=="Inbound-identity") = "Stationary state dir 2";
    states(states=="Inbound-uniform") = "Discontinuous state dir 2";
    states(states=="Outbound-empirical_movement") = "Movement state dir 1";
    states(states=="Outbound-identity") = "Stationary state dir 1";
    states(states=="Outbound-uniform") = "Discontinuous state dir 1";
    
    %% legend
    if ii_ex==1
        if exist('panels_legend','var')
            delete(panels_legend);
        end
        panels_legend = axes('position', [13.5 24.6 6 0.8]);
        cla
        hold on
        lw=2;
        t = [0 0.05];
        x = [.62 .62 .62 0 0 0];
        y = [1 .5 0 1 .5 0];
        hax = gca;
        hax.ColorOrder = states_clrs;
        for ii=1:size(states_clrs,1)
            plot(x(ii)+t,y([ii ii]),'-','LineWidth',lw,'Clipping','off');
            text(x(ii)+t(2)+0.02, y(ii), states(ii),'FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
        end
        ha = annotation("arrow");
        ha.Parent = panels_legend;
        ha.HeadStyle = 'vback1';
        ha.HeadLength = 5;
        ha.HeadWidth = 5;
        ha.Units = 'normalized';
        ha.X = [0 0]-.05;
        ha.Y = [0 1];
        ha = annotation("arrow");
        ha.Parent = panels_legend;
        ha.HeadStyle = 'vback1';
        ha.HeadLength = 5;
        ha.HeadWidth = 5;
        ha.Units = 'normalized';
        ha.X = [0 0]+.58;
        ha.Y = [1 0];
        hax.Clipping = 'off';

        x = [.0 .2];
        y = -0.7;
        plot(x, [1 1].*y, 'r--','Clipping','off')
        text(x(2)+0.02, y, 'Threshold probability = 0.8','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
        

        xlim([0 1])
        ylim([0 1])
        axis off

        panel = panels_legend;
        colormap(cmap);
        hcb = colorbar('southoutside');
        hcb.Units = 'centimeters';
        hcb.Position = [18.8 23.9 1 0.25];
        hcb.Label.Rotation = 0;
%         hcb.Label.Position(1) = 1.5;
        hcb.Label.String = 'Probability';
        hcb.Label.Position(2) = -0.2;
        hcb.Ticks = [];
        text(0.85, y, '0','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
        text(1.07, y, 'Max','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
%         text(cb_offset_x_middle,1,'Max','Units','normalized','FontSize',7,'HorizontalAlignment','center');
%         text(cb_offset_x_middle,0,'0','Units','normalized','FontSize',7,'HorizontalAlignment','center');
    end

    %% plot MUA
    axes(panels_ex{ii_ex}(5));
    cla
    hold on
    IX = get_data_in_ti(exp.MUA.t,ti);
    x = exp.MUA.t(IX);
    y = exp.MUA.FR(IX);
    area(x,y,'FaceColor','k');
    xlim(ti)
    hax=gca;
    m = roundToNearestValue(max(y), 20, @ceil)
    hax.XTick = [];
    hax.YTick = [0 m];
    hax.YLim = [0 m]
    rescale_plot_data('x',[1e-6 t0]);
    % axis off
    if ismember(ii_ex,[1 5])
        ylabel({'Multiunit';'firing-rate (Hz)'}, 'Units','normalized', 'Position',[-0.2 .42]);
    end
    title_ypos = [1.15 1.15 1.15 1.15 1.3 1.3 1.3 1.35];
    text(0.5, title_ypos(ii_ex), titles_str{ii_ex}, ...
        'Units','normalized','FontWeight','normal','FontSize',8,'HorizontalAlignment','center');
%     text(0.5, 1.45, sprintf('%d_%s_%s_%d',ex_num,epoch_type,exp_ID,event_num),'Interpreter','none', ...
%         'Units','normalized','FontWeight','normal','FontSize',6,'HorizontalAlignment','center');
    
    %% plot LFP (ripple-band)
    axes(panels_ex{ii_ex}(4));
    cla
    hold on
    plot(LFP_ripples.ts, LFP_ripples.avg_signal,'k');
    xlim(ti)
    xticks([])
    yticks([])
    rescale_plot_data('x',[1e-6 t0]);
    axis off
    if ismember(ii_ex,[1 5])
        text(-0.06, .5, 'Ripple', 'Units','normalized','FontSize',8.25,'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    
    %% plot LFP (raw)
    axes(panels_ex{ii_ex}(3));
    cla
    hold on
    plot(LFP.ts, LFP.avg_signal,'k');
    xlim(ti)
    xticks([])
    yticks([])
    rescale_plot_data('x',[1e-6 t0]);
    axis off
    if ismember(ii_ex,[1 5])
        text(-0.06, .5, 'SWR', 'Units','normalized','FontSize',8.25,'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
    end

    %% plot posterior (state)
    axes(panels_ex{ii_ex}(2));
    cla
    hold on
    hax = gca;
    hax.ColorOrder = states_clrs;
    IX = get_data_in_ti(decode.time, ti);
    prob_t = decode.time(IX);
    h=plot(decode.time(IX), decode.posterior_state(:,IX),'LineWidth',1);
    h(event.state_num).LineWidth=2.5; % highlight the event state
    yline(0.8,'r--')
    hax = gca;
    hax.XLim = prob_t([1 end]);
    hax.YLim = [0 1];
    box off
    colormap(cmap);
    hax.XLim = seq_ti+[-1 1].*0.2*range(seq_ti);
    hax.XTick = xtick;
    hax.XTickLabel=[];
    hax.YTick = [0 1];
    hax.TickDir = 'out';
    hax.TickLength = [0.02 0.02];
    hax.XRuler.TickLabelGapOffset = -4;
    hax.XRuler.TickLabelGapOffset = -1;
    rescale_plot_data('x',[1e-6 t0]);
    if ismember(ii_ex,[1 5])
        ylabel({'State';'prob.'}, 'Units','normalized', 'Position',[-0.12 .5]);
    end
    
    %% plot posterior (position)
    axes(panels_ex{ii_ex}(1));
    cla
    hold on
    IX = get_data_in_ti(decode.time, ti);
    prob_t = decode.time(IX);
    prob_pos = squeeze(decode.posterior(:,event.state_num,IX));
    imagesc(prob_t, decode.pos, prob_pos);
    plot([seq.start_ts seq.end_ts],[seq.start_pos seq.end_pos],'-r','LineWidth',1);
    hax = gca;
    clim_prctiles = [1 99];
    hax.CLim = prctile(prob_pos(:),clim_prctiles);
    hax.XLim = prob_t([1 end]);
    hax.YLim = [min(decode.pos) max(decode.pos)] + [-1 1].*median(diff(decode.pos))*1;
    box on
    colormap(cmap);
    rescale_plot_data('x',[1e-6 t0]);
    hax.XLim = [-1 1].*win_s;
    hax.XTick = xtick;
    hax.XTickLabelRotation = 0;
%     hax.YTick = [];
    hax.TickDir = 'out';
    hax.TickLength = [0.02 0.02];
    hax.XRuler.TickLabelGapOffset = 0;
    xlabel('Time (s)', 'Units','normalized', 'Position',[0.5 -0.15]);
    if ismember(ii_ex,[1 5])
        ylabel('Replay position (m)', 'Units','normalized', 'Position',[-0.23 .5]);
    end
    
    %% workaround to fix the image occluding the axes
    xlim(xlimits)
    plot(hax.XLim([1 1]),hax.YLim,'k-')
    plot(hax.XLim([2 2]),hax.YLim,'k-')

    %% link x axes
    linkaxes(panels_ex{ii_ex}(:),'x');
    xlim(xlimits)

end

%% add panel letters
font_size = 11;
axes(panels_ex{1}(end))
text(-0.13,1.6, 'A', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels_ex{2}(end))
text(-0.13,1.6, 'B', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels_ex{3}(end))
text(-0.13,1.6, 'C', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels_ex{4}(end))
text(-0.13,1.6, 'D', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels_ex{5}(end))
text(-0.13,1.6, 'E', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels_ex{6}(end))
text(-0.13,1.6, 'F', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels_ex{7}(end))
text(-0.13,1.6, 'G', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels_ex{8}(end))
text(-0.13,1.6, 'H', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = fig_name_str;
fig_name = [fig_name sprintf('_ex_%.2d_%d_%d_%d',replay_ex_opt,replay_examples_options(replay_ex_opt,:))];
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%
