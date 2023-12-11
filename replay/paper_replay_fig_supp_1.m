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
];
replay_ex_opt = 1; 
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

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Extended_Data_Fig_2';
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
offsets_y = linspace(2,14,2)+0.2;
offsets_y = flip(offsets_y);
clear panels
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
        y=y+H;       H=H1; panels{jj,ii}(1) = axes('position', [offset+[0 y] W H]);
        y=y+H+ymarg; H=H1; panels{jj,ii}(2) = axes('position', [offset+[0 y] W H]);
        y=y+H+ymarg; H=H2; panels{jj,ii}(3) = axes('position', [offset+[0 y] W H]);
        y=y+H+ymarg; H=H2; panels{jj,ii}(4) = axes('position', [offset+[0 y] W H]);
        y=y+H+ymarg; H=H2; panels{jj,ii}(5) = axes('position', [offset+[0 y] W H]);

%         y=y+0; panels{jj,ii}(1) = axes('position', [offset+[0 2*H1+2*H2] W 1]);
%         y=y+0; panels{jj,ii}(2) = axes('position', [offset+[0 2*H1+H2] W 1]);
%         y=y+0; panels{jj,ii}(3) = axes('position', [offset+[0 2*H1] W 1]);
%         y=y+0; panels{jj,ii}(4) = axes('position', [offset+[0  H1] W H1]);
%         y=y+0; panels{jj,ii}(5) = axes('position', [offset+[0  0] W H1]);
    end
end
panels = {panels{:}}';

%%
% W = 5;
% ii=1
% x = 2;
% panels{ii}(1) = axes('position', [x 23.8 W 1.5]);
% panels{ii}(2) = axes('position', [x 22.1 W 1.5]);
% panels{ii}(3) = axes('position', [x 20.5 W 1.5]);
% panels{ii}(4) = axes('position', [x 13.0 W 6]);
% panels{ii}(5) = axes('position', [x  5.5 W 6]);
% ii=2
% x = 8.5;
% panels{ii}(1) = axes('position', [x 23.8 W 1.5]);
% panels{ii}(2) = axes('position', [x 22.1 W 1.5]);
% panels{ii}(3) = axes('position', [x 20.5 W 1.5]);
% panels{ii}(4) = axes('position', [x 13.0 W 6]);
% panels{ii}(5) = axes('position', [x  5.5 W 6]);
% ii=3
% x = 15;
% panels{ii}(1) = axes('position', [x 23.8 W 1.5]);
% panels{ii}(2) = axes('position', [x 22.1 W 1.5]);
% panels{ii}(3) = axes('position', [x 20.5 W 1.5]);
% panels{ii}(4) = axes('position', [x 13.0 W 6]);
% panels{ii}(5) = axes('position', [x  5.5 W 6]);
% 
% for ii=1:length(panels)
%     for jj=1:length(panels{ii})
%         panels{ii}(jj).Position(2) = panels{ii}(jj).Position(2) - 2;
%     end
% end

%% load data (temp)
% ii_ex=1;
%     ex = replay_examples(ii_ex);
%     addFieldsToWorkspace(ex);
%     decode = decoding_load_data(exp_ID, epoch_type, params_opt );
%     exp = exp_load_data(exp_ID,'details','path','MUA','ripples');
%     events = decoding_load_events_quantification(exp_ID,epoch_type,params_opt,"posterior");
%     event = events([events.num] ==event_num);
%     seq = event.seq_model;
%     seq_ti = [event.start_ts event.end_ts];
%     t0 = mean(seq_ti);
%     ti = t0 + [-1 1].*win_s*1e6;
%     
%     TT = exp.ripples.stats.best_TT;
%     [LFP.signal, LFP.ts, LFP.fs, LFP.params] = LFP_load(exp_ID,TT,'band','ripple','limits_ts',ti);
%     LFP.avg_signal = nanmean(LFP.signal,[2 3]);


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
    [LFP.signal, LFP.ts, LFP.fs, LFP.params] = LFP_load(exp_ID,TT,'band','ripple','limits_ts',ti);
    LFP.avg_signal = nanmean(LFP.signal,[2 3]);
    
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
        panels_legend = axes('position', [13.5 26.6 6 0.8]);
        cla
        hold on
        lw=2;
        t = [0 0.05];
        x = [.62 .62 .62 0 0 0];
        y = [1 .5 0 1 .5 0];
        for ii=1:6
            plot(x(ii)+t,y([ii ii]),'-','LineWidth',lw,'Clipping','off');
            text(x(ii)+t(2)+0.02, y(ii), states(ii),'FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
        end
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
        hcb.Position = [18.8 25.9 1 0.25];
        hcb.Label.Rotation = 0;
%         hcb.Label.Position(1) = 1.5;
        hcb.Label.String = 'Probability';
        hcb.Ticks = [];
        text(0.85, y, '0','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
        text(1.07, y, 'Max','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','middle');
%         text(cb_offset_x_middle,1,'Max','Units','normalized','FontSize',7,'HorizontalAlignment','center');
%         text(cb_offset_x_middle,0,'0','Units','normalized','FontSize',7,'HorizontalAlignment','center');
    end

    %% plot MUA
    axes(panels{ii_ex}(5));
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
        ylabel({'Multiunit';'firing-rate';'(Hz)'}, 'Units','normalized', 'Position',[-0.12 .42]);
    end
    text(0.5, 1.15, titles_str{ii_ex}, ...
        'Units','normalized','FontWeight','normal','FontSize',8,'HorizontalAlignment','center');
    text(0.5, 1.45, sprintf('%d_%s_%s_%d',ex_num,epoch_type,exp_ID,event_num),'Interpreter','none', ...
        'Units','normalized','FontWeight','normal','FontSize',6,'HorizontalAlignment','center');
    
    %% plot LFP
    axes(panels{ii_ex}(4));
    cla
    hold on
    plot(LFP.ts, LFP.avg_signal,'k');
    xlim(ti)
    xticks([])
    yticks([])
    rescale_plot_data('x',[1e-6 t0]);
    axis off
    if ismember(ii_ex,[1 5])
        text(-0.06, .5, 'Ripple', 'Units','normalized','FontSize',8.25,'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    
    %% plot posterior (state)
    axes(panels{ii_ex}(3));
    cla
    hold on
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
    axes(panels{ii_ex}(2));
    cla
    hold on
    IX = get_data_in_ti(decode.time, ti);
    prob_t = decode.time(IX);
    prob_pos = decode.posterior_pos(:,IX);
    imagesc(prob_t, decode.pos, prob_pos);
    plot([seq.start_ts seq.end_ts],[seq.start_pos seq.end_pos],'-r','LineWidth',1);
    hax = gca;
    hax.CLim = quantile(prob_pos(:),[0.01 0.99]);
    hax.XLim = prob_t([1 end]);
    hax.YLim = [min(decode.pos) max(decode.pos)] + [-1 1].*median(diff(decode.pos))*1;
    box on
    colormap(cmap);
    alpha = 0.05;
    % area([hax.XLim(1) event.start_ts],hax.YLim([2 2]), 'FaceColor','r','FaceAlpha',alpha);
    % area([event.end_ts hax.XLim(2)],hax.YLim([2 2]),   'FaceColor','r','FaceAlpha',alpha);
    % area([event.start_ts event.end_ts],hax.YLim([2 2]),   'FaceColor','g','FaceAlpha',alpha);
    rescale_plot_data('x',[1e-6 t0]);
    hax.XLim = [-1 1].*win_s;
    hax.XTick = xtick;
    hax.XTickLabel=[];
    hax.YTick=[];
    hax.TickDir = 'out';
    hax.TickLength = [0.02 0.02];
    hax.XRuler.TickLabelGapOffset = -4;
    if ismember(ii_ex,[1 5])
        ylabel('Replay position (m)', 'Units','normalized', 'Position',[-0.07 .5]);
    end
    
    %% plot posterior (position)
    axes(panels{ii_ex}(1));
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
    alpha = 0.05;
    % area([hax.XLim(1) event.start_ts],hax.YLim([2 2]), 'FaceColor','r','FaceAlpha',alpha);
    % area([event.end_ts hax.XLim(2)],hax.YLim([2 2]),   'FaceColor','r','FaceAlpha',alpha);
    % area([event.start_ts event.end_ts],hax.YLim([2 2]),   'FaceColor','g','FaceAlpha',alpha);
    rescale_plot_data('x',[1e-6 t0]);
    hax.XLim = [-1 1].*win_s;
    hax.XTick = xtick;
    hax.XTickLabelRotation = 0;
    hax.YTick = [];
    hax.TickDir = 'out';
    hax.TickLength = [0.02 0.02];
    hax.XRuler.TickLabelGapOffset = 0;
    xlabel('Time (s)', 'Units','normalized', 'Position',[0.5 -0.15]);
    if ismember(ii_ex,[1 5])
        ylabel('Replay position (m)', 'Units','normalized', 'Position',[-0.07 .5]);
    end
    
    %% link x axes
    linkaxes(panels{ii_ex}(:),'x');
    xlim(xlimits)

end

%% add colorbar
% for ii=[4 5]
%     panel = panels{4}(ii);
%     axes(panel)
%     hcb = colorbar('southoutside');
%     hcb.Units = 'centimeters';
%     cb_offset_x = 1.05;
%     hcb.Position(1) = panel.Position(1) + panel.Position(3)*cb_offset_x;
%     cb_offset_x_middle = (hcb.Position(1)+hcb.Position(3)/2-panel.Position(1))/panel.Position(3);
%     hcb.Label.Rotation = -90;
%     hcb.Label.Position(1) = 1.5;
%     hcb.Label.String = 'Probability';
%     hcb.Ticks = [];
%     text(cb_offset_x_middle,1,'Max','Units','normalized','FontSize',7,'HorizontalAlignment','center');
%     text(cb_offset_x_middle,0,'0','Units','normalized','FontSize',7,'HorizontalAlignment','center');
% %     text(cb_offset_middle,1,clim_prctiles(2)+"%",'Units','normalized','FontSize',7,'HorizontalAlignment','center');
% %     text(cb_offset_middle,0,clim_prctiles(1)+"%",'Units','normalized','FontSize',7,'HorizontalAlignment','center');
% end

%% add panel letters
% font_size = 11;
% axes(panels{1}(1))
% text(-0.13,1.25, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
% axes(panels{1}(2))
% text(-0.13,0.9, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);
% axes(panels{1}(3))
% text(-0.13,1.1, 'c', 'Units','normalized','FontWeight','bold','FontSize',font_size);
% axes(panels{1}(4))
% text(-0.13,1.0, 'd', 'Units','normalized','FontWeight','bold','FontSize',font_size);
% axes(panels{1}(5))
% text(-0.13,1.0, 'e', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = fig_name_str;
fig_name = [fig_name sprintf('_ex_%.2d_%d_%d_%d',replay_ex_opt,replay_examples_options(replay_ex_opt,:))];
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%
