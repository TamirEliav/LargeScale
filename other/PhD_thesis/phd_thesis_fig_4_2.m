%% PhD thesis figure 5.2 - replay examples

%%
close all
clear 
clc

%% plotting options

%% define output files
res_dir = 'E:\Tamir\PhD\Thesis\resources\ch_4_seq';
mkdir(res_dir)
fig_name_str = 'Fig_5_2_replay_examples';
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
figure ;
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
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

% create panels
panel_A_sizes = [
    5 4;
    5 .8;
    5 .8;
    5 .8;
    ];
panel_A_vmargin = 0.1;
panel_A_pos = [2 2];
panel_A = [];
for ii=1:3
    for jj=1:3
        offset_x = (ii-1)*6;
        offset_y = (jj-1)*7.5;
        offset = panel_A_pos + [offset_x offset_y];
        for ii_sub = 1:size(panel_A_sizes,1)
            panel_A(ii,jj,ii_sub) = axes('position', [offset+[0 sum(panel_A_sizes(1:(ii_sub-1),2))+(ii_sub-1)*panel_A_vmargin] panel_A_sizes(ii_sub,:) ]);
        end
    end
end
panel_A = panel_A(:,3:-1:1,:);
panel_A = reshape(panel_A,[9 4]);
% panel_A_legend = axes('position', [2.2 24.7 0.5 0.15]);

%% examples list
examples_list = [
    struct('exp_ID','b0184_d191130', 'epoch_type','sleep', 'event_num',139),
    struct('exp_ID','b0184_d191130', 'epoch_type','sleep', 'event_num',143),
    struct('exp_ID','b0184_d191130', 'epoch_type','sleep', 'event_num',279),
    struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',10),
    struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',124),
    struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',204),
    struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',206),
    struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',209),
    struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',26),
    struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',154),
    struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',159),
    struct('exp_ID','b0184_d191202', 'epoch_type','sleep', 'event_num',234),
    struct('exp_ID','b0184_d191205', 'epoch_type','sleep', 'event_num',98),
    struct('exp_ID','b0184_d191208', 'epoch_type','sleep', 'event_num',51),
    struct('exp_ID','b0184_d191208', 'epoch_type','sleep', 'event_num',64),
    struct('exp_ID','b0184_d191208', 'epoch_type','sleep', 'event_num',70),
    struct('exp_ID','b0184_d191209', 'epoch_type','sleep', 'event_num',78),
    struct('exp_ID','b9861_d180526', 'epoch_type','sleep', 'event_num',49),
    ];
% examples_IX = [1:9];
% examples_IX = [10:18];
% examples_IX = [1 3 4 6 14 8 9 13 18];
% examples_IX = [1 4 6 7 8 11 12 14 18];
examples_IX = [1 4 6 13 8 11 12 14 18]; % chosen examples
examples_list = examples_list(examples_IX);

%% graphics
event_margins = 0.2; % in %

%% load data for panels A and B
for ii_example = 1:size(panel_A,1)
    %% load data
    exp_ID = examples_list(ii_example).exp_ID;
    epoch_type = examples_list(ii_example).epoch_type;
    event_num = examples_list(ii_example).event_num;
    decode = decoding_load_data(exp_ID, 'sleep', 11);
    exp = exp_load_data(exp_ID,'details','path','MUA','ripples');
    TT = exp.ripples.stats.best_TT;
    [LFP.signal, LFP.ts, LFP.fs, LFP.params] = LFP_load(exp_ID,TT,'band','ripple');
    LFP.avg_signal = nanmean(LFP.signal,[2 3]);
    events = decoding_load_events_quantification(exp_ID,epoch_type,11,"posterior");
    event = events(event_num);
    ti = [event.start_ts event.end_ts];
    ti = ti + [-1 1].*event_margins*range(ti);

    %% plot LFP
    axes(panel_A(ii_example, 4));
    cla
    hold on
    plot(LFP.ts, LFP.avg_signal,'k');
    xlim(ti)
    xticks([])
    yticks([])
    rescale_plot_data('x',[1e-6 ti(1)]);
    axis off
    text(0.5,1.1,sprintf('%s_event_%d',exp_ID,event_num),'HorizontalAlignment','center','FontSize',7,'Units','normalized','Interpreter','none');

    %% plot MUA
    axes(panel_A(ii_example, 3));
    cla
    hold on
%     plot(exp.MUA.t, exp.MUA.zFR, 'k');
    area(exp.MUA.t, exp.MUA.FR,'FaceColor','k');
    xlim(ti)
    xticks([])
%     yticks(round(max()))
    rescale_plot_data('x',[1e-6 ti(1)]);

    %% plot posterior (states)
    axes(panel_A(ii_example, 2));
    cla
    hold on
    plot(decode.time, decode.posterior_state,'LineWidth',2);
    xlim(ti)
    ylim([0 1])
    xticks([])
    yticks([0 1])
    rescale_plot_data('x',[1e-6 ti(1)]);

    %% plot posterior (position)
    axes(panel_A(ii_example, 1));
    cla
    hold on
    IX = get_data_in_ti(decode.time, ti);
    prob_t = decode.time(IX);
    prob_pos = decode.posterior_pos(:,IX);
    imagesc(prob_t, decode.pos, prob_pos);
    hax = gca;
    hax.TickDir = 'out';
    hax.TickLength(1) = 0.008;
    hax.CLim = quantile(prob_pos(:),[0.01 0.99]);
    hax.XLim = prob_t([1 end]);
    hax.YLim = [min(decode.pos) max(decode.pos)] + [-1 1].*median(diff(decode.pos))*1;
    box on
    hax.XRuler.TickLabelGapOffset = -4;
    cmap = bone;
    cmap = flipud(cmap);
    colormap(cmap);
    rescale_plot_data('x',[1e-6 ti(1)]);
    hax.XLim = [0 range(ti)*1e-6];

    %% link x axes
    linkaxes(panel_A(ii_example, :),'x')
end

%% add labels
axes(panel_A(1, 1));
ylabel('Position (m)')
axes(panel_A(4, 1));
ylabel('Position (m)')
axes(panel_A(7, 1));
ylabel('Position (m)')
xlabel('Time (s)')
axes(panel_A(8, 1));
xlabel('Time (s)')
axes(panel_A(9, 1));
xlabel('Time (s)')
axes(panel_A(1, 2));
ylabel('Prob.')
axes(panel_A(1, 3));
ylabel('MUA (Hz)')
axes(panel_A(4, 2));
ylabel('Prob.')
axes(panel_A(4, 3));
ylabel('MUA (Hz)')
axes(panel_A(7, 2));
ylabel('Prob.')
axes(panel_A(7, 3));
ylabel('MUA (Hz)')

%% save fig
fig_name_out = fullfile(res_dir, [fig_name_str '_' char(strjoin(string(examples_IX),'_'))]);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
disp('figure was successfully saved to pdf/tiff/fig formats');


