%% Replay - Fig S1 - many replay examples
%%
clear 
clc
close all

%% replay examples - new format
replay_examples_list_filename = "L:\Analysis\Code\inclusion_lists\replay_examples.xlsx";
replay_examples_list = table2struct(readtable(replay_examples_list_filename));
replay_examples_options = [
500 328 348 555 359 378 ...
540 394 543 460 557 450 ...
561 339 463 364 493 288 ...
168  58  76  15  39 112 ... 
177 109 113  62  51  41 ...
238 228 216 257 253  83 ...
;
500 328 348 555 359 378 ...
516 394 543 460 557 450 ...
561 339 463 364 493 288 ...
168  58  76  15  39 112 ... 
177 109 113  62  51  41 ...
193 228 131 257 253  83 ...
];
replay_ex_opt = 2; 
replay_examples = replay_examples_list(replay_examples_options(replay_ex_opt,:))

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Extended_Data_Fig_1';
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
close all
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
clear panels
nrow = 6;
ncol = 6;
offsets_x = linspace(0,15,ncol)+2;
offsets_y = linspace(3,22,nrow);
offsets_y(nrow/2+1:end) = offsets_y(nrow/2+1:end)+1;
offsets_y = flip(offsets_y);
offsets_y = offsets_y - 0.5;
nrow = length(offsets_y);
ncol = length(offsets_x);
for ii=1:nrow
    for jj=1:ncol
        offset_x = offsets_x(jj);
        offset_y = offsets_y(ii);
        offset = [offset_x offset_y];
        panels{1}(jj,ii,1) = axes('position', [offset+[0 0] 2.3 2.8]);
        panels{1}(jj,ii,2) = axes('position', [offset+[0 2.8] 2.3 .3]);
%         panels{1}(jj,ii,3) = axes('position', [offset+[0 4] 2.3 .5]);
    end
end

total_offset = [.5 -0.8];
for ii = 1:length(panels)
    subpanels = panels{ii};
    subpanels = subpanels(:);
    for jj = 1:length(subpanels)
        subpanels(jj).Position([1 2]) = subpanels(jj).Position([1 2]) + total_offset;
    end
end


%%
% ii_ex=1;
% win_s = 0.5;
% ex = replay_examples(ii_ex);
% addFieldsToWorkspace(ex);
% decode = decoding_load_data(exp_ID, epoch_type, params_opt);
% exp = exp_load_data(exp_ID,'details','path','MUA','ripples');
% events = decoding_load_events_quantification(exp_ID,epoch_type,params_opt,"posterior");
% event = events(event_num);
% seq = event.seq_model;
% seq_ti = [event.start_ts event.end_ts];
% t0 = mean(seq_ti);
% ti = t0 + [-1 1].*win_s*1e6;
% 
% TT = exp.ripples.stats.best_TT;
% [LFP.signal, LFP.ts, LFP.fs, LFP.params] = LFP_load(exp_ID,TT,'band','ripple','limits_ts',ti);
% LFP.avg_signal = nanmean(LFP.signal,[2 3]);


%% plot replay examples
win_s = 0.5;
panels_ex = panels{1}
panels_ex = reshape(panels_ex, size(panels_ex,1)*size(panels_ex,2), size(panels_ex,3))
cmap = bone;
cmap = flipud(cmap);
for ii_ex = 1:size(panels_ex,1)
    
    %% load data
    ex = replay_examples(ii_ex);
    addFieldsToWorkspace(ex);
    decode = decoding_load_data(exp_ID, epoch_type, params_opt);
    exp = exp_load_data(exp_ID,'details','path','MUA','ripples');
    events = decoding_load_events_quantification(exp_ID,epoch_type,params_opt,"posterior");
    event = events(event_num);
    seq = event.seq_model;
    seq_ti = [event.start_ts event.end_ts];
    t0 = mean(seq_ti);
    ti = t0 + [-1 1].*win_s*1e6;
    
    TT = exp.ripples.stats.best_TT;
    [LFP.signal, LFP.ts, LFP.fs, LFP.params] = LFP_load(exp_ID,TT,'band','ripple','limits_ts',ti);
    LFP.avg_signal = nanmean(LFP.signal,[2 3]);

    %% plot LFP
%     axes(panels_ex(ii_ex,3));
%     cla
%     hold on
%     plot(LFP.ts, LFP.avg_signal,'k');
%     xlim(seq_ti+[-1 1].*0.2*range(seq_ti))
%     xticks([])
%     yticks([])
%     rescale_plot_data('x',[1e-6 seq_ti(1)]);
%     axis off
%     text(0.5,1.07,sprintf('%d_%s_%s_%d',ex_num,epoch_type,exp_ID,event_num),'units','normalized','Interpreter','none','FontWeight','normal','FontSize',5,'HorizontalAlignment','center');
    
    %% plot posterior (state)
    axes(panels_ex(ii_ex,2));
    cla
    hold on
    IX = get_data_in_ti(decode.time, ti);
    prob_t = decode.time(IX);
    prob_state = squeeze(decode.posterior_state(event.state_num,IX));
    plot(prob_t, prob_state, 'k','LineWidth',2);
    hax = gca;
    hax.XLim = prob_t([1 end]);
    hax.YLim = [0 1];
    hax.YTick = [0 1];
    hax.XTickLabel = [];
    box on
    colormap(cmap);
    hax.XLim = seq_ti+[-1 1].*0.2*range(seq_ti);
    hax.TickDir = 'out';
    hax.TickLength = [0.02 0.02];
    hax.XRuler.TickLabelGapOffset = -4;
    rescale_plot_data('x',[1e-6 seq_ti(1)]);
    text(0.5,1.01,sprintf('%d_%s_%s_%d',ex_num,epoch_type,exp_ID,event_num),'units','normalized','Interpreter','none','FontWeight','normal','FontSize',5,'HorizontalAlignment','center','VerticalAlignment','bottom');
    
    %% plot posterior (position)
    axes(panels_ex(ii_ex,1));
    cla
    hold on
    IX = get_data_in_ti(decode.time, ti);
    prob_t = decode.time(IX);
    prob_pos = squeeze(decode.posterior(:,event.state_num,IX));
    imagesc(prob_t, decode.pos, prob_pos);
    plot([seq.start_ts seq.end_ts],[seq.start_pos seq.end_pos],'-r','LineWidth',0.8);
    hax = gca;
    clim_prctiles = [1 99];
    hax.CLim = prctile(prob_pos(:),clim_prctiles);
%     hax.CLim = [0 max(prob_pos(:))];
    hax.XLim = prob_t([1 end]);
    hax.YLim = [min(decode.pos) max(decode.pos)] + [-1 1].*median(diff(decode.pos))*1;
    box on
    colormap(cmap);
    hax.XLim = seq_ti+[-1 1].*0.2*range(seq_ti);
    hax.TickDir = 'out';
    hax.TickLength = [0.02 0.02];
    hax.XRuler.TickLabelGapOffset = -2;
    hax.YRuler.TickLabelGapOffset = -1;
    if ismember(ii_ex,[[1:ncol]+12 [1:ncol]+30]);
        xlabel('Time (s)','Units','normalized','Position',[0.5 -0.105]);
    end
    rescale_plot_data('x',[1e-6 seq_ti(1)]);
    
    %% link x axes
    linkaxes(panels_ex(ii_ex,:),'x');

    %% add colorbar
    if ii_ex==size(panels{1},1)
        hcb = colorbar('east');
        hcb.Units = 'centimeters';
        cb_offset = 1.2;
        hcb.Position(1) = panels_ex(ii_ex,1).Position(1) + panels_ex(ii_ex,1).Position(3)*cb_offset;
        cb_offset_middle = (hcb.Position(1)+hcb.Position(3)/2-panels_ex(ii_ex,1).Position(1))/panels_ex(ii_ex,1).Position(3);
        hcb.Label.Rotation = -90;
        hcb.Label.Position(1) = 1.5;
        hcb.Label.String = 'Probability';
        hcb.Ticks = [];
        text(cb_offset_middle,1,'Max','Units','normalized','FontSize',7,'HorizontalAlignment','center');
        text(cb_offset_middle,0,'0','Units','normalized','FontSize',7,'HorizontalAlignment','center');
%         text(cb_offset_middle,1,clim_prctiles(2)+"%",'Units','normalized','FontSize',7,'HorizontalAlignment','center');
%         text(cb_offset_middle,0,clim_prctiles(1)+"%",'Units','normalized','FontSize',7,'HorizontalAlignment','center');
%         hcb.Ticks = hcb.Limits;
%         hcb.TickLabels = {'0','Max'};
    end
end

%%
for ii_row = 1:nrow
    axes(panels{1}(1,ii_row,1))
    ylabel('Position (m)','Units','normalized','Position',[-0.3 0.5]);
    axes(panels{1}(1,ii_row,2))
    ylabel('Prob.','Units','normalized','Position',[-0.3 0.5]);
end

%% add panel letters
font_size = 11;
axes(panels{1}(1,1,end));
text(-0.5,3, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
text(-0.05,3, 'Sleep', 'Units','normalized','FontWeight','normal','FontSize',9);
axes(panels{1}(1,4,end));
text(-0.5,3, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);
text(-0.05,3, 'Awake', 'Units','normalized','FontWeight','normal','FontSize',9);

%% print/save the figure
fig_name_out = fullfile(res_dir, sprintf('%s',fig_name_str));

print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% exportgraphics(gcf,fullfile(res_dir,'testexport.pdf'),'BackgroundColor','none','ContentType','vector');

disp('figure was successfully saved to pdf/tiff/fig formats');

%%

