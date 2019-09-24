%% Large Scale - fig 2 - Behavioral and neural recordings from bats fliying over large spatial scales.

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_3';
fig_caption_str = 'Place cells are distributed uniformly across the tunnel, with slight over-representation at particular landmarks';
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
% =========================================================================
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
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none');
pause(0.2); % workaround to solve matlab automatically changing the axes positions...

% create panels
panel_A_size = [7 3];
panel_B_size = [3 3];
panel_C_size = [7 3];
panel_D_size = [3 3];
panel_A(1) = axes('position', [ 2 22 panel_A_size]);
panel_A(2) = axes('position', [ 2 18 panel_A_size]);
panel_B(1) = axes('position', [11 22 panel_B_size]);
panel_B(2) = axes('position', [11 18 panel_B_size]);
panel_C(1) = axes('position', [ 2 13.5 panel_C_size]);
panel_C(2) = axes('position', [ 2  9.5 panel_C_size]);
panel_D(1) = axes('position', [11 13.5 panel_D_size]);
panel_D(2) = axes('position', [11  9.5 panel_D_size]);


%% load population data
% =========================================================================
prm = PARAMS_GetAll();
cells_t = DS_get_cells_summary();
bats = [79,148,34,9861,2289];
cells_t(~ismember(cells_t.bat, bats ),:) = [];
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.details];
cells(~contains({cells.brain_area}, 'CA1')) = [];
cells(~ismember([cells.ClusterQuality], [2])) = [];
cells = cellfun(@(c)(cell_load_data(c,'details','stats')), {cells.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];
cells_details = [cells.details];
cells_ID = {cells_details.cell_ID};
stats = [cells.stats];
stats = [stats.all];
cells_ID([stats.meanFR_all]>prm.inclusion.interneuron_FR_thr)=[];
clear cells stats cells_details cells_t
cells = cellfun(@(c)(cell_load_data(c,'details','stats','meanFR','stats','inclusion','signif','fields','FR_map')), cells_ID, 'UniformOutput',0);
cells = [cells{:}];
whos cells

% load LM data
exp_ID = 'b2289_d180615';
exp = exp_load_data(exp_ID,'LM');
LM = exp.LM;
LM([1 2 end]) = [];

%% panel A - field pos CDF
% =========================================================================
% figure
for ii_dir = 1:2
%     subplot(1,2,ii_dir)
    axes(panel_A(ii_dir));
    cla
    hold on
    
    % arrange data
    signif = cat(1,cells.signif);
    cells_dir = cells([signif(:,ii_dir).TF]);
    fields = cellfun(@(x)(x{ii_dir}), {cells_dir.fields},'UniformOutput',0);
    fields =[fields{:}];
    fields( [fields.in_low_speed_area] ) = []; % TODO: decide if we want those fields near the landing balls
    
    % plot cdf
    h = cdfplot([fields.loc]);
    h.Color = prm.graphics.colors.flight_directions{ii_dir};
    h.LineWidth = 0.5;
    title('')
    ha=gca;
    ha.GridLineStyle = 'none';
    
    % plot LM
    ypos = 0.02;
    angle = 45;
    ylimits = get(gca,'ylim');
    for ii_LM=1:length(LM)
        x = LM(ii_LM).pos_proj;
        name = LM(ii_LM).name;
        h=xline(x, '-', 'color', 0.5.*[1 1 1], 'LineWidth',1);
%         text(x, ylimits(2)+ypos*diff(ylimits), LM(ii_LM).name, 'Rotation', 45, 'FontSize',8);
    end
    
    % labels & graphics
    xlabel('Position (m)', 'Units','normalized','Position',[0.5 -0.17]);
    ylabel('CDF', 'Units','normalized','Position',[-0.05 0.5]);
    ha= gca;
    ha.TickDir='out';
    ha.TickLength = [0.015 0.015];
    ha.YTick = ha.YLim;
    ha.XRuler.TickLabelGapMultiplier = -0.3;
    ha.YRuler.TickLabelGapMultiplier = 0.001;
end

axes(panel_A(1));
text(-0.13,1.1, 'A', 'Units','normalized','FontWeight','bold');

%% add arrows indicating over-represented LM
if 0
axes(panel_A(1));
xa = [70 80] - 1;
ya = [0.45 0.4] + 0.01;
[xaf,yaf] = ds2nfu(xa,ya);
h=annotation('arrow', xaf,yaf);
h.HeadWidth=5;
h.HeadLength = 5;
h.LineStyle='none';

axes(panel_A(2));
xa = [70 80] - 1;
ya = [0.45 0.4] + 0.04;
[xaf,yaf] = ds2nfu(xa,ya);
h=annotation('arrow', xaf,yaf);
h.HeadWidth=5;
h.HeadLength = 5;
h.LineStyle='none';
end



%% panel B - field vs. LM position
% =========================================================================
% figure
for ii_dir = 1:2
%     subplot(1,2,ii_dir)
    axes(panel_B(ii_dir));
    cla
    hold on
    
    %% arrange data - same as in panel C
    signif = cat(1,cells.signif);
    cells_dir = cells([signif(:,ii_dir).TF]);
    fields = cellfun(@(x)(x{ii_dir}), {cells_dir.fields},'UniformOutput',0);
    fields =[fields{:}];
    switch ii_dir
        case 1
            start_ind = 1;
            end_ind = 2;
            interp_next = 'next';
            interp_prev = 'previous';
        case 2
            start_ind = 2;
            end_ind = 1;
            interp_next = 'previous';
            interp_prev = 'next';
    end
    edges = cat(1,fields.edges_prc);
    [fields.start] = disperse(edges(:,start_ind));
    [fields.end]   = disperse(edges(:,end_ind  ));
  
    [fields.LM_nearest_by_peak]  = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.loc]  , 'nearest'));
    [fields.LM_nearest_by_start] = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.start] , 'nearest'));
    [fields.LM_nearest_by_end]   = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.end]   , 'nearest'));
    [fields.LM_next_by_peak]  = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.loc]  , interp_next));
    [fields.LM_next_by_start] = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.start] , interp_next));
    [fields.LM_next_by_end]   = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.end]   , interp_next));
    [fields.LM_prev_by_peak]  = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.loc]  , interp_prev));
    [fields.LM_prev_by_start] = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.start] , interp_prev));
    [fields.LM_prev_by_end]   = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.end]   , interp_prev));
    
    % remove fields near balls
    fields([fields.in_low_speed_area])=[];
    fields(isnan([fields.LM_nearest_by_peak])) = [];
    
    %% shuffle!
    rng(0);
    shuffle = struct();
    shuffle.n = 10000;
    limits = [min([LM.pos_proj]) max([LM.pos_proj])];
%     limits = prm.fields.valid_speed_pos;
    M = [[fields.loc];[fields.start];[fields.end]];
    shifts = rand(size(M,2), shuffle.n);
    shifts = shifts .* diff(limits);
    shifts = shiftdim(shifts,-1);
    M = M + shifts;
    M = reshape(M,size(M,1), []); % pool shuffles x fields together
    % circularly correct according to the peak
    peaks_corrected = limits(1) + mod(M(1,:), diff(limits));
    shift_correction = peaks_corrected - M(1,:);
    M = M + shift_correction;
    % remove shuffled fields same as we removed real fields 
    invalid = M(1,:)<limits(1) | M(1,:)>limits(2);
    M(:,invalid) = [];
    shuffle.peaks = M(1,:);
    shuffle.start = M(2,:);
    shuffle.end   = M(3,:);
    % get relevant landmarks
    shuffle.LM_nearest_by_peak  = interp1([LM.pos_proj], [LM.pos_proj], shuffle.peaks, 'nearest');
    shuffle.LM_nearest_by_start = interp1([LM.pos_proj], [LM.pos_proj], shuffle.start, 'nearest');
    shuffle.LM_nearest_by_end   = interp1([LM.pos_proj], [LM.pos_proj], shuffle.end,   'nearest');
    shuffle.LM_next_by_peak  = interp1([LM.pos_proj], [LM.pos_proj], shuffle.peaks, interp_next);
    shuffle.LM_next_by_start = interp1([LM.pos_proj], [LM.pos_proj], shuffle.start, interp_next);
    shuffle.LM_next_by_end   = interp1([LM.pos_proj], [LM.pos_proj], shuffle.end,   interp_next);
    shuffle.LM_prev_by_peak  = interp1([LM.pos_proj], [LM.pos_proj], shuffle.peaks, interp_prev);
    shuffle.LM_prev_by_start = interp1([LM.pos_proj], [LM.pos_proj], shuffle.start, interp_prev);
    shuffle.LM_prev_by_end   = interp1([LM.pos_proj], [LM.pos_proj], shuffle.end,   interp_prev);
        
    %% plot 
    c = prm.graphics.colors.flight_directions{ii_dir};
    x1 = [fields.loc] - [fields.LM_nearest_by_peak];
    h1 = histogram( x1 );
    h1.Normalization = 'pdf';
%     h1.NumBins = h1.NumBins * 2;
    h1.FaceColor = c;
    x2 = shuffle.peaks - shuffle.LM_nearest_by_peak;
    h2 = histogram( x2 );
%     h2.BinEdges = h1.BinEdges;
    h2.Normalization = 'pdf';
    h2.DisplayStyle = 'stairs';
    h2.EdgeColor = 'k';
    h2.LineWidth = 2;
    [H,P,KSSTAT] = kstest2(x1,x2);
    text(1,1,sprintf('P_{KS}=%.2f',P),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');
    xline(0,'-','LineWidth',2);
    
    % labels & graphics
    ha= gca;
    ha.TickDir='out';
    ha.TickLength = [0.03 0.03];
    ha.YTick = ha.YLim;
    ha.XRuler.TickLabelGapMultiplier = -0.3;
    ha.YRuler.TickLabelGapMultiplier = 0.001;
    xlabel('Distance from landmark (m)', 'Units','normalized','Position',[0.5 -0.13])
    ylabel('PDF', 'Units','normalized','Position',[-0.08 0.5])
    
    % labels & graphics
    ha= gca;
    ha.TickDir='out';
    ha.TickLength = [0.03 0.03];
    ha.YTick = ha.YLim;
    ha.XRuler.TickLabelGapMultiplier = -0.3;
    ha.YRuler.TickLabelGapMultiplier = 0.001;
    xlabel({'Distance to nearest landmark (m)'}, 'Units','normalized','Position',[0.5 -0.13])
    ylabel('PDF', 'Units','normalized','Position',[-0.08 0.5])
end

axes(panel_B(1));
text(-0.3,1.1, 'B', 'Units','normalized','FontWeight','bold');


%% panel C - field size vs. pos
% =========================================================================
% figure
for ii_dir = 1:2
%     subplot(1,2,ii_dir)
    axes(panel_C(ii_dir));
    cla
    hold on
    
    % arrange data
    signif = cat(1,cells.signif);
    cells_dir = cells([signif(:,ii_dir).TF]);
    fields = cellfun(@(x)(x{ii_dir}), {cells_dir.fields},'UniformOutput',0);
    fields =[fields{:}];
    fields( [fields.in_low_speed_area] ) = [];
    
    % plot cdf
    c = prm.graphics.colors.flight_directions{ii_dir};
    x = [fields.loc];
    y = [fields.width_prc];
    y_clipping = 20;
    y(y>y_clipping) = y_clipping;
%     plot(x,y,'.','MarkerSize',4, 'Color',c);
    plot(x(y< y_clipping), y(y< y_clipping),'.','MarkerSize',4, 'Color',c);
    plot(x(y>=y_clipping), y(y>=y_clipping),'o','MarkerSize',4, 'Color',c);
%     yline(20)
    
    % plot LM
    ypos = 0.02;
    angle = 45;
    ylimits = get(gca,'ylim');
    for ii_LM=1:length(LM)
        x = LM(ii_LM).pos_proj;
        name = LM(ii_LM).name;
        h=xline(x, '-', 'color', 0.9.*[1 1 1], 'LineWidth',0.5);
%         text(x, ylimits(2)+ypos*diff(ylimits), LM(ii_LM).name, 'Rotation', 45, 'FontSize',8);
    end
    
    % labels & graphics
    xlabel('position (m)', 'Units','normalized','Position',[0.5 -0.17]);
    ylabel('Field size (m)', 'Units','normalized','Position',[-0.05 0.5]);
    ha= gca;
    ha.TickDir='out';
    ha.TickLength = [0.015 0.015];
    ha.XTick = [0:50:200];
    ha.YLim = [0 y_clipping];
    ha.YTick = ha.YLim;
    ha.XRuler.TickLabelGapMultiplier = -0.3;
    ha.YRuler.TickLabelGapMultiplier = 0.001;
end

axes(panel_C(1));
text(-0.13,1.1, 'C', 'Units','normalized','FontWeight','bold');


%% panel D - field size - near vs. far from LM
% =========================================================================
% figure
axes(panel_D(1));
text(-0.3,1.1, 'D', 'Units','normalized','FontWeight','bold');
for ii_dir = 1:2
%     subplot(1,2,ii_dir)
    axes(panel_D(ii_dir));
    cla
    hold on
    
    %% arrange data
    signif = cat(1,cells.signif);
    cells_dir = cells([signif(:,ii_dir).TF]);
    fields = cellfun(@(x)(x{ii_dir}), {cells_dir.fields},'UniformOutput',0);
    fields =[fields{:}];
    switch ii_dir
        case 1
            start_ind = 1;
            end_ind = 2;
            interp_next = 'next';
            interp_prev = 'previous';
        case 2
            start_ind = 2;
            end_ind = 1;
            interp_next = 'previous';
            interp_prev = 'next';
    end
    edges = cat(1,fields.edges_prc);
    [fields.start] = disperse(edges(:,start_ind));
    [fields.end]   = disperse(edges(:,end_ind  ));
  
    [fields.LM_nearest_by_peak]  = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.loc]  , 'nearest'));
    [fields.LM_nearest_by_start] = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.start] , 'nearest'));
    [fields.LM_nearest_by_end]   = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.end]   , 'nearest'));
    [fields.LM_next_by_peak]  = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.loc]  , interp_next));
    [fields.LM_next_by_start] = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.start] , interp_next));
    [fields.LM_next_by_end]   = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.end]   , interp_next));
    [fields.LM_prev_by_peak]  = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.loc]  , interp_prev));
    [fields.LM_prev_by_start] = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.start] , interp_prev));
    [fields.LM_prev_by_end]   = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.end]   , interp_prev));
    
    % remove fields near balls
    fields([fields.in_low_speed_area])=[];
    fields(isnan([fields.LM_nearest_by_peak])) = [];
    
    %% plot 
    c = prm.graphics.colors.flight_directions{ii_dir};
    x = abs([fields.loc] - [fields.LM_nearest_by_peak]);
    y = [fields.width_prc];
%     thr = median(x);
    thr = 5;
    y1 = y(x<thr);
    y2 = y(x>=thr);
    
    fprintf('y1: n=%d/%d (%2.0f%%)\r',length(y1),length(y),100*length(y1)/length(y));
    fprintf('y2: n=%d/%d (%2.0f%%)\r',length(y2),length(y),100*length(y2)/length(y));
    
    nBins = 25;
    h1 = histogram( y1 );
    h1.NumBins = nBins;
    h1.Normalization = 'pdf';
    h1.DisplayStyle = 'stairs';
    h1.EdgeColor = c;
    h1.LineWidth = 1;
    h2 = histogram( y2 );
    h2.NumBins = nBins;
    h2.Normalization = 'pdf';
    h2.DisplayStyle = 'stairs';
    h2.EdgeColor = c;
    h2.LineWidth = 2;
    [H,P,KSSTAT] = kstest2(y1,y2);
    text(1,1,sprintf('P_{KS}=%.2f',P),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');
    
    % labels & graphics
    ha= gca;
    ha.TickDir='out';
    ha.TickLength = [0.03 0.03];
    ha.YTick = ha.YLim;
    ha.XLim = [0 20];
    ha.XRuler.TickLabelGapMultiplier = -0.3;
    ha.YRuler.TickLabelGapMultiplier = 0.001;
    xlabel('Field size (m)', 'Units','normalized','Position',[0.5 -0.13]);
    ylabel('PDF', 'Units','normalized','Position',[-0.08 0.5]);
%     legend_pos = [ha.Position([1 2])+[2.2 0.5].*ha.Position([3 4]) 0.1 0.03];
%     legend({'<thr';'>thr'},'Location','eastoutside','Units','centimeters','Position',legend_pos);
%     text(2.2,0.15, {"thr="+thr;...
%                     sprintf('y1: n=%d/%d (%2.0f%%)',length(y1),length(y),100*length(y1)/length(y));...
%                     sprintf('y2: n=%d/%d (%2.0f%%)',length(y2),length(y),100*length(y2)/length(y))
%         }, 'Units','normalized','HorizontalAlignment','center');
    
    % add legend
    leg_ax = axes('position', [ha.Position([1 2])+[0.65 0.5].*ha.Position([3 4]) [0.2 0.3].*ha.Position([3 4])]);
    hold on
    plot([1 2],[1 1], 'Color',c, 'LineWidth',1);
    plot([1 2],[2 2], 'Color',c, 'LineWidth',2);
    xlim([0 3]);
    ylim([0 3]);
    text(3,1, "Distance$<$"   +thr+"m",'FontSize',7,'Interpreter','latex')
    text(3,2, "Distance$\geq$"+thr+"m",'FontSize',7,'Interpreter','latex')
    set(gca,'Visible','off');
end


%% add direction arrows
arrow_x = 0.1 +[0 0.05];
arrow_y = repelem(0.96,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.008,  'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y      ,  'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);






%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');







%%

