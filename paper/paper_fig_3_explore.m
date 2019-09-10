%%
clear
clc

%%
dir_out = 'L:\paper_figures';

%% load population data
prm = PARAMS_GetAll();
dir_colors = prm.graphics.colors.flight_directions;
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

%% load LM info
exp_ID = 'b2289_d180615';
exp = exp_load_data(exp_ID,'LM');

%% plot population fields distribution - lines and edges
% =========================================================================
figure('Units','normalized','Position',[0 0 1 1])
pnl = panel();
pnl.pack('h',2);
pnl(1).pack('v',[30 70]);
pnl(2).pack('v',[30 70]);
pnl.margin = 30;
pnl.de.margin = 15;
h=pnl.title('population fields distribution - peak, width and edges');
h.FontSize = 16;
h.Position = [0.5 1.07];
vline_w = 1;
vline_h = 1*[-1 1];
color_edges_L = 'g';
color_edges_R = 'r';
color_peaks = 'm';
for ii_dir = 1:2
    
    %% field line/edges
    pnl(ii_dir,2).select();
    hold on
    signif = cat(1,cells.signif);
    cells_dir = cells([signif(:,ii_dir).TF]);
    stats = cat(1,cells_dir.stats);
    stats = cat(1,stats.dir);
    SI = [stats(:,ii_dir).SI_bits_spike];
    [SI,sort_IX] = sort(SI,'descend');
    cells_dir = cells_dir(sort_IX);
    fields_edges = [];
    fields_peaks = [];
    for ii_cell = 1:length(cells_dir)
        cell = cells_dir(ii_cell);
        fields = cell.fields{ii_dir};
        x = cat(1,fields.edges_prc);
        fields_edges = [fields_edges ; x];
        y = ii_cell;
        plot( x',[y y],'-','Color',prm.graphics.colors.flight_directions{ii_dir},'LineWidth',0.5);
        plot( x(:,[1 1])', y+vline_h,'-','Color',color_edges_L,'LineWidth',vline_w);
        plot( x(:,[2 2])', y+vline_h,'-','Color',color_edges_R,'LineWidth',vline_w);
        x = cat(1,fields.loc);
        fields_peaks = [fields_peaks ; x];
        plot( [x x]', y+vline_h,'-','Color',color_peaks,'LineWidth',vline_w);
    end
    ylim([0 length(cells_dir)+1])
    xlabel('Position (m)')
    ylabel('Cell no.')
    
    %% plot peak/edges distribution
    pnl(ii_dir,1).select();
    hold on
    hold on
    x = 0:0.2:200;
    [f,xi] = ksdensity(fields_peaks,x,'Kernel','normal','Bandwidth',1.5);
    plot(xi,f, '-', 'LineWidth',0.01, 'Color',color_peaks);
    [f,xi] = ksdensity(fields_edges(:,1),x,'Kernel','normal','Bandwidth',1.5);
    plot(xi,f, '-', 'LineWidth',0.01, 'Color',color_edges_L);
    [f,xi] = ksdensity(fields_edges(:,2),x,'Kernel','normal','Bandwidth',1.5);
    plot(xi,f, '-', 'LineWidth',0.01, 'Color',color_edges_R);
    ha = gca;
    ha.XLim = [0 200];
    ha.XTickLabel = {};
    xlabel('Position (m)')
    ylabel('density')
    title("direction" +ii_dir,'Units','normalized','Position',[0.5 1.15])
    hl=columnlegend(3,{' peak','left ','right '});
    hl.Position = ha.Position.*[1 1 0.3 0.1] + ha.Position([3 4 1 2]).*[0.3 0.8 0 0];
    % add LM
    plot_LM(exp.LM,0.01);
end
% save figure
figname = 'fields_distribution';
file_out = fullfile(dir_out,figname);
saveas(gcf, file_out, 'tif')
% saveas(gcf, file_out, 'pdf')


%% plot FR maps ordered by field location (heatmap)
% =========================================================================
figure('Units','normalized','Position',[0 0 1 1])
pnl = panel();
pnl.pack('h',2);
pnl.margin = 30;
h=pnl.title('FR maps ordered by field location and normalized to field peak');
h.FontSize = 16;
h.Position = [0.5 1.07];
for ii_dir = 1:2
    %% arrange data
    signif = cat(1,cells.signif);
    cells_dir = cells([signif(:,ii_dir).TF]);
    pop_data = struct('field',{}, 'PSTH',{}, 'fields',{});
    ii_field = 1;
    for ii_cell = 1:length(cells_dir)
        cell = cells_dir(ii_cell);
        fields = cell.fields{ii_dir};
        for jj_field = 1:length(fields)
            pop_data(ii_field).field = fields(jj_field);
            pop_data(ii_field).PSTH = cell.FR_map(ii_dir).all.PSTH;
            pop_data(ii_field).fields = fields;
            ii_field = ii_field + 1;
        end
    end
    locs = arrayfun(@(x)(x.loc), [pop_data.field]);
    [locs, sort_IX] = sort(locs,'ascend');
    pop_data = pop_data(sort_IX);
    
    PSTH_pos_bins = cell.FR_map(ii_dir).all.bin_centers;
    M = nan(length(pop_data),length(PSTH_pos_bins));
    for ii_field = 1:length(pop_data)
        PSTH = pop_data(ii_field).PSTH;
        PSTH = PSTH ./ pop_data(ii_field).field.peak;
        % TODO: consider normalizding to the href value instead of the peak
        % value, so we can visualize the field width
        M(ii_field, :) = PSTH;
    end

    %% plot
    pnl(ii_dir).select();
    hold on
    m = 1.5;
    M(M>m) = m;
    imagesc(PSTH_pos_bins,1:size(M,1),M)
    
    colorbar
    set(gca,'CLim',[0 m]);
    colormap parula
%     colormap jet
    ylim([0 size(M,1)+1]);
    xlabel('Position (m)')
    ylabel('Field no.')
    title("direction" +ii_dir,'Units','normalized','Position',[0.5 1.05])
    % add LM
%     plot_LM(exp.LM,0.01);
    
    
end
% save figure
figname = 'FR_map_heatmap_by_fields';
file_out = fullfile(dir_out,figname);
saveas(gcf, file_out, 'tif')
% saveas(gcf, file_out, 'pdf')

%% plot fields (lines) ordered by field location
% =========================================================================
figure('Units','normalized','Position',[0 0 1 1])
pnl = panel();
pnl.pack('h',2);
pnl.margin = 30;
h=pnl.title('Fields ordered by field location');
h.FontSize = 16;
h.Position = [0.5 1.07];
for ii_dir = 1:2
    %% arrange data
    signif = cat(1,cells.signif);
    cells_dir = cells([signif(:,ii_dir).TF]);
    pop_data = struct('field',{}, 'PSTH',{}, 'fields',{});
    ii_field = 1;
    for ii_cell = 1:length(cells_dir)
        cell = cells_dir(ii_cell);
        fields = cell.fields{ii_dir};
        for jj_field = 1:length(fields)
            pop_data(ii_field).field = fields(jj_field);
            pop_data(ii_field).PSTH = cell.FR_map(ii_dir).all.PSTH;
            pop_data(ii_field).fields = fields;
            ii_field = ii_field + 1;
        end
    end
    locs = arrayfun(@(x)(x.loc), [pop_data.field]);
    [locs, sort_IX] = sort(locs,'ascend');
    pop_data = pop_data(sort_IX);
    sdf=[pop_data.field];
    pop_data([sdf.in_low_speed_area]) = [];
    
    %% plot
    c = prm.graphics.colors.flight_directions{ii_dir};
    pnl(ii_dir).select(); 
    hold on
    for ii_field = 1:length(pop_data)
        fields = pop_data(ii_field).fields;
        fields([fields.in_low_speed_area])=  [];
        x = cat(1,fields.edges_prc)';
        y = ones(size(x)) .* ii_field;
        plot(x,y,'Color',c);
    end
    ylim([0 length(pop_data)+1]);
    xlabel('Position (m)')
    ylabel('Field no.')
    title("direction" +ii_dir,'Units','normalized','Position',[0.5 1.05])
    % add LM
    plot_LM(exp.LM,0.01);
    
end
% save figure
figname = 'fields_by_fields_location_with_LM';
% figname = 'fields_by_fields_location_no_LM';
% figname = 'fields_by_fields_location';
file_out = fullfile(dir_out,figname);
saveas(gcf, file_out, 'tif')
% saveas(gcf, file_out, 'pdf')



%% Fields size/position vs. landmarks
% =========================================================================
figure('Units','normalized','Position',[0 0 1 1])
pnl = panel();
pnl.pack('h',2);
pnl(1).pack('v',2);
pnl(2).pack('v',2);
pnl(1,1).pack('h',3);
pnl(1,2).pack('h',3);
pnl(2,1).pack('h',3);
pnl(2,2).pack('h',3);
pnl.margin = 30;
h=pnl.title('Fields size/position vs. landmarks');
h.FontSize = 16;
h.Position = [0.5 1.07];
% pnl.select('all');
% pnl.identify();

LM = exp.LM;
LM( contains({LM.name},{'ball','enter'}) ) = [];
% LM( contains({LM.name},{'enter'}) ) = [];
for ii_dir = 1:2
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
    
    %% shuffle!
%     shuffle.n = 20000;
%     shuffle.pos = linspace(prm.fields.valid_speed_pos(1), prm.fields.valid_speed_pos(end), shuffle.n);
%     shuffle.LM_nearest = interp1([LM.pos_proj], [LM.pos_proj], shuffle.pos, 'nearest');
%     shuffle.LM_next = interp1([LM.pos_proj], [LM.pos_proj], shuffle.pos, interp_next);
%     shuffle.LM_prev = interp1([LM.pos_proj], [LM.pos_proj], shuffle.pos, interp_prev);
    
    %% shuffle!
    shuffle.n = 1000;
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
    
    %% plot distance to LM (histograms)
    c = prm.graphics.colors.flight_directions{ii_dir};
    
    pnl(ii_dir,1,1).select(); hold on
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
    text(1,1,sprintf('ks: p=%.3f',P),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');
    xline(0,'-','LineWidth',2);
    xlabel('Distance from landmark (m)')
    title('field peak to nearest LM')
    
    pnl(ii_dir,1,2).select(); hold on
    x1 = [fields.start] - [fields.LM_prev_by_peak];
    h1 = histogram( x1 );
    h1.Normalization = 'pdf';
%     h1.NumBins = h1.NumBins * 2;
    h1.FaceColor = c;
    x2 = shuffle.start - shuffle.LM_prev_by_peak;
    h2 = histogram( x2 );
%     h2.BinEdges = h1.BinEdges;
    h2.Normalization = 'pdf';
    h2.DisplayStyle = 'stairs';
    h2.EdgeColor = 'k';
    h2.LineWidth = 2;
    [H,P,KSSTAT] = kstest2(x1,x2);
    text(1,1,sprintf('ks: p=%.3f',P),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');
    xline(0,'-','LineWidth',2);
    xlabel('Distance from landmark (m)')
    title('field start to previous LM')
    
    pnl(ii_dir,1,3).select(); hold on
    x1 = [fields.end] - [fields.LM_next_by_peak];
    h1 = histogram( x1 );
    h1.Normalization = 'pdf';
%     h1.NumBins = h1.NumBins * 2;
    h1.FaceColor = c;
    x2 = shuffle.end - shuffle.LM_next_by_peak;
    h2 = histogram( x2 );
%     h2.BinEdges = h1.BinEdges;
    h2.Normalization = 'pdf';
    h2.DisplayStyle = 'stairs';
    h2.EdgeColor = 'k';
    h2.LineWidth = 2;
    [H,P,KSSTAT] = kstest2(x1,x2);
    text(1,1,sprintf('ks: p=%.3f',P),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');
    xline(0,'-','LineWidth',2);
    xlabel('Distance from landmark (m)')
    title('field end to next LM')
    
    %% plot size vs. distance to LM (scatter)
    c = prm.graphics.colors.flight_directions{ii_dir};
    
    pnl(ii_dir,2,1).select(); hold on
    x = [fields.loc] - [fields.LM_nearest_by_peak];
    y = [fields.width_prc];
    ksdensity([x;y]','PlotFcn','contour');
    view([0 90]);
    plot( x, y, '.', 'Color', c);
    ylim([0 10]);
    xline(0,'-','LineWidth',1);
    xlabel('Distance from landmark (m)');
    ylabel('Field width (m)');
    title('field peak to nearest LM');
    
    pnl(ii_dir,2,2).select(); hold on
    x = [fields.start] - [fields.LM_prev_by_peak];
    y = [fields.width_prc];
    ksdensity([x;y]','PlotFcn','contour');
    view([0 90]);
    plot( x, y, '.', 'Color', c);
    ylim([0 10]);
    xline(0,'-','LineWidth',1);
    xlabel('Distance from landmark (m)');
    ylabel('Field width (m)');
    title('field start to previous LM');
    
    pnl(ii_dir,2,3).select(); hold on
    x = [fields.end] - [fields.LM_next_by_peak];
    y = [fields.width_prc];
    ksdensity([x;y]','PlotFcn','contour');
    view([0 90]);
    plot( x, y, '.', 'Color', c);
    ylim([0 10]);
    xline(0,'-','LineWidth',1);
    xlabel('Distance from landmark (m)');
    ylabel('Field width (m)');
    title('field end to next LM');
    
    %% link axes
    linkaxes([pnl(ii_dir,1,1).axis pnl(ii_dir,2,1).axis], 'x');
    linkaxes([pnl(ii_dir,1,2).axis pnl(ii_dir,2,2).axis], 'x');
    linkaxes([pnl(ii_dir,1,3).axis pnl(ii_dir,2,3).axis], 'x');
    
end
% save figure
figname = 'fields_size_position_vs_landmarks_with_correct_shuffle';
file_out = fullfile(dir_out,figname);
saveas(gcf, file_out, 'tif')
% saveas(gcf, file_out, 'pdf')


%% Fields size/position vs. landmarks
% =========================================================================
figure('Units','normalized','Position',[0 0 1 1])
pnl = panel();
pnl.pack('h',2);
pnl(1).pack('v',2);
pnl(2).pack('v',2);
pnl.margin = 30;
h=pnl.title('Fields size vs. inter-landmarks-distance');
h.FontSize = 16;
h.Position = [0.5 1.07];
LM = exp.LM;
LM( contains({LM.name},{'ball','enter'}) ) = [];
% LM( contains({LM.name},{'enter'}) ) = [];
for ii_dir = 1:2
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
    
    %% scatter plot with linear fit
    c = prm.graphics.colors.flight_directions{ii_dir};
    pnl(ii_dir,1).select()
    x = abs([fields.LM_prev_by_peak]-[fields.LM_next_by_peak]);
    y = [fields.width_prc];
    lm = fitlm(x,y);
    plot(x, y, '.', 'Color',c)
    h=plot(lm);
    h(1).Color = c;
    h(1).Marker = '.';
    text(0.1,0.9,sprintf('R^2=%f, pval=%f',lm.Rsquared.Ordinary,lm.Coefficients.pValue(2)),'Units','normalized');
    xlabel('Inter Landmark distance (m)')
    ylabel('Field Size (m)')
    title("direction"+ii_dir,'Units','normalized','Position',[0.5 1.07])
    
    %% violin plots
    pnl(ii_dir,2).select()
    violinplot(y,x)
    xlabel('Inter Landmark distance (m)')
    ylabel('Field Size (m)')
end
% save figure
figname = 'fields_size_vs_inter_LM_dist';
file_out = fullfile(dir_out,figname);
saveas(gcf, file_out, 'tif')


%% Fields locations cdf
% =========================================================================
figure('Units','normalized','Position',[0 0 1 1])
pnl = panel();
pnl.pack('h',2);
pnl.margin = 30;
h=pnl.title('Fields locations cdf');
h.FontSize = 16;
h.Position = [0.5 1.07];
LM = exp.LM;
% LM( contains({LM.name},{'ball','enter'}) ) = [];
% LM( contains({LM.name},{'enter'}) ) = [];
for ii_dir = 1:2
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
  
    % remove fields near balls
%     fields([fields.in_low_speed_area])=[];
%     fields(isnan([fields.LM_nearest_by_peak])) = [];
    
    %% plot cdf
    c = prm.graphics.colors.flight_directions{ii_dir};
    pnl(ii_dir).select();
    hold on
    h = cdfplot([fields.loc]);
    h.Color = c;
    h.LineWidth = 2;
    h = cdfplot([fields.start]);
    h.Color = c;
    h.LineWidth = 1;
    h.LineStyle = ':';
    h = cdfplot([fields.end]);
    h.Color = c;
    h.LineWidth = 1;
    h.LineStyle = '--';
    ha=gca;
    ha.GridLineStyle = 'none';
    plot_LM(LM)
    legend("field " + {"peak","start","end"} + " cdf",'Location','northwest');
    xlabel('Position (m)')
    ylabel('cdf')
    title("direction"+ii_dir,'Units','normalized','Position',[0.5 1.07])
end
% save figure
figname = 'fields_locations_cdf';
file_out = fullfile(dir_out,figname);
saveas(gcf, file_out, 'tif')

%% local function to plot lines for landmarks
function plot_LM(LM,ypos,angle)
if nargin<2
    ypos = 0.02;
    angle = 45;
end
if nargin<3
    angle = 45;
end
ylimits = get(gca,'ylim');
for ii_LM=1:length(LM)
    x = LM(ii_LM).pos_proj;
    name = LM(ii_LM).name;
    xline(x, '-', 'color', 0.9.*[1 1 1], 'LineWidth',0.5);
    text(x, ylimits(2)+ypos*diff(ylimits), LM(ii_LM).name, 'Rotation', 45, 'FontSize',8);
end
end




%%
