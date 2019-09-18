%% Large Scale - fig 2 - Behavioral and neural recordings from bats fliying over large spatial scales.

%%
% TODO: make sure to use only included cells
% TODO: check conditions for which data we are taking in the histograms
% panels (B->end)

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_2';
fig_caption_str = 'Multi-scale spatial coding with many fields in individual dorsla-CA1 neurons';
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

pause(0.2); % workaround to solve matlab automatically changing the axes positions...

% create panels
panel_A_size_raster = [5 1];
panel_A_size_FR_map = [5 1];
panel_A_pos = [2 14];
panel_A = [];
for ii=1:3
    for jj=1:3
        offset_x = (ii-1)*6;
        offset_y = (jj-1)*4;
        offset = panel_A_pos + [offset_x offset_y];
        panel_A(ii,jj,1) = axes('position', [offset+[0 2.25] panel_A_size_FR_map]);
        panel_A(ii,jj,2) = axes('position', [offset+[0 1   ] panel_A_size_raster]);
        panel_A(ii,jj,3) = axes('position', [offset+[0 0   ] panel_A_size_raster]);
    end
end
% panel_A = permute(panel_A,[2 1 3]);
panel_A = panel_A(:,3:-1:1,:);
panel_A = reshape(panel_A,[9 3]);

panel_BCDE_size = [2 2];
panel_B = axes('position', [ 2  10.5 panel_BCDE_size ]);
panel_C = axes('position', [ 5  10.5 panel_BCDE_size ]);
panel_D = axes('position', [ 8  10.5 panel_BCDE_size ]);
panel_E = axes('position', [ 11 10.5 panel_BCDE_size ]);

panel_FGH_size = [2 2];
panel_F = axes('position', [2   7  panel_FGH_size ]);
panel_G = axes('position', [5   7  panel_FGH_size ]);
panel_H = axes('position', [8   7  panel_FGH_size ]);
panel_I = axes('position', [2 2.5 3 3]);
panel_J(1) = axes('position', [6.5 2.5 3.0 3.0]);
panel_J(2) = axes('position', [8.0 4.0 1.5 1.5]);


%%
prm = PARAMS_GetAll();

%% FR map + rasters - 9 examples
cell_examples = {
'b0079_d160928_TT2_SS01';
'b0034_d180312_TT4_SS01';
'b0034_d180312_TT3_SS02';
'b0148_d170718_TT4_SS03';
'b0148_d170608_TT4_SS01';
'b0079_d160925_TT2_SS03';
'b0148_d170626_TT4_SS01';
'b0034_d180310_TT4_SS04';
'b0034_d180310_TT4_SS04';
};
for ii_cell = 1:9
    cell_ID = cell_examples{ii_cell};
    cell = cell_load_data(cell_ID,'details','FR_map','fields','stats','FE');
    c = prm.graphics.colors.flight_directions;
    
    % map+fields
    axes(panel_A(ii_cell, 1));
    cla
    maps=[cell.FR_map.all];
    x = maps(1).bin_centers;
    y = cat(1,maps.PSTH);
    h=plot(x,y);
    [h.Color] = disperse(c);
    box off
    h=gca;
    h.TickDir = 'out';
    h.XTick = [];
    m = round(max(y(:)));
    h.YTick = [0 m];
    h.YLim = [0 m+1];
    h.XLim = [0 200];
    
    % fields
    for ii_dir=1:2
        fields = cell.fields{ii_dir};
        if isfield(fields,'in_low_speed_area')
            fields([fields.in_low_speed_area])=[];
        end
        for ii_field = 1:length(fields)
            dir_offsets = [-0.1 -0.2];
            [xaf,yaf] = ds2nfu(fields(ii_field).edges_prc, repelem(dir_offsets(ii_dir)*range(h.YLim),2));
            annotation('line',xaf,yaf,'Linewidth', 2, 'Color', c{ii_dir});
        end
    end
    
    % cell details
    cell_num_str_pos_list = [repelem(1.05,3) repelem(0.9,6)];
    cell_num_str_pos = cell_num_str_pos_list(ii_cell);
    text(0.5, cell_num_str_pos, "cell "+ii_cell,...
        'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
    text(1,1.00,...
        sprintf('%.1f    %.1f',...
                cell.stats.all.field_largest,...
                cell.stats.all.field_smallest),...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',7);
    
    % rasters
    FEs = [cell.FE];
    for ii_dir=1:2
        axes(panel_A(ii_cell, ii_dir+1));
        cla
        FE = FEs{ii_dir};
        x = [FE.spikes_pos];
        y = [FE.spikes_ts];
%         [FE.number2] = disperse(1:length(FE));
%         y = arrayfun(@(FE)(FE.number2*ones(1,FE.num_spikes)),FE,'UniformOutput',0);
%         y = [y{:}];
        plot(x,y,'.','Color',c{ii_dir},'MarkerSize',0.05);
        behave_session_ts = exp_get_sessions_ti(cell.details.exp_ID,'Behave');
        t0 = behave_session_ts(1);
        rescale_plot_data('y', [1e-6/60 t0]);
        box off
        h=gca;
        m = ceil(max((y-t0)*1e-6/60));
        h.YTick = [0 m];
        h.XLim = [0 200];
        h.YLim = [0 m]+[-1 1];
        switch ii_dir
            case 1
                h.XTick = [];
                h.YTickLabel = {'',num2str(m)};
            case 2
                h.XTick = 0:50:200;
                h.XRuler.TickLabelGapOffset = -2;
                h.YTickLabel = {'0',num2str(m)};
                h.TickDir = 'out';
        end
        
    end
end

%% add direction arrows
arrow_x = 0.1 +[0 0.05];
arrow_y = repelem(0.965,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.008,  'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y      ,  'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);

%% add x/y labels for specific panels
for ii = [7 8 9]
    axes(panel_A(ii, 3));
    xlabel('Position (m)', 'Units','normalized','Position',[0.5 -0.35]);
end
for ii = [1 4 7]
    axes(panel_A(ii, 3));
    ylabel('Time (min)',   'Units','normalized','Position',[-0.1 1]);
end
axes(panel_A(1, 1));
text(-0.15,1.15, 'A', 'Units','normalized','FontWeight','bold');



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% load population data
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

%% panel B - field count histogram
% figure
axes(panel_B);
hold on
text(-0.15,1.15, 'B', 'Units','normalized','FontWeight','bold');
nFields = nan(2,length(cells));
for ii_dir = 1:2
    for ii_cell = 1:length(cells)
        cell = cells(ii_cell);
        % TODO: check signif definition!
        if ~cell.signif(ii_dir).TF % check signif per direction
            continue;
        end
        % TODO: this include the fields in the low speed areas, make sure
        % we want this!
        nFields(ii_dir,ii_cell) = cell.stats.dir(ii_dir).field_num;
    end
    h = histogram(nFields(ii_dir,:));
    h.FaceColor = prm.graphics.colors.flight_directions{ii_dir};
end
xlabel({'No. of fields';'per direction'},'Units','normalized','Position',[0.5 -0.22]);
ylabel('Count')
ha = gca;
% ha.XLim = [0 35];
% ha.YLim = [0 40];
% ha.XTick = [0:5:35];
% ha.YTick = [0 40];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

%% panel C - field size histogram
% figure
axes(panel_C);
hold on
text(-0.15,1.15, 'C', 'Units','normalized','FontWeight','bold');
fields_size = [];
for ii_dir = 1:2
    for ii_cell = 1:length(cells)
        cell = cells(ii_cell);
        if ~cell.signif(ii_dir).TF % check signif per direction
            continue;
        end
        fields = cell.fields{ii_dir};
        fields([fields.in_low_speed_area]) = []; % remove fields in low speed area
        fields_size = [fields_size fields.width_prc];
    end
end
h = histogram(fields_size);
h.FaceColor = 0.5*[1 1 1];
ha=gca;
ha.YScale = 'log';
xlabel('Field Size (m)')
ylabel('Count')
ha = gca;
% ha.XLim = [0 35];
% ha.YLim = [0 40];
% ha.XTick = [0:5:35];
ha.YTick = [1 10 100];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

%% panel D - smallest / largest field size
% figure
axes(panel_D);
cla
hold on
text(-0.15,1.15, 'D', 'Units','normalized','FontWeight','bold');
LS_field_size = nan(2,length(cells));
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    if ~all([cell.signif.TF]) % at least one direction is significant
        continue;
    end
    % TODO: what to do if we have only 1 field?!?!
    % TODO: note that for this panel we pool directions, but for the next
    % one we also separate directions
    LS_field_size(1,ii_cell) = cell.stats.all.field_smallest;
    LS_field_size(2,ii_cell) = cell.stats.all.field_largest;
end
% LS_field_size(:,any(isnan(LS_field_size))) = [];
rng(1); % fix seed for violin plot randomization
hv = violinplot(LS_field_size');
[hv.ViolinColor] = disperse(repelem({'none'},length(hv)));
hs=[hv.ScatterPlot];
[hs.Marker]          = disperse(repelem({'.'},length(hs)));
[hs.SizeData]        = disperse(repelem(35,length(hs)));
[hs.MarkerFaceColor] = disperse(repelem({0.5*[1 1 1]},length(hs)));
[hs.MarkerEdgeColor] = disperse(repelem({0.5*[1 1 1]},length(hs)));
hm=[hv.MedianPlot];
% remove default median plot (circle), and add a red bar instead
[hm.Marker] = disperse(repelem({'none'},2));
plot([hm.XData]+0.075*[-1 1]', repelem([hm.YData],2,1), 'r', 'LineWidth',1.5);
ha=gca;
ylimits = ha.YLim;
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XTick = [1 2];
ha.XTickLabel = {};
text([1 2],repelem(ylimits(1)-0.04*diff(ylimits),2),{{'Smallest';'field'},{'Largest';'field'}},...
    'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7);
% xlabel('')
ylabel('Field size (m)')

%% panel E - field ratio (largest/smallest)
% figure
axes(panel_E);
cla
hold on
text(-0.15,1.15, 'E', 'Units','normalized','FontWeight','bold');
LS_field_ratio_all = nan(1,length(cells));
LS_field_ratio_dir = nan(2,length(cells));
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    if any([cell.signif.TF]) % check at least one direction is signif
        LS_field_ratio_all(ii_cell) = cell.stats.all.field_ratio_LS;
    end
    for ii_dir = 1:2
        if cell.signif(ii_dir).TF % check signif per direction
            LS_field_ratio_dir(:,ii_cell) = cell.stats.dir(ii_dir).field_ratio_LS;
        end
    end
end
LS_field_ratio_dir = LS_field_ratio_dir(:);
LS_field_ratio_dir(isnan(LS_field_ratio_dir)) = [];
LS_field_ratio_all(isnan(LS_field_ratio_all)) = [];
edges = linspace(1,ceil(max(max(LS_field_ratio_all))),9);
clear h
h(1) = histogram(LS_field_ratio_dir);
h(2) = histogram(LS_field_ratio_all);
h(1).BinEdges = edges;
h(2).BinEdges = edges;
h(1).FaceColor = 'g';
h(2).FaceColor = 0.5*[1 1 1];
% h(1).BinLimits = [1 h(1).BinLimits([2])];
% h(2).BinLimits = [1 h(2).BinLimits([2])];
% h(1).NumBins = 15;
% h.FaceColor = prm.graphics.colors.flight_directions{ii_dir};
% legend({'per cell';'per direction'})
ha=gca;
% ha.YScale = 'log';
ha.YScale = 'linear';
% ha.XLim = [0 35];
ha.YLim = [0 10];
% ha.YLim = [0 260];
% ha.XTick = [0:5:35];
ha.YTick = [1 10 100];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel({'Fields ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.17]);
ylabel('Count')

%% arragne population SI/sparsity
signif = arrayfun(@(x)(x.TF), cat(1,cells.signif));
% TODO: check why SI/sparsity can be nan in the original cell calculation?
SI = arrayfun(@(x)([x.dir.SI_bits_spike]), cat(1,cells.stats),'UniformOutput',0);
sparsity = arrayfun(@(x)([x.dir.sparsity]), cat(1,cells.stats),'UniformOutput',0);
SI = cat(1,SI{:});
SI(~signif) = nan;
SI = SI(:);
SI(isnan(SI)) = [];
sparsity = cat(1,sparsity{:});
sparsity(~signif) = nan;
sparsity = sparsity(:);
sparsity(isnan(sparsity)) = [];

%% panel F - spatial info histogram
axes(panel_F);
cla
hold on
text(-0.15,1.15, 'F', 'Units','normalized','FontWeight','bold');
h = histogram(SI);
h.NumBins = 12;
h.FaceColor = 'k';
ha=gca;
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel({'Spatial information';'(bits/spike)'}, 'Units','normalized','Position',[0.5 -0.17]);
ylabel('Count')

%% panel G - sparsity histogram
axes(panel_G);
cla
hold on
text(-0.15,1.15, 'G', 'Units','normalized','FontWeight','bold');
h = histogram(sparsity);
h.NumBins = 10;
h.FaceColor = 'k';
ha=gca;
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel('Sparsity', 'Units','normalized','Position',[0.5 -0.17])
ylabel('Count', 'Units','normalized','Position',[-0.2 0.5])

%% panel H - map correlations histogram
axes(panel_H);
cla
hold on
text(-0.15,1.15, 'H', 'Units','normalized','FontWeight','bold');
maps_dir_corr = [];
for ii_cell = 1:length(cells)
        cell = cells(ii_cell);
        maps_dir_corr(ii_cell) = corr(...
             cell.FR_map(1).all.PSTH',...
             cell.FR_map(2).all.PSTH','rows', 'complete');
end
% maps_dir_corr(~all(signif')) = [];
h = histogram(maps_dir_corr);
h.FaceColor = 'k';
ha= gca;
ha.XLim = [-1 1];
ha.XTick = -1:0.5:1;
ha.TickDir = 'out';
ha.TickLength = [0.03 0.03];
ha=gca;
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel('Map correlation', 'Units','normalized','Position',[0.5 -0.17])
ylabel('Count', 'Units','normalized','Position',[-0.2 0.5])

%% panels I&J - prepare data
distances_all = [];
field_size_diff_all = [];
field_size_all = [];
field_vel_all = [];
field_vel2_all = [];
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    for ii_dir = 1:2
        if ~cell.signif(ii_dir).TF
            continue;
        end
        fields = cell.fields{ii_dir};
        distances_all = [distances_all diff([fields.loc])];
        field_size_diff_all = [field_size_diff_all diff([fields.width_prc])];
        field_size_all = [field_size_all [fields.width_prc]];
        field_vel_all = [field_vel_all [fields.vel]];
        field_vel2_all = [field_vel2_all [fields.vel2]];
    end
end

%% panel I - field size vs. speed
axes(panel_I);
cla
hold on
text(-0.2,1.1, 'I', 'Units','normalized','FontWeight','bold');
plot(abs(field_vel_all) , field_size_all ,'.k', 'MarkerSize',3)
% plot(abs(field_vel2_all) , field_size_all ,'.k')
ha = gca;
ha.XLim = [4 10];
ha.YLim = [0 35];
ha.XTick = [5 7 9];
ha.YTick = [0 10 20 30];
ha.TickDir='out';
ha.TickLength = [0.02 0.02];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel('Speed (m/s)', 'Units','normalized', 'Position',[0.5 -0.12]) 
ylabel('Field size (m)', 'units','normalized', 'Position',[-0.15 0.5])

%% panel J - field diff(size) vs. distance
axes(panel_J(1));
cla
hold on
text(-0.2,1.1, 'J', 'Units','normalized','FontWeight','bold');
plot(distances_all, abs(field_size_diff_all), 'k.' , 'MarkerSize',4);
ha = gca;
ha.TickDir='out';
ha.TickLength = [0.02 0.02];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel('Distance between fields (m)', 'Units','normalized', 'Position',[0.5 -0.12]) 
ylabel('{\Delta} Field size (m)', 'units','normalized', 'Position',[-0.15 0.5])

% panel J inset
axes(panel_J(2));
cla
hold on
plot(distances_all, abs(field_size_diff_all), 'k.' , 'MarkerSize',3);
ha = gca;
ha.XLim = [0 40];
ha.YLim = [0 15];
ha.XTick = ha.XLim;
ha.YTick = ha.YLim;
ha.TickDir='out';
ha.TickLength = [0.02 0.02];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');





%%








%%
