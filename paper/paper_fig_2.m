%% Large Scale - fig 2 - Multi-scale spatial coding with many fields in individual dorsla-CA1 neurons

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'Fig_2';
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
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

pause(0.2); % workaround to solve matlab automatically changing the axes positions...

% create panels
panel_A_size_raster = [5 1];
panel_A_size_FR_map = [5 1];
panel_A_pos = [2 12];
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

panel_BCD_size = [2 2];
panel_B = axes('position', [ 2.0  8.5  panel_BCD_size          ]);
panel_C = axes('position', [ 5.3  8.5  panel_BCD_size          ]);
panel_D = axes('position', [ 8.6  8.5  panel_BCD_size.*[1.3 1] ]);

panel_EFG_size = [2 2];
panel_E = axes('position', [ 2.0  5 panel_EFG_size]          );
panel_F = axes('position', [ 5.3  5 panel_EFG_size.*[1.4 1] ]);
panel_G = axes('position', [ 9.4  5 panel_EFG_size]          );

panel_H = axes('position', [13.0  8.5 panel_EFG_size.*[1.2 1] ]);


panel_I = axes('position', [13.0  4 3 3]);

%%
prm = PARAMS_GetAll();

%% FR map + rasters - 9 examples
cell_examples = {
433; 56; 51;
609; 67; 477;
419; 628; 337;
};
for ii_cell = 1:length(cell_examples)
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
            dir_offsets = [-0.1 -0.17]+0.015;
            [xaf,yaf] = ds2nfu(fields(ii_field).edges_prc, repelem(dir_offsets(ii_dir)*range(h.YLim),2));
            annotation('line',xaf,yaf,'Linewidth', 2, 'Color', c{ii_dir});
        end
    end
    
    % cell details
    cell_num_str_pos_x   = [0.50 0.50 0.50 0.50 0.50 0.45 0.50 0.50 0.50];
    cell_num_str_pos_y   = [1.05 1.05 1.05 0.85 0.90 0.90 0.90 0.90 0.90];
    cell_stats_str_pos_x = [0.80 0.95 0.80 0.80 0.20 0.80 0.80 0.50 0.85];
    cell_stats_str_pos_y = [1.20 1.05 1.15 0.90 0.95 1.10 1.10 0.90 0.90]+0.05;
    text(cell_num_str_pos_x(ii_cell), cell_num_str_pos_y(ii_cell), "cell "+ ii_cell,...
        'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
    switch ii_cell
        case {1,2,3,4,5,6,7,8}
            cell_stats_str = {  sprintf('max=%.1fm', cell.stats.all.field_largest);...
                                sprintf('min=%.1fm', cell.stats.all.field_smallest);...
                                sprintf('ratio=%.1f', cell.stats.all.field_ratio_LS);...
                             };
            text(cell_stats_str_pos_x(ii_cell), cell_stats_str_pos_y(ii_cell)-0*0.23, cell_stats_str{1},...
                'Units','normalized','HorizontalAlignment','center','VerticalAlignment','Top','FontSize',6);
            text(cell_stats_str_pos_x(ii_cell), cell_stats_str_pos_y(ii_cell)-1*0.23, cell_stats_str{2},...
                'Units','normalized','HorizontalAlignment','center','VerticalAlignment','Top','FontSize',6);
            text(cell_stats_str_pos_x(ii_cell), cell_stats_str_pos_y(ii_cell)-2*0.23, cell_stats_str{3},...
                'Units','normalized','HorizontalAlignment','center','VerticalAlignment','Top','FontSize',6);
        case 9
            cell_stats_str = {  sprintf('single field=%.1fm', cell.fields{1}.width_prc) };
            text(cell_stats_str_pos_x(ii_cell), cell_stats_str_pos_y(ii_cell)-0*0.23, cell_stats_str{1},...
                'Units','normalized','HorizontalAlignment','center','VerticalAlignment','Top','FontSize',6);
    end
    
    % rasters
    FEs = [cell.FE];
    for ii_dir=1:2
        axes(panel_A(ii_cell, ii_dir+1));
        cla
        FE = FEs{ii_dir};
        x = [FE.spikes_pos];
        [FE.number2] = disperse(1:length(FE));
        y = arrayfun(@(FE)(FE.number2*ones(1,FE.num_spikes)),FE,'UniformOutput',0);
        y = [y{:}];
        plot(x,y,'.','Color',c{ii_dir},'MarkerSize',0.05);
        box off
        h=gca;
        m = length(FE);
        h.YTick = [1 m];
        h.XLim = [0 200];
        h.YLim = [0 m+1];
        switch ii_dir
            case 1
                h.XTick = [];
                h.YTickLabel = {'',num2str(m)};
            case 2
                h.XTick = 0:50:200;
                h.XRuler.TickLabelGapOffset = -2;
                h.YTickLabel = {'1',num2str(m)};
                h.TickDir = 'out';
        end
        
    end
end

%% add zoom in panel
% choose cell/dir/field to zoom
ii_cell = 3;
ii_dir = 1;
ii_field = 2;
cell_ID = cell_examples{ii_cell};
cell = cell_load_data(cell_ID,'details','FR_map','fields','stats','FE');
c = prm.graphics.colors.flight_directions;
% add zoom ("out") lines
axes(panel_A(ii_cell,1));
x=[];
y=[];
x(1,:) = cell.fields{ii_dir}(ii_field).edges_href;
x(2,:) = cell.fields{ii_dir}(ii_field).edges_href;
% x(3,:) = [15;50];
x(3,:) = cell.fields{ii_dir}(ii_field).edges_href + 10*[-1 1];
y(1,:) =  4*[1;1];
y(2,:) = 7*[1;1];
y(3,:) = 14*[1;1];
[xf, yf] = ds2nfu(x,y);
for jj = 1:2
    for ii = 1:2
        hl = annotation('line');
%         hl.Parent = gca;
        hl.X = xf(ii+[0 1],jj);
        hl.Y = yf(ii+[0 1],jj);
        hl.LineWidth = 0.5;
        hl.Color = 0.5*[1 1 1];
    end
end
% create zoom panel
POSf = ds2nfu([x(end,1)-0.4 y(end,1)-2 diff(x(end,:)) 20]);
panel_A_zoom = axes('Units','normalized', 'position', POSf);
cla
hold on
set(gca,'visible','off');
FE = cell.FE{ii_dir};
x = [FE.spikes_pos];
[FE.number2] = disperse(1:length(FE));
y = arrayfun(@(FE)(FE.number2*ones(1,FE.num_spikes)),FE,'UniformOutput',0);
y = [y{:}];
xlimits = cell.fields{ii_dir}(ii_field).edges_href;
IX = find( x>xlimits(1) & x<xlimits(end) );
x=x(IX);
y=y(IX);
ylimits = [min(y) max(y)]+[-1 3];
plot(x,y,'.','Color',c{ii_dir},'MarkerSize',3);
% plot field size bar
x = cell.fields{ii_dir}(ii_field).edges_prc;
y = ylimits([end end]);
plot(x,y,'-','Color',c{ii_dir},'LineWidth',1);
text(mean(x),mean(y)+1, sprintf('%.1fm',diff(x)), 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',6);
box off
h=gca;
h.XLim = xlimits;
h.YLim = ylimits;


%% add direction arrows
arrow_x = 0.1 +[0 0.05];
arrow_y = repelem(0.8933,2);
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
%     ylabel('Time (min)',   'Units','normalized','Position',[-0.1 1]);
    ylabel('Flight no.',   'Units','normalized','Position',[-0.1 1]);
    axes(panel_A(ii, 1));
    ylabel({'Firing rate';'(Hz)'},   'Units','normalized','Position',[-0.07 0.42]);
end
axes(panel_A(1, 1));
text(-0.25,1.6, 'A', 'Units','normalized','FontWeight','bold');



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

%% panel E - field count histogram
% figure
axes(panel_E);
cla
hold on
text(-0.45,1.15, 'E', 'Units','normalized','FontWeight','bold');
nFields = nan(2,length(cells));
for ii_dir = 1:2
    for ii_cell = 1:length(cells)
        cell = cells(ii_cell);
        % check signif per direction
        if ~cell.signif(ii_dir).TF
            continue;
        end
        nFields(ii_dir,ii_cell) = cell.stats.dir(ii_dir).field_num;
    end
    h = histogram(nFields(ii_dir,:));
    h.FaceColor = prm.graphics.colors.flight_directions{ii_dir};
    nBinEdges = 17;
    h.BinEdges = linspace(0,35,nBinEdges);
end
xlabel({'No. of fields';'per direction'},'Units','normalized','Position',[0.5 -0.18]);
ylabel('No. of cells')
ha = gca;
ha.YScale = 'log';
ha.YLim = [0.7 100];
ha.XLim = [0 33];
ha.XTick = [0:10:30];
ha.YTick = [1 10 100];
ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
% ha.XLim = [0 35];
% ha.YLim = [0 40];
% ha.XTick = [0:5:35];
% ha.YTick = [0 40];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

%% panel F - field size histogram
% figure
axes(panel_F);
cla
hold on
text(-0.35,1.15, 'F', 'Units','normalized','FontWeight','bold');
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
% h.NumBins = 33;
h.BinWidth = 1;
ha=gca;
ha.YScale = 'log';
xlabel('Field Size (m)')
ylabel('Counts','Units','normalized','Position',[-0.23 0.5])
ha = gca;
% ha.XLim = [0 35];
% ha.YLim = [0 40];
% ha.XTick = [0:5:35];
ha.YTick = [1 10 100];
ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

%% panel G - smallest / largest field size
% figure
axes(panel_G);
cla
hold on
text(-0.45,1.15, 'G', 'Units','normalized','FontWeight','bold');
LS_field_size = nan(2,length(cells));
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    if ~all([cell.signif.TF]) % at least one direction is significant
        continue;
    end
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
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.001;
text([1 2],repelem(ylimits(1)-0.04*diff(ylimits),2),{{'Smallest';'field'},{'Largest';'field'}},...
    'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7);
% xlabel('')
ylabel('Field size (m)','Units','normalized','Position',[-0.21 0.5])



%% panel H - field ratio (largest/smallest)
% figure
axes(panel_H);
cla
hold on
text(-0.4,1.15, 'H', 'Units','normalized','FontWeight','bold');
LS_field_ratio_all = nan(1,length(cells));
LS_field_ratio_dir = nan(2,length(cells));
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    % pooled stats - check at least one direction is signif
    if any([cell.signif.TF])
        LS_field_ratio_all(ii_cell) = cell.stats.all.field_ratio_LS;
    end
    % per dir stats - check signif per direction
    for ii_dir = 1:2
        if cell.signif(ii_dir).TF 
            LS_field_ratio_dir(ii_dir,ii_cell) = cell.stats.dir(ii_dir).field_ratio_LS;
        end
    end
end
LS_field_ratio_dir = LS_field_ratio_dir(:);
% LS_field_ratio_dir(isnan(LS_field_ratio_dir)) = [];
% LS_field_ratio_all(isnan(LS_field_ratio_all)) = [];
% edges = linspace(1,ceil(max(max(LS_field_ratio_all))),9);
nBinEdges = 9;
edges = logspace(0,log10(25),nBinEdges);
clear h
% h(1) = histogram(LS_field_ratio_dir);
h(2) = histogram(LS_field_ratio_all);
% h(1).BinEdges = edges;
h(2).BinEdges = edges;
% h(1).FaceColor = 'g';
h(2).FaceColor = 0.5*[1 1 1];
% h(1).BinLimits = [1 h(1).BinLimits([2])];
% h(2).BinLimits = [1 h(2).BinLimits([2])];
% h(1).NumBins = 15;
% h.FaceColor = prm.graphics.colors.flight_directions{ii_dir};
% legend({'per cell';'per direction'})
ha=gca;
ha.YScale = 'log';
% ha.YScale = 'linear';
ha.XScale = 'log';
ha.XLim = [0 27];
% ha.YLim = [0 10];
ha.YLim = [7e-1 260];
% ha.XTick = [0:5:35];
ha.YTick = [1 10 100];
ha.XTick = [1 2 5 10 20];
ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.17]);
ylabel('No. of cells','Units','normalized','Position',[-0.24 0.5])

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

%% panel B - spatial info histogram
axes(panel_B);
cla
hold on
text(-0.45,1.15, 'B', 'Units','normalized','FontWeight','bold');
h = histogram(SI);
h.NumBins = 12;
h.FaceColor = 0.5*[1 1 1];
ha=gca;
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.1;
xlabel({'Spatial information';'(bits/spike)'}, 'Units','normalized','Position',[0.5 -0.17]);
ylabel('No. of cells')

%% panel C - sparsity histogram
axes(panel_C);
cla
hold on
text(-0.45,1.15, 'C', 'Units','normalized','FontWeight','bold');
h = histogram(sparsity);
h.NumBins = 15;
h.FaceColor = 0.5*[1 1 1];
ha=gca;
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.1;
xlabel('Sparsity', 'Units','normalized','Position',[0.5 -0.17])
ylabel('No. of cells', 'Units','normalized','Position',[-0.2 0.5])

%% panel D - map correlations histogram
% figure
axes(panel_D);
cla
hold on
text(-0.45,1.15, 'D', 'Units','normalized','FontWeight','bold');

% arrange data
signif = arrayfun(@(x)(x.TF), cat(1,cells.signif));
FR_maps_all = cat(1,cells.FR_map);
FR_maps_all = reshape([FR_maps_all.all],size(FR_maps_all,1),size(FR_maps_all,2),[]);
M = cat(1,FR_maps_all.PSTH);
M = reshape(M,size(FR_maps_all,1),size(FR_maps_all,2),[]);
signif = repmat(signif,1,1,size(M,3));
M(~signif) = nan;
ccc = corr(squeeze(M(:,1,:))', squeeze(M(:,2,:))' ,'rows', 'pairwise');
data = diag(ccc);
switch 2
    case 1 % compare all PSTH from other cells/dir
        M2 = reshape(M,size(M,1)*size(M,2),[]);
        ccc_shuffle = corr(M2' ,'rows', 'pairwise');
        mask = tril(true(size(ccc_shuffle)),-1);
        shuffle = ccc_shuffle(mask);
    case 2 % compare only different cells betweeb different directions
        mask = tril(true(size(ccc)),-1);
        shuffle = ccc(mask);
end

% plot
nBinEdges = 21;
edges = linspace(-1,1,nBinEdges);
histogram(data,    'Normalization','pdf','BinEdges',edges,'FaceColor', 0.5*[1 1 1]);
histogram(shuffle, 'Normalization','pdf','BinEdges',edges,'DisplayStyle','stairs','EdgeColor','k','LineWidth',1.5);
[~,P_KS] = kstest2(data, shuffle);
P_RankSum = ranksum(data, shuffle);
text(1,0.9, sprintf('P_{KS}=%.02f',P_KS),'Units','normalized','FontSize',7,'HorizontalAlignment','right');
% text(1,0.9, sprintf('P=%.02f',P_RankSum),'Units','normalized','FontSize',7,'HorizontalAlignment','right');

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
ylabel('Probability', 'Units','normalized','Position',[-0.2 0.5])

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

%% panel I - field size ratio vs. speed ratio (direct control for speed!)
% figure
axes(panel_I);
cla
hold on
text(-0.27,1.1, 'I', 'Units','normalized','FontWeight','bold');
% arrange data
cells_signif = cat(1,cells.signif);
cells_signif = arrayfun(@(x)(x.TF), cells_signif);
signif_IX = any(cells_signif,2);
stats = [cells(signif_IX).stats];
stats_all = [stats.all];
% plot 
x = abs([stats_all.field_ratio_LS_vel]);
y = [stats_all.field_ratio_LS];
[r,pval] = corr(x',y','rows','pairwise');
plot(x, y, '.k');
text(0.1,1, {sprintf('r=%.2f',r);sprintf('P=%.2f',pval)}, ...
    'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);
xlabel('Speed ratio', 'Units','normalized', 'Position',[0.5 -0.12]);
ylabel('Field size ratio', 'units','normalized', 'Position',[-0.15 0.5]);
set(gca,'yscale','log');
set(gca,'ytick',[1 2 3 5 10 15 20]);
ylim([0.9 max(y)+1]);
xlim([0.7 1.3])
ha = gca;
ha.TickDir='out';
ha.TickLength = [0.02 0.02];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.1;


%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');






%%











%%











%%
