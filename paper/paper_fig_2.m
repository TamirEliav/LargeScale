%% Large Scale - fig 2 - Multi-scale spatial coding with many fields in individual dorsla-CA1 neurons

%%
clear 
clc

%% plotting options
field_speed_opt = 1;
% corr_type = 'pearson';
corr_type = 'spearman';

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
        offset_y = (jj-1)*4.3;
        offset = panel_A_pos + [offset_x offset_y];
        panel_A(ii,jj,1) = axes('position', [offset+[0 2.25] panel_A_size_FR_map]);
        panel_A(ii,jj,2) = axes('position', [offset+[0 1   ] panel_A_size_raster]);
        panel_A(ii,jj,3) = axes('position', [offset+[0 0   ] panel_A_size_raster]);
    end
end
% panel_A = permute(panel_A,[2 1 3]);
panel_A = panel_A(:,3:-1:1,:);
panel_A = reshape(panel_A,[9 3]);
% panel_A_legend = axes('position', [2.2 24.7 0.5 0.15]);

panels_size = [2 2];
panel_B    = axes('position', [ 2.0  8.5  panels_size           ]);
panel_C    = axes('position', [ 5.3  8.5  panels_size           ]);
panel_D(1) = axes('position', [ 8.6  8.5  panels_size.*[1 0.9]  ]);
panel_D(2) = axes('position', [ 8.6  8.5  panels_size.*[1 0.9]  ]);
panel_E = axes('position', [12.2  8.5  panels_size.*[1.3 1]  ]);
panel_F = axes('position', [16.5  8.5  panels_size.*[1.3 1]  ]);
panel_G = axes('position', [ 2.0  5 panels_size          ]);
panel_H = axes('position', [ 5.3  5 panels_size.*[1.4 1] ]);
panel_I = axes('position', [ 9.1  5 panels_size]          );
panel_J = axes('position', [12.5  5 panels_size.*[1.2 1] ]);
panel_K = axes('position', [16.2  4.4 3 3]);

%%
prm = PARAMS_GetAll();

%%
% delete(panel_A_legend);
% panel_A_legend = axes('position', [2.2 24.7 0.5 0.15]);

% axes(panel_A_legend);
% cla
% hold on
% box off
% set(gca,'Visible','off');
% plot([0 1],[1 1],'Color',prm.graphics.colors.flight_directions{1},'LineWidth',2);
% plot([0 1],[0 0 ],'Color',prm.graphics.colors.flight_directions{2},'LineWidth',2);
% text(1.2, 0.5, 'Detected fields');

%%
bat_ID_num_map = containers.Map([34 79 148 2289 9861],...
                                 1:5);

%% FR map + rasters - 9 examples
cell_examples = {
433;  56;  51;
609; 419; 477;
 57; 628; 337;
};
% other options: 57 474 658
if 1
for ii_cell = 1:length(cell_examples)
    %%
    cell_ID = cell_examples{ii_cell};
    cell = cell_load_data(cell_ID,'details','FR_map','fields','stats','FE');
    c = prm.graphics.colors.flight_directions;
    
    % map+fields
    axes(panel_A(ii_cell, 1));
    cla
    hold on
    maps=[cell.FR_map.all];
    x = maps(1).bin_centers;
    y = cat(1,maps.PSTH);
    m = round(max(y(:)));
    ylimits = [0 m+1];
    % low-speed area
    area_lowerval = ylimits(1) - 0.23*range(ylimits);
    area_upperval = ylimits(1) - 1e-2*range(ylimits);
    area([4 prm.fields.valid_speed_pos(1)]    , repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
    area([  prm.fields.valid_speed_pos(2) 194], repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
    h=plot(x,y);
    [h.Color] = disperse(c);
    box off
    h=gca;
    h.TickDir = 'out';
    h.XTick = [];
    h.YTick = [0 m];
    h.YLim = ylimits;
    h.XLim = [0 200];
    
    % fields
    dir_offsets = [-0.1 -0.17]+0.015;
    for ii_dir=1:2
        fields = cell.fields{ii_dir};
        if isfield(fields,'in_low_speed_area')
            fields([fields.in_low_speed_area])=[];
        end
        for ii_field = 1:length(fields)
            
            [xaf,yaf] = ds2nfu(fields(ii_field).edges_prc, repelem(dir_offsets(ii_dir)*range(h.YLim),2));
            annotation('line',xaf,yaf,'Linewidth', 2, 'Color', c{ii_dir});
        end
    end
    
    
    % cell details
    cell_num_str_pos_x   = [0.50 0.50 0.50 0.50 0.45 0.35 0.45 0.50 0.30];
    cell_num_str_pos_y   = [1.05 1.05 1.05 0.85 0.90 0.90 0.90 0.90 0.90];
    cell_stats_str_pos_x = [0.80 0.95 0.80 0.80 0.80 0.80 0.12 0.50 0.85];
    cell_stats_str_pos_y = [1.20 1.05 1.15 0.90 1.10 1.10 1.10 0.90 0.90]+0.05;
    text(cell_num_str_pos_x(ii_cell),...
         cell_num_str_pos_y(ii_cell),...
         sprintf('Cell %d',ii_cell),...
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
    % fields num (here to be above the 
    for ii_dir=1:2
        fields = cell.fields{ii_dir};
        if isfield(fields,'in_low_speed_area')
            fields([fields.in_low_speed_area])=[];
        end
        dir_offsets =2.05+[0.18 0];
        text(1.05, dir_offsets(ii_dir), num2str(length(fields)), 'Units','normalized', 'Color', c{ii_dir},...
            'HorizontalAlignment','right', 'VerticalAlignment','middle','FontSize',6);
    end
end

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

%% add zoom in panel (cell 3)
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

%% add zoom in panel (cell 7)
% choose cell/dir/field to zoom
ii_cell = 7;
ii_dir = 2;
ii_field = 4;
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
y(1,:) = 2*[1;1];
y(2,:) = 3*[1;1];
y(3,:) = 6*[1;1];
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
POSf = ds2nfu([x(end,1) 4 diff(x(end,:)) 8]);
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
arrow_y = repelem(0.9111,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.008,  'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y      ,  'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);

end







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
cells = cellfun(@(c)(cell_load_data(c,'details','meanFR')), {cells.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];
cells_details = [cells.details];
cells_ID = {cells_details.cell_ID};
meanFR = [cells.meanFR];
cells_ID([meanFR.all]>prm.inclusion.interneuron_FR_thr)=[];
clear cells stats cells_details cells_t meanFR
cells = cellfun(@(c)(cell_load_data(c,'details','stats','meanFR','stats','inclusion','signif','fields','FR_map','FE')), cells_ID, 'UniformOutput',0);
% cells = cellfun(@(c)(cell_load_data(c,'details')), cells_ID, 'UniformOutput',0);
cells = [cells{:}];

%% save loaded population data
% save( fullfile(res_dir,'Fig2_cells_data'), 'cells');


%% arragne population SI/sparsity/coverage (panels B,C,D)
signif = arrayfun(@(x)(x.TF), cat(1,cells.signif));
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

% count total fields coverage per cell per direction
total_area = nan(length(cells),2);
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    for ii_dir = 1:2
        if ~cell.signif(ii_dir).TF
            continue;
        end
        fields = cell.fields{ii_dir};
        fields([fields.in_low_speed_area]) = [];
        total_area(ii_cell, ii_dir) = sum([fields.width_prc]);
    end
end

%% panel B - spatial info histogram
axes(panel_B);
cla
hold on
text(-0.6,1.15, 'B', 'Units','normalized','FontWeight','bold');
h = histogram(SI);
h.FaceColor = 0.5*[1 1 1];
h.BinEdges = 0:0.5:6;
xlim([0 6])
ylim([0.7e0 1.4e2])
ha=gca;
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.YTick = [1 10 100];
ha.YTickLabel = {'10^{0}';'10^{1}';'10^{2}'};
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.1;
% ha.YScale = 'linear';
ha.YScale = 'log';
xlabel({'Spatial information';'(bits/spike)'}, 'Units','normalized','Position',[0.5 -0.17]);
ylabel('No. of cells', 'Units','normalized','Position',[-0.32 0.5])

x = SI;
hl=xline(nanmean(x)); hl.Color='r';
m = ha.YLim(2) + 0.15*range(ha.YLim);
plot(prctile(x,[25 75]), [m m], 'r-','LineWidth',1   ,'Clipping','off');
plot(prctile(x,[50]),    m    , 'r.','MarkerSize',10 ,'Clipping','off');

%% panel C - sparsity histogram
axes(panel_C);
cla
hold on
text(-0.45,1.15, 'C', 'Units','normalized','FontWeight','bold');
h = histogram(sparsity);
h.NumBins = 12;
h.FaceColor = 0.5*[1 1 1];
xlim([0 1])
ylim([0.7e0 1.4e2])
ha=gca;
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.YTick = [1 10 100];
ha.YTickLabel = {'10^{0}';'10^{1}';'10^{2}'};
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.1;
% ha.YScale = 'linear';
ha.YScale = 'log';
xlabel('Sparsity', 'Units','normalized','Position',[0.5 -0.17])
ylabel('No. of cells', 'Units','normalized','Position',[-0.32 0.5])

x = sparsity;
hl=xline(nanmean(x)); hl.Color='r';
m = ha.YLim(2) + 0.15*range(ha.YLim);
plot(prctile(x,[25 75]), [m m], 'r-','LineWidth',1   ,'Clipping','off');
plot(prctile(x,[50]),    m    , 'r.','MarkerSize',10 ,'Clipping','off');

%% panel D - Total area histogram
axes(panel_D(1));
cla
hold on
text(-0.45,1.275, 'D', 'Units','normalized','FontWeight','bold');

exp=exp_load_data(cell.details.exp_ID);
IX=find(contains({exp.LM.name},'ball'));
LM_locs = [exp.LM.pos_proj];
ball2ball_dist = diff(LM_locs(IX));
% total_area_L = ball2ball_dist;
total_area_L = diff(prm.fields.valid_speed_pos);

h = histogram(total_area(:));
h.NumBins = 14;
h.FaceColor = 0.5*[1 1 1];
xlim([0 0.9*total_area_L])
ylim([0.7e0 1.4e2])
ha=gca;
ha.XTick = [0:50:150];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.YTick = [1 10 100];
ha.YTickLabel = {'10^{0}';'10^{1}';'10^{2}'};
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.1;
% ha.YScale = 'linear';
ha.YScale = 'log';
xlabel('Coverage (m)', 'Units','normalized','Position',[0.5 -0.17])
ylabel('No. of cells', 'Units','normalized','Position',[-0.32 0.5])
ha.XLim(1) = 0;
ha.YLim(2) = 2e2;

x = total_area(:);
hl=xline(nanmean(x)); hl.Color='r';
m = ha.YLim(2) - 0.27*range(ha.YLim);
plot(prctile(x,[25 75]), [m m], 'r-','LineWidth',1   ,'Clipping','off');
plot(prctile(x,[50]),    m    , 'r.','MarkerSize',10 ,'Clipping','off');

fprintf( 'Average total area in meters : %.4g\n\r',nanmean(total_area(:)) )
fprintf( 'Average total area in prc (%%): %.4g\n\r',100*nanmean(total_area(:)) / total_area_L )
fprintf( 'Median total area in meters : %.4g\n\r',nanmedian(total_area(:)) )
fprintf( 'Median total area in prc (%%): %.4g\n\r',100*nanmedian(total_area(:)) / total_area_L )

% add normalized x-axis
axes(panel_D(2));
cla
hold on
hax = gca;
hax.XLim = 100 * panel_D(1).XLim / total_area_L;
hax.XTick = [0:30:90];
hax.XAxisLocation = 'top';
hax.YAxisLocation = 'right';
hax.Color = 'none';
hax.XColor = 'k';
hax.YColor = 'none';
box off
hax.XRuler.TickLabelGapMultiplier = -0.35;
hax.TickLength = panel_D(1).TickLength;
xlabel('Coverage (%)', 'Units','normalized','Position',[0.45 1.2]);

%% panel E - stability
axes(panel_E);
cla
hold on
text(-0.4,1.15, 'E', 'Units','normalized','FontWeight','bold');
% arrange data
signif = arrayfun(@(x)(x.TF), cat(1,cells.signif));
stats = [cells.stats];
stats_dir = cat(1,stats.dir);
stats_dir = stats_dir(signif);
[r c ~]=find(signif);
stats_all_dir = [stats(r).all];
% plot
x = [stats_dir.corr_odd_even];
h = histogram(x);
panel_E_num_bins = 20;
h.BinLimits = [-1 1];
h.NumBins = panel_E_num_bins;
% h.BinEdges = [-1:0.15:1];
% h.BinEdges = linspace(-1,1,27);
h.FaceColor = 0.5*[1 1 1];
h.Normalization = 'count';

ha= gca;
ha.XLim = [-1 1];
ha.XTick = -1:0.5:1;
ha.TickDir = 'out';
ha.TickLength = [0.03 0.03];
ha.YLim = [0 200];
ha.YLim = [0 1.1*max(h.Values)];
ha.YTick = [0:50:150];
ha=gca;
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel({'Map correlation'}, 'Units','normalized','Position',[0.5 -0.17]);
ylabel('No. of cells', 'Units','normalized','Position',[-0.24 0.5]);
title('Stability');

hl=xline(nanmean(x)); hl.Color='r';
m = ha.YLim(2) + 0.05*range(ha.YLim);
plot(prctile(x,[25 75]), [m m], 'r-','LineWidth',1   ,'Clipping','off');
plot(prctile(x,[50]),    m    , 'r.','MarkerSize',10 ,'Clipping','off');
fprintf('\tMap correlations odd vs. even flights: mean=%.2f median=%.2f IQR=%.2f-%.2f\n', nanmean(x), prctile(x,[50]), prctile(x,[25 75]) );

%%
% figure
% subplot(121); hold on
% histogram([stats_dir.corr_begin_end])
% histogram([stats_dir.corr_odd_even])
% legend('Corr (begin,end)','Corr (odd,even)','Location','northwest')
% xlabel('Correlation')
% ylabel('Counts')
% title('stability hist')
% subplot(122); hold on
% plot([stats_dir.corr_begin_end],[stats_all_dir.field_ratio_LS],'.')
% plot([stats_dir.corr_odd_even],[stats_all_dir.field_ratio_LS],'.')
% legend('Corr (begin,end)','Corr (odd,even)','Location','northwest')
% xlabel('Stability')
% ylabel('ratio (L/S)')
% title('Ratio L/S vs. stability')

%% panel F - map correlations histogram
% figure
axes(panel_F);
cla
hold on
text(-0.45,1.15, 'F', 'Units','normalized','FontWeight','bold');

% arrange data
signif = arrayfun(@(x)(x.TF), cat(1,cells.signif));
FR_maps_all = cat(1,cells.FR_map);
FR_maps_all = reshape([FR_maps_all.all],size(FR_maps_all,1),size(FR_maps_all,2),[]);
M = cat(1,FR_maps_all.PSTH);
M = reshape(M,size(FR_maps_all,1),size(FR_maps_all,2),[]);
signif = repmat(signif,1,1,size(M,3));
M(~signif) = nan;
pos_bins = cells(1).FR_map(1).all.bin_centers;
invalid_pos_IX = pos_bins<=prm.fields.valid_speed_pos(1) | pos_bins>=prm.fields.valid_speed_pos(2);
M(:,:,invalid_pos_IX) = nan;
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

data(isnan(data))=[];
shuffle(isnan(shuffle))=[];

% plot
nBinEdges = 21;
edges = linspace(-1,1,nBinEdges);
histogram(data,    'Normalization','pdf','BinEdges',edges,'FaceColor', 0.5*[1 1 1]);
histogram(shuffle, 'Normalization','pdf','BinEdges',edges,'DisplayStyle','stairs','EdgeColor','k','LineWidth',1.5);
[~,P_KS,KS_stat] = kstest2(data, shuffle);
P_RankSum = ranksum(data, shuffle);
text(1.1,0.8, sprintf('P_{KS} = %.02f',P_KS),'Units','normalized','FontSize',7,'HorizontalAlignment','right');
% text(1,0.9, sprintf('P=%.02f',P_RankSum),'Units','normalized','FontSize',7,'HorizontalAlignment','right');

fprintf('panel F (directionality)\n');
fprintf('two-sample KS test comparing map correlation (data vs shuffle)\n');
fprintf('P = %.2f, KS_stat=%.2f (n_data=%d,n_shuffle=%d)\n',P_KS,KS_stat,length(data),length(shuffle));

ha= gca;
ha.XLim = [-1 1];
ha.YLim = [0 5.2];
ha.XTick = -1:0.5:1;
ha.YTick = [0 2 4];
ha.TickDir = 'out';
ha.TickLength = [0.03 0.03];
ha=gca;
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel('Map correlation', 'Units','normalized','Position',[0.5 -0.17])
ylabel({'Probability';'density function'}, 'Units','normalized','Position',[-0.17 0.5])
title('Directionality');

%% panel G - field count histogram
% figure
axes(panel_G);
cla
hold on
text(-0.6,1.15, 'G', 'Units','normalized','FontWeight','bold');
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
%     h = histogram(nFields(ii_dir,:));
%     h.FaceColor = prm.graphics.colors.flight_directions{ii_dir};
%     nBinEdges = 17;
%     h.BinEdges = linspace(0,35,nBinEdges);
end
h = histogram(nFields(:));
h.FaceColor = 0.5*[1 1 1];
% nBinEdges = 12;
% h.BinEdges = linspace(0,35,nBinEdges);
h.BinEdges = 0.5+[0:20];
h.Data(h.Data > h.BinLimits(2)) = h.BinLimits(2);
xlabel({'No. of fields per direction'},'Units','normalized','Position',[0.5 -0.18]);
ylabel('No. of cells','Units','normalized','Position',[-0.32 0.5])
ha = gca;
ha.YScale = 'log';
ha.YLim = [0.7e0 1.4e2];
ha.XLim = [0 h.BinLimits(2)];
ha.XTick = [0:10:30];
ha.YTick = [1 10 100];
ha.YTickLabel = {'10^{0}';'10^{1}';'10^{2}'};
% ha.XLim = [0 35];
% ha.YLim = [0 40];
% ha.XTick = [0:5:35];
% ha.YTick = [0 40];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

x = nFields(:);
hl=xline(nanmean(x)); hl.Color='r';
m = ha.YLim(2) + 0.15*range(ha.YLim);
plot(prctile(x,[25 75]), [m m], 'r-','LineWidth',1   ,'Clipping','off');
plot(prctile(x,[50]),    m    , 'r.','MarkerSize',10 ,'Clipping','off');
fprintf('\tNo. of fields: mean=%.1f median=%.1f IQR=%.1f-%.1f\n', nanmean(x), prctile(x,[50]), prctile(x,[25 75]) );

%% panel H - field size histogram
% figure
axes(panel_H);
cla
hold on
text(-0.36,1.15, 'H', 'Units','normalized','FontWeight','bold');
fields_size = [];
fields_size_smallest1 = [];                     % per-direction
fields_size_smallest2 = nan(length(cells),2);   % per-cell
for ii_dir = 1:2
    for ii_cell = 1:length(cells)
        cell = cells(ii_cell);
        if ~cell.signif(ii_dir).TF % check signif per direction
            continue;
        end
        fields = cell.fields{ii_dir};
        fields([fields.in_low_speed_area]) = []; % remove fields in low speed area
        fields_size = [fields_size fields.width_prc];
        fields_size_smallest1 = [fields_size_smallest1 min([fields.width_prc])];
        fields_size_smallest2(ii_cell,ii_dir) = min([fields.width_prc]);
    end
end
h = histogram(fields_size);
h.FaceColor = 0.5*[1 1 1];
% h.NumBins = 33;
h.BinWidth = 1;
ha=gca;
ha.YScale = 'log';
xlabel('Field size (m)')
ylabel('No. of fields','Units','normalized','Position',[-0.23 0.5])
ha = gca;
% ha.XLim = [0 35];
ha.YLim = [0.8 350];
% ha.XTick = [0:5:35];
ha.YTick = [1 10 100];
ha.YTickLabel = {'10^{0}';'10^{1}';'10^{2}'};
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

x = fields_size;
hl=xline(nanmean(x)); hl.Color='r';
m = ha.YLim(2) + 0.15*range(ha.YLim);
plot(prctile(x,[25 75]), [m m], 'r-','LineWidth',1   ,'Clipping','off');
plot(prctile(x,[50]),    m    , 'r.','MarkerSize',10 ,'Clipping','off');
fprintf('\t Fields sizes: mean=%.1fm median=%.1fm IQR=%.1f-%.1fm\n', nanmean(x), prctile(x,[50]), prctile(x,[25 75]) );

save( fullfile(res_dir,'pop_dist_fields_size'), 'fields_size');

% report smallest field distribution (including single-field cells)
fields_size_smallest2 = min(fields_size_smallest2,[],2);
fields_size_smallest2(isnan(fields_size_smallest2))=[];
[fields_size_smallest_gamma_fit1] = gamfit(fields_size_smallest1);
[fields_size_smallest_gamma_fit2] = gamfit(fields_size_smallest2);
mu1 = prod(fields_size_smallest_gamma_fit1);
mu2 = prod(fields_size_smallest_gamma_fit2);
CV1 = 1/sqrt(fields_size_smallest_gamma_fit1(1));
CV2 = 1/sqrt(fields_size_smallest_gamma_fit2(1));
fprintf('smallest fields gamma distribution fit (per-direction):\n')
fprintf('\t Shape parameter kappa = %.2f\n', fields_size_smallest_gamma_fit1(1));
fprintf('\t Scale parameter theta = %.2fm\n', fields_size_smallest_gamma_fit1(2));
fprintf('\t mu = %.2f\n', mu1);
fprintf('\t CV = %.2f\n', CV1);
fprintf('smallest fields gamma distribution fit (per-cell):\n')
fprintf('\t Shape parameter kappa = %.2f\n', fields_size_smallest_gamma_fit2(1));
fprintf('\t Scale parameter theta = %.2fm\n', fields_size_smallest_gamma_fit2(2));
fprintf('\t mu = %.2f\n', mu2);
fprintf('\t CV = %.2f\n', CV2);

% figure
% hold on
% xxx = linspace(0,15,100);
% plot(xxx,gampdf(xxx,fields_size_smallest_gamma_fit1(1),fields_size_smallest_gamma_fit1(2)))
% plot(xxx,gampdf(xxx,fields_size_smallest_gamma_fit2(1),fields_size_smallest_gamma_fit2(2)))
% legend({sprintf('per-direction (CV=%.2f)',CV1);...
%         sprintf('per-cell      (CV=%.2f)',CV2)});


%% panel I - smallest & largest field size
% figure
axes(panel_I);
cla
hold on
text(-0.45,1.15, 'I', 'Units','normalized','FontWeight','bold');
LS_field_size = nan(2,length(cells));
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    if ~any([cell.signif.TF]) % at least one direction is significant
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
ylabel('Field size (m)','Units','normalized','Position',[-0.25 0.5])

%% panel J - field ratio (largest/smallest)
axes(panel_J);
cla
hold on
text(-0.4,1.15, 'J', 'Units','normalized','FontWeight','bold');
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
% LS_field_ratio_dir = LS_field_ratio_dir(:);
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
ha.XScale = 'log';
ha.XLim = [1 27];
% ha.YLim = [0 10];
% ha.YLim = [7e-1 120];
ha.YLim = [0.7 max(h(2).Values)*1.1];
% ha.XTick = [0:5:35];
ha.YTick = [1 10 100];
ha.XTick = [1 2 5 10 20];
ha.YTickLabel = {'10^{0}';'10^{1}';'10^{2}'};
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.17]);
ylabel('No. of cells','Units','normalized','Position',[-0.24 0.5])

x = LS_field_ratio_all;
hl=xline(nanmean(x)); hl.Color='r';
m = ha.YLim(2) + 0.15*range(ha.YLim);
plot(prctile(x,[25 75]), [m m], 'r-','LineWidth',1   ,'Clipping','off');
plot(prctile(x,[50]),    m    , 'r.','MarkerSize',10 ,'Clipping','off');
fprintf('\tFields size ratio (L/S): mean=%.1f median=%.1f IQR=%.1f-%.1f\n', nanmean(x), prctile(x,[50]), prctile(x,[25 75]) );

%% panel K - prepare data
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

%% panel K - field size ratio vs. speed ratio (direct control for speed!)
% figure
axes(panel_K);
cla
hold on
text(-0.3,0.9667, 'K', 'Units','normalized','FontWeight','bold');
% arrange data
cells_signif = cat(1,cells.signif);
cells_signif = arrayfun(@(x)(x.TF), cells_signif);
signif_IX = any(cells_signif,2);
stats = [cells(signif_IX).stats];
stats_all = [stats.all];
% plot 
switch field_speed_opt 
    case 1 % median speed @ field peak
        x = abs([stats_all.field_ratio_LS_vel]);
    case 2 % median of spikes speed
        x = abs([stats_all.field_ratio_LS_vel2]);
    case 3 % mean of spikes speed
        x = abs([stats_all.field_ratio_LS_vel3]);
end
y = [stats_all.field_ratio_LS];
[r,r_pval] = corr(x',y','rows','pairwise','type','Pearson');
[rho,rho_pval] = corr(x',y','rows','pairwise','type','Spearman');
fprintf('field size ratio vs speed ratio corr (pearson): r=%.2f p=%.2f df=%d \n',r,r_pval,sum(~isnan(y))-2);
fprintf('field size ratio vs speed ratio corr (spearman): rho=%.2f p=%.2f df=%d \n',rho,rho_pval,sum(~isnan(y))-2);
[~,ttest_P,~,ttest_stats] = ttest(x,1);
[signtest_P,~,signtest_stats] = signtest(x,1);
fprintf('LS speed ratio: (ttest) P=%.2g t=%.2g df=%d sd=%.3g\n',ttest_P,ttest_stats.tstat,ttest_stats.df,ttest_stats.sd);
fprintf('LS speed ratio: (signtest) P=%.2g sign=%d zval=%.3g\n',signtest_P,signtest_stats.sign,signtest_stats.zval);
plot(x, y, '.k');
switch corr_type
    case 'pearson'
        text(0.1,1, {sprintf('r = %.2f',r);sprintf('P = %.2f',r_pval)}, ...
            'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);
    case 'spearman'
        text(0.52,1.06, {['{\rho}' sprintf(' = %.2f',rho)];sprintf('P = %.2f',rho_pval)}, ...
            'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);
end
xlabel('Speed ratio', 'Units','normalized', 'Position',[0.5 -0.12]);
ylabel('Field size ratio', 'units','normalized', 'Position',[-0.15 0.5]);
set(gca,'yscale','log');
set(gca,'ytick',[1 2 3 5 10 15 20]);
ylim([0.9 max(y)+1]);
xlim([0.7 1.3])
% xlim([min(x) max(x)])
% xlim([0 2])
ha = gca;
ha.TickDir='out';
ha.TickLength = [0.02 0.02];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.1;


%% print/save the figure
fig_name_out = fullfile(res_dir, sprintf('%s__nbin=%d',fig_name_str,panel_E_num_bins));
% fig_name_out = fullfile(res_dir, sprintf('%s__corr_%s_%d',fig_name_str,corr_type,field_speed_opt));
% fig_name_out = fullfile(res_dir, sprintf('%s__corr_%s_%d_paramset_%d',fig_name_str,corr_type,field_speed_opt,prm.parmaset));
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('======================================================');
disp('figure was successfully saved to pdf format');
disp('======================================================');








%% calc some more stats... (old version = pooling directions)
single_PF_cells_IX = find(nansum(nFields)==1);
single_PF_cells = cells(single_PF_cells_IX);

fields=cat(1,single_PF_cells.fields);
fields=fields(~isnan(nFields(:,single_PF_cells_IX)'));
fields=[fields{:}];
fields([fields.in_low_speed_area])=[];

single_PF_field_size_mean = mean([fields.width_prc]);
single_PF_field_size_std  = std([fields.width_prc]);

% print stats
disp('Old version = pooling directions')
fprintf('no. of single field cells = %d/%d (%.2f%%),\n with mean field size = %.1fm std=%.1fm\n',...
    length(single_PF_cells), ...
    length(cells), ...
    100 * length(single_PF_cells) / length(cells),...
    single_PF_field_size_mean,...
    single_PF_field_size_std);

%% calc some more stats... (new version = per direction)
single_PF_cells_IX = (nFields==1)';
signif_cells_IX = ~isnan(nFields)';
fields=cat(1,cells.fields);
fields = fields(single_PF_cells_IX);
fields=[fields{:}];
fields([fields.in_low_speed_area])=[];

single_PF_field_size_mean = mean([fields.width_prc]);
single_PF_field_size_std  = std([fields.width_prc]);

% print stats
disp('New version = per direction')
fprintf('no. of single field cells x directions = %d/%d (%.2f%%),\n with mean field size = %.1fm std=%.1fm\n',...
    sum(single_PF_cells_IX,'all'), ...
    sum(signif_cells_IX,'all'), ...
    100 * sum(single_PF_cells_IX,'all') / sum(signif_cells_IX,'all'),...
    single_PF_field_size_mean,...
    single_PF_field_size_std);

%% figure for Liora
if 0
cells_details=[cells.details];
cells_anatomy_PD_prc = [cells_details.TT_pos_proximodistal_prc];

figure('Units','centimeters', 'position', [2 2 20 20])
subplot(2,2,1)
hold on
nBinEdges = 9;
edges = logspace(0,log10(25),nBinEdges);
h=histogram(LS_field_ratio_all,edges); h.FaceColor = 0.5*[1 1 1];
h=histogram(LS_field_ratio_dir,edges); h.FaceColor = 'g';
ha=gca;
ha.XScale = 'log';
ha.YScale = 'log';
ha.XLim = [0 27];
ha.YLim(1) = 7e-1;
ha.YTick = [1 10 100];
ha.XTick = [1 2 5 10 20];
ha.YTickLabel = {'10^{0}';'10^{1}';'10^{2}'};
ha.TickDir='out';
xlabel('Largest/Smallest ratio')
ylabel('No. cells')
legend({'dir pooled';'per dir'})
subplot(2,2,2)
hold on
plot(LS_field_ratio_dir(1,:), nFields(1,:),'.', 'Color',prm.graphics.colors.flight_directions{1});
plot(LS_field_ratio_dir(2,:), nFields(2,:),'.', 'Color',prm.graphics.colors.flight_directions{2});
xlabel('Largest/Smallest ratio')
ylabel('No. fields')
legend({'dir1';'dir2'},'Location','northwest')
subplot(2,2,3)
plot(LS_field_ratio_all, sum(nFields,1), '.k');
legend('dir pooled')
xlabel('Largest/Smallest ratio')
ylabel('No. fields')
subplot(2,2,4)
plot(total_area, cells_anatomy_PD_prc+0.01*randn(size(cells_anatomy_PD_prc)), '.')
xlabel('Total area (m)')
ylabel('Proximo-distal axis (%)');
legend({'dir1';'dir2'},'Location','southeast')

fig_name_out = fullfile(res_dir, 'Fig_2__for_Liora');
saveas(gcf, fig_name_out, 'fig');
saveas(gcf, fig_name_out, 'pdf');
saveas(gcf, fig_name_out, 'tif');
% close(gcf)
end




%%
