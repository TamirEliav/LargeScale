%% Large Scale - fig 6 - Lab born vs Wild bats
%%
clear 
clc

%% plotting options
field_speed_opt = 1;
% corr_type = 'pearson';
corr_type = 'spearman';
match_TT_pos = false;

% examples_option = 5; %1:5
hist_normalization_option = 'probability'; %count
hist_Yscale_option = 'log';

color_Lab = [0    0.75    0];%[0.4940    0.1840    0.5560];
color_Wild = [0.3 0.25 0.2];%color_wild = [0.65 0.6 0.65];
color_wild_text = [0.45 0.4 0.45];
Lab_alpha = 0.35;
Wild_alpha = 0.6;

%% add data folders to path
addpath('L:\processed_data_structs');

%% define output files
res_dir =  'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'Fig_6';
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

pause(0.2); % workaround to solve matlab automatically changing the axes positions...

% create panels

% panel_A    = axes('position', [ 1.5  22.3  12  2]); %before comments from the authors
panel_A    = axes('position', [ -0.3  22.3  16  2.3]); 

panel_B_size_raster = [5 1];
panel_B_size_FR_map = [5 1];
panel_B_pos = [2 13.5];
panel_B = [];
for ii=1:3
    for jj=1:2
        offset_x = (ii-1)*6;
        offset_y = (jj-1)*4.3;
        offset = panel_B_pos + [offset_x offset_y];
        panel_B(ii,jj,1) = axes('position', [offset+[0 2.25] panel_B_size_FR_map]);
        panel_B(ii,jj,2) = axes('position', [offset+[0 1   ] panel_B_size_raster]);
        panel_B(ii,jj,3) = axes('position', [offset+[0 0   ] panel_B_size_raster]);
    end
end
% panel_A = permute(panel_A,[2 1 3]);
panel_B = panel_B(:,2:-1:1,:);
panel_B = reshape(panel_B,[6 3]);

lab_wild_legend = axes('position', [1.55 11.7 1.4 0.7]);

panels_size = [2 2];
panel_C    = axes('position', [ 2    9.5  panels_size           ]);
panel_D    = axes('position', [ 6.2  9.5  panels_size.*[1.4 1]  ]);
panel_E    = axes('position', [ 11.2    9.5  panels_size           ]);
panel_F = axes('position', [ 2.0  5.9 panels_size          ]);
panel_G = axes('position', [ 6.5  5.9 panels_size          ]);
panel_H = axes('position', [ 11.2  5.9 panels_size          ]);
panel_I = axes('position', [ 16  8.5  3  3 ]);
panel_I_legend = axes('position', [18.5 11.5 1 1]);
% panel_I_legend = axes('position', [19.4 10.5 1 1]);
% panel_I = axes('position', [ 1.6  2.5   panels_size.*[1.4 1]          ]);
% panel_J = axes('position', [ 6.25  2.5   panels_size.*[1.4 1]          ]);
% panel_K = axes('position', [ 10.8  2.5  panels_size.*[1.4 1]          ]);

inset_size = [1.1 1.1];
panel_C_inset    = axes('position', [ 3.2    11  inset_size           ]);
panel_D_inset    = axes('position', [ 7.6    11  inset_size.*[1.5 1]  ]);
panel_E_inset    = axes('position', [ 13   11  inset_size.*[1.2 1]           ]);

%%
prm = PARAMS_GetAll();

%% Block Diagram

% Block_diagram_filename = 'Block_diagram_revised_lioraVersion.PNG';
Block_diagram_filename = 'L:\resources\Block_diagram_revised_lioraVersion.PNG';
axes(panel_A);
image = imread(Block_diagram_filename);
imshow(image(:,:,1:3));
text(-0.035,0.9, 'A', 'Units','normalized','FontWeight','bold');


%% load population data - directly from the "per-cell" structs
% prm = PARAMS_GetAll();
% cells_t = DS_get_cells_summary();
% bats = [9845, 102, 194];
% cells_t(~ismember(cells_t.bat, bats ),:) = [];
% cells_t([cells_t.remove]==1,:) = [];%cells chosen to remove because of repititions between days
% cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
% cells = [cells{:}];
% cells = [cells.details];
% cells(~contains({cells.brain_area}, 'CA1')) = [];
% cells(~ismember([cells.ClusterQuality], [2])) = [];
% % remove cells based on FR:
% cells = cellfun(@(c)(cell_load_data(c,'details','meanFR')), {cells.cell_ID}, 'UniformOutput',0);
% cells = [cells{:}];
% cells_details = [cells.details];
% cells_ID = {cells_details.cell_ID};
% meanFR = [cells.meanFR];
% cells_ID([meanFR.all]>prm.inclusion.interneuron_FR_thr)=[];
% % remove non-signif cells:
% cells = cellfun(@(c)(cell_load_data(c,'details','signif')), cells_ID, 'UniformOutput',0);
% cells = [cells{:}];
% signif = cat(1,cells.signif);
% signif = arrayfun(@(x)(x.TF), signif);
% signif = any(signif,2);
% cells_ID(~signif) = [];
% 
% clear cells stats cells_details cells_t
% cells = cellfun(@(c)(cell_load_data(c,'details','stats','meanFR','inclusion','signif','fields','FR_map','FE')), cells_ID, 'UniformOutput',0);
% cells = [cells{:}];
% save('Cells_lab.mat','cells');

%% load the data - from the total pre-compiled big struct for all cells
cells_lab = load('cells_lab.mat');
cells_lab = cells_lab.cells;
cells_lab_IDs = arrayfun(@(x)(x.cell_ID), cat(1,cells_lab.details), 'UniformOutput', false);

cells_details = load('cells_details_updated_TT_pos.mat');
cells_details = cells_details.cells;
cells = load('cells_bat_200m.mat');
cells = cells.cells;
for ii_cell = 1:length(cells)
    cells(ii_cell).details = cells_details(ii_cell).details; %updated TT pos
end
cells_wild = cells;

if match_TT_pos
    % remove cells based on TT position:
    details_lab = cat(1,cells_lab.details);
    details_wild = cat(1,cells_wild.details);
    TT_pos_proximodistal_prc_lab = arrayfun(@(x)(x.TT_pos_proximodistal_prc), details_lab);
    TT_pos_proximodistal_prc_wild = arrayfun(@(x)(x.TT_pos_proximodistal_prc), details_wild);
    
    max_lab = prctile(TT_pos_proximodistal_prc_lab,98);
    min_lab = prctile(TT_pos_proximodistal_prc_lab,2);
    max_wild = prctile(TT_pos_proximodistal_prc_wild,98);
    min_wild = prctile(TT_pos_proximodistal_prc_wild,2);
    
    TT_desired_pos = [max(min_wild,min_lab),min(max_wild,max_lab)];
    
    IX_desired_pos_lab = and(TT_pos_proximodistal_prc_lab>=TT_desired_pos(1), ...
        TT_pos_proximodistal_prc_lab<=TT_desired_pos(2));
%     cells_ID(~IX_desired_pos_lab) = [];
%     cells = cellfun(@(c)(cell_load_data(c,'details','stats','meanFR','stats','inclusion','signif','fields','FR_map','FE')), cells_ID, 'UniformOutput',0);
%     cells_lab = [cells{:}];
    cells_lab(~IX_desired_pos_lab) = [];

    IX_desired_pos_wild = and(TT_pos_proximodistal_prc_wild>=TT_desired_pos(1), ...
        TT_pos_proximodistal_prc_wild<=TT_desired_pos(2));
    cells_wild(~IX_desired_pos_wild) = [];
end

%% FR map + rasters - 6 examples
cell_examples = {
            'b0194_d180516_TT4_SS04';
            'b0194_d180508_TT1_SS01';
            'b0194_d180611_TT4_SS01';
            'b0194_d180513_TT2_SS02';
            'b0102_d171012_TT4_SS01';
            'b0194_d180612_TT4_SS02'};

if 1
for ii_cell = 1:length(cell_examples)
    %%
    example_IX = find(strcmp(cells_lab_IDs,cell_examples{ii_cell}));
    cell = cells_lab(example_IX);
    
%     cell_ID = cell_examples{ii_cell};
%     cell = cell_load_data(cell_ID,'details','FR_map','fields','stats','FE');
    c = prm.graphics.colors.flight_directions;
    
    % map+fields
    axes(panel_B(ii_cell, 1));
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
    area([4 cell.stats.all.valid_speed_pos(1)]    , repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
    area([cell.stats.all.valid_speed_pos(2) 194], repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
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
    cell_num_str_pos_x   = [0.50 0.50 0.50 0.50 0.50 0.50];
    cell_num_str_pos_y   = [1.05 1.05 1.05 1.05 1.05 1.05];
    cell_stats_str_pos_x = [0.95 0.95 0.95 0.97 0.95 0.95];
    cell_stats_str_pos_y = [1.2 1.05 1.07 1.15  1.06 1.05]+0.05;
%     h = text(cell_num_str_pos_x(ii_cell), cell_num_str_pos_y(ii_cell), cell_ID,...
%         'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
%     h.Interpreter='none';
    text(cell_num_str_pos_x(ii_cell), cell_num_str_pos_y(ii_cell), "Cell "+ ii_cell,...
        'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
    switch ii_cell
        case {1,2,3,4,5,6,7,8}
            cell_stats_str = {  sprintf('max=%.1fm', cell.stats.all.field_largest);...
                                sprintf('min=%.1fm', cell.stats.all.field_smallest);...
                                sprintf('ratio=%.1f', cell.stats.all.field_ratio_LS);...
%                                 sprintf('Iso dis=%.1f', cell.stats.all.IsoDist);...
                             };
            text(cell_stats_str_pos_x(ii_cell), cell_stats_str_pos_y(ii_cell)-0*0.23, cell_stats_str{1},...
                'Units','normalized','HorizontalAlignment','center','VerticalAlignment','Top','FontSize',6);
            text(cell_stats_str_pos_x(ii_cell), cell_stats_str_pos_y(ii_cell)-1*0.23, cell_stats_str{2},...
                'Units','normalized','HorizontalAlignment','center','VerticalAlignment','Top','FontSize',6);
            text(cell_stats_str_pos_x(ii_cell), cell_stats_str_pos_y(ii_cell)-2*0.23, cell_stats_str{3},...
                'Units','normalized','HorizontalAlignment','center','VerticalAlignment','Top','FontSize',6);
%             text(cell_stats_str_pos_x(ii_cell), cell_stats_str_pos_y(ii_cell)-3*0.23, cell_stats_str{4},...
%                 'Units','normalized','HorizontalAlignment','center','VerticalAlignment','Top','FontSize',6);
        case 9
            cell_stats_str = {  sprintf('single field=%.1fm', cell.fields{1}.width_prc) };
            text(cell_stats_str_pos_x(ii_cell), cell_stats_str_pos_y(ii_cell)-0*0.23, cell_stats_str{1},...
                'Units','normalized','HorizontalAlignment','center','VerticalAlignment','Top','FontSize',6);
    end
    
    % rasters
    FEs = [cell.FE];
    for ii_dir=1:2
        axes(panel_B(ii_cell, ii_dir+1));
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
for ii = [4 5 6]
    axes(panel_B(ii, 3));
    xlabel('Position (m)', 'Units','normalized','Position',[0.5 -0.35]);
end
for ii = [1 4]
    axes(panel_B(ii, 3));
%     ylabel('Time (min)',   'Units','normalized','Position',[-0.1 1]);
    ylabel('Flight no.',   'Units','normalized','Position',[-0.1 1]);
    axes(panel_B(ii, 1));
    ylabel({'Firing rate';'(Hz)'},   'Units','normalized','Position',[-0.07 0.42]);
end
axes(panel_B(1, 1));
text(-0.25,1.6, 'B', 'Units','normalized','FontWeight','bold');


%% add direction arrows
arrow_x = 0.1 +[0 0.05];
arrow_y = repelem(0.815,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.008,  'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y      ,  'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);

end







% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% legend:
axes(lab_wild_legend);
cla
box off
hax = gca;
h_app = annotation('rectangle',[0 0.3 0.25 0.25],'color', 'k', 'FaceColor',color_Lab, 'FaceAlpha',Lab_alpha, 'LineWidth', 0.01);
h_app.Parent = hax;
h_dis = annotation('rectangle',[0 0.7 0.25 0.25],'color','k', 'FaceColor',0.5*[1 1 1], 'LineWidth', 0.01);
h_dis.Parent = hax;
text(0.32,0.45,'Lab', 'HorizontalAlignment','left','FontSize',7);
text(0.32,0.85,'Wild', 'HorizontalAlignment','left','FontSize',7);
hax.Visible = 'off';

%% panel C - field count histogram
axes(panel_C);
cla
hold on
text(-0.6,1.25, 'C', 'Units','normalized','FontWeight','bold');
nFields_lab = nan(2,length(cells_lab));
for ii_dir = 1:2
    for ii_cell = 1:length(cells_lab)
        cell = cells_lab(ii_cell);
        % check signif per direction
        if ~cell.signif(ii_dir).TF
            continue;
        end
        nFields_lab(ii_dir,ii_cell) = cell.stats.dir(ii_dir).field_num;
    end
end
nFields_lab = nFields_lab(:);
nFields_lab(isnan(nFields_lab)) = [];

nFields_wild = nan(2,length(cells_wild));
for ii_dir = 1:2
    for ii_cell = 1:length(cells_wild)
        cell = cells_wild(ii_cell);
        % check signif per direction
        if ~cell.signif(ii_dir).TF
            continue;
        end
        nFields_wild(ii_dir,ii_cell) = cell.stats.dir(ii_dir).field_num;
    end
end
nFields_wild = nFields_wild(:);
nFields_wild(isnan(nFields_wild)) = [];

h = histogram(nFields_wild(:));
h.FaceColor = 0.5*[1 1 1];
h.BinEdges = 0.5+[0:20];
h.Data(h.Data > h.BinLimits(2)) = h.BinLimits(2);
h.Normalization = hist_normalization_option;

h = histogram(nFields_lab(:));
h.FaceColor = color_Lab;
h.FaceAlpha = Lab_alpha;
h.BinEdges = 0.5+[0:20];
h.Normalization = hist_normalization_option;
h.Data(h.Data > h.BinLimits(2)) = h.BinLimits(2);
xlabel({'No. of fields per direction'},'Units','normalized','Position',[0.5 -0.18]);
ha = gca;
ha.YScale = hist_Yscale_option;
ha.XLim = [0 h.BinLimits(2)];
ha.XTick = [0:10:30];
% % ha.XLim = [0 35];
% % ha.YLim = [0 40];
% % ha.XTick = [0:5:35];
% % ha.YTick = [0 40];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
switch hist_normalization_option
    case 'count'
        ylim([0.7e0 1.4e2])
        ha.YTick = [1 10 100];
        ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
        ylabel('No. of cells', 'Units','normalized','Position',[-0.28 0.5])
    case 'probability'
        if strcmp(hist_Yscale_option,'log')
            ylim([0.001 0.5])
            ha.YTick = [0.01 0.1];
            ha.YTickLabel = {'10^{-2}';'10^{-1}'};
        else
            ylim([0 0.2])
            ha.YTick = [0 0.2];
        end
        ylabel('Fraction of cells', 'Units','normalized','Position',[-0.33 0.5])
end
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

%% panel C Inset
axes(panel_C_inset);
cla
hold on
h = histogram(nFields_wild(:));
h.FaceColor = 0.5*[1 1 1];
h.BinEdges = 0.5+[0:20];
h.Data(h.Data > h.BinLimits(2)) = h.BinLimits(2);
h.Normalization = hist_normalization_option;

h = histogram(nFields_lab(:));
h.FaceColor = color_Lab;
h.FaceAlpha = Lab_alpha;
h.BinEdges = 0.5+[0:20];
h.Normalization = hist_normalization_option;
h.Data(h.Data > h.BinLimits(2)) = h.BinLimits(2);
ha = gca;
ha.YScale = 'Linear';
ha.XLim = [0 h.BinLimits(2)];
ha.XTick = [0 10 20];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ylim([0 0.2])
ha.YTick = [0 0.1 0.2];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

%% panel F - field count bars
axes(panel_F);
cla
hold on
text(-0.6,1.15, 'F', 'Units','normalized','FontWeight','bold');

% [~,P_KS] = kstest2(nFields_lab(:), nFields_wild(:));
% text(1.5,0.8, sprintf('P_{KS} = %.02f',P_KS),'Units','normalized','FontSize',7,'HorizontalAlignment','right');
% [~,P_ttest] = ttest2(nFields_lab(:), nFields_wild(:));
[P_wilcoxon,~,wilcoxon_stats] = ranksum(nFields_lab(:), nFields_wild(:));

nFields_lab_noNans = nFields_lab(:); nFields_lab_noNans(isnan(nFields_lab_noNans)) = [];
nFields_wild_noNans = nFields_wild(:); nFields_wild_noNans(isnan(nFields_wild_noNans)) = [];
% text(1,0.8, sprintf('N = %d',length(nFields_lab_noNans)),'Units','normalized','FontSize',7,'HorizontalAlignment','right','color',color_Lab);
% text(1,0.7, sprintf('N = %d',length(nFields_wild_noNans)),'Units','normalized','FontSize',7,'HorizontalAlignment','right','color',0.5*[1 1 1]);

% % bar version:
% x = [1,2];
% y = [mean(nFields_lab_noNans(:)),mean(nFields_wild_noNans(:))];
% hb = bar(x,y);
% hb.FaceColor = 'flat';
% hb.CData = [color_Lab; 0.5*[1 1 1]];
% sem_lab = std(nFields_lab_noNans(:))/sqrt(length(nFields_lab_noNans(:)));
% sem_wild = std(nFields_wild_noNans(:))/sqrt(length(nFields_wild_noNans(:)));
% 
% he=errorbar(x,y,[sem_lab sem_wild]);
% he.CapSize = 2;
% he.LineStyle = 'none';
% he.Color = 'k';
% h=gca;
% h.XTick = [1 2];
% h.XTickLabels = {'Lab','Wild'};
% h.TickLength = [0.03 0.03];

% boxplot version:
plot_boxplot(nFields_lab_noNans,nFields_wild_noNans,gca,[0 0.75 0],[0.5 0.5 0.5])
h = gca;
h.XTickLabels = {'Lab','Wild'};
h.TickLength = [0.03 0.03];


ylabel({'No. of fields', 'per direction'}, 'Units','normalized','Position',[-0.22 0.5])
h.XRuler.TickLabelGapMultiplier = -0.2;
h.YRuler.TickLabelGapMultiplier = -0.05;

text(1.4,0.85, sprintf('P = %.02f',P_wilcoxon),'Units','normalized','FontSize',7,'HorizontalAlignment','right');

%% panel D - field size histogram
axes(panel_D);
cla
hold on
text(-0.37,1.25, 'D', 'Units','normalized','FontWeight','bold');
fields_size_lab = [];
for ii_dir = 1:2
    for ii_cell = 1:length(cells_lab)
        cell = cells_lab(ii_cell);
        if ~cell.signif(ii_dir).TF % check signif per direction
            continue;
        end
        fields = cell.fields{ii_dir};
        fields([fields.in_low_speed_area]) = []; % remove fields in low speed area
        fields_size_lab = [fields_size_lab fields.width_prc];
    end
end
fields_size_wild = [];
for ii_dir = 1:2
    for ii_cell = 1:length(cells_wild)
        cell = cells_wild(ii_cell);
        if ~cell.signif(ii_dir).TF % check signif per direction
            continue;
        end
        fields = cell.fields{ii_dir};
        fields([fields.in_low_speed_area]) = []; % remove fields in low speed area
        fields_size_wild = [fields_size_wild fields.width_prc];
    end
end
h = histogram(fields_size_wild);
h.FaceColor = 0.5*[1 1 1];
% h.NumBins = 33;
h.BinWidth = 1;
h.Normalization = hist_normalization_option;
h = histogram(fields_size_lab);
h.FaceColor = color_Lab;
h.FaceAlpha = Lab_alpha;
% h.NumBins = 33;
h.BinWidth = 1;
h.Normalization = hist_normalization_option;
ha=gca;
ha.YScale = hist_Yscale_option;
xlabel('Field size (m)','Units','normalized','Position',[0.5 -0.18]);
ha = gca;
% % ha.XLim = [0 35];
% ha.YLim = [0.8 300];
% % ha.XTick = [0:5:35];
switch hist_normalization_option
    case 'count'
        ylim([0.8 300])
        ha.YTick = [1 10 100];
        ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
        ylabel('Counts', 'Units','normalized','Position',[-0.28 0.5])
    case 'probability'
        if strcmp(hist_Yscale_option,'log')
            ylim([0.0005 0.5])
            ha.YTick = [0.001 0.01 0.1];
            ha.YTickLabel = {'10^{-3}';'10^{-2}';'10^{-1}'};
        else
            ylim([0 0.25])
            ha.YTick = [0 0.25];
        end
        ylabel('Fraction of fields', 'Units','normalized','Position',[-0.23 0.5])
end
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

%% panel D Inset
axes(panel_D_inset);
cla
hold on

h = histogram(fields_size_wild);
h.FaceColor = 0.5*[1 1 1];
h.BinWidth = 1;
h.Normalization = hist_normalization_option;
h = histogram(fields_size_lab);
h.FaceColor = color_Lab;
h.FaceAlpha = Lab_alpha;
h.BinWidth = 1;
h.Normalization = hist_normalization_option;
ha=gca;
ha.YScale = 'linear';
ha = gca;
ylim([0 0.3])
ha.XLim = [0 20];
ha.XTick = [0 10 20];
ha.YTick = [0 0.15 0.3];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

%% panel G - field size bars
axes(panel_G);
cla
hold on
text(-0.45,1.15, 'G', 'Units','normalized','FontWeight','bold');

% % bar plot version:
% x = [1,2];
% y = [mean(fields_size_lab),mean(fields_size_wild)];
% hb = bar(x,y);
% hb.FaceColor = 'flat';
% hb.CData = [color_Lab; 0.5*[1 1 1]];
% sem_lab = std(fields_size_lab)/sqrt(length(fields_size_lab));
% sem_wild = std(fields_size_wild)/sqrt(length(fields_size_wild));
% 
% he=errorbar(x,y,[sem_lab sem_wild]);
% he.CapSize = 2;
% he.LineStyle = 'none';
% he.Color = 'k';
% h=gca;
% h.XTick = [1 2];
% h.XTickLabels = {'Lab','Wild'};
% h.TickLength = [0.03 0.03];
% h.YLim = [0,7];

% boxplot version:
plot_boxplot(fields_size_lab,fields_size_wild,gca,[0 0.75 0],[0.5 0.5 0.5])
h = gca;
h.XTickLabels = {'Lab','Wild'};
h.TickLength = [0.03 0.03];
h.YLim = [0,10];

ylabel('Field size (m)', 'Units','normalized','Position',[-0.22 0.5])

% [~,P_KS] = kstest2(fields_size_lab, fields_size_wild);
% text(1.7,0.8, sprintf('P_{KS} = %0.2e',P_KS),'Units','normalized','FontSize',7,'HorizontalAlignment','right');
% text(1,0.8, sprintf('N = %d',length(fields_size_lab)),'Units','normalized','FontSize',7,'HorizontalAlignment','right','color',color_Lab);
% text(1,0.7, sprintf('N = %d',length(fields_size_wild)),'Units','normalized','FontSize',7,'HorizontalAlignment','right','color',0.5*[1 1 1]);
% [~,P_ttest] = ttest2(fields_size_lab, fields_size_wild);
% text(1.2,1, sprintf('P_{T} = %0.00e',P_ttest),'Units','normalized','FontSize',7,'HorizontalAlignment','right');
[P_wilcoxon,~,wilcoxon_stats] = ranksum(fields_size_lab, fields_size_wild);
text(1.5,0.85, sprintf('P = %0.00e',P_wilcoxon),'Units','normalized','FontSize',7,'HorizontalAlignment','right');
% [~,P_ttest_log] = ttest2(log(fields_size_lab), log(fields_size_wild));
% text(1.7,0.6, sprintf('P_{T(log)} = %0.00e',P_ttest_log),'Units','normalized','FontSize',7,'HorizontalAlignment','right');

h.XRuler.TickLabelGapMultiplier = -0.2;
h.YRuler.TickLabelGapMultiplier = -0.05;

%% panel E - field ratio (largest/smallest)
axes(panel_E);
cla
hold on
text(-0.55,1.25, 'E', 'Units','normalized','FontWeight','bold');
LS_field_ratio_all_lab = nan(1,length(cells_lab));
LS_field_ratio_dir = nan(2,length(cells_lab));
for ii_cell = 1:length(cells_lab)
    cell = cells_lab(ii_cell);
    % pooled stats - check at least one direction is signif
    if any([cell.signif.TF])
        LS_field_ratio_all_lab(ii_cell) = cell.stats.all.field_ratio_LS;
    end
    % per dir stats - check signif per direction
    for ii_dir = 1:2
        if cell.signif(ii_dir).TF 
            LS_field_ratio_dir(ii_dir,ii_cell) = cell.stats.dir(ii_dir).field_ratio_LS;
        end
    end
end
LS_field_ratio_all_lab(isnan(LS_field_ratio_all_lab)) = [];

LS_field_ratio_all_wild = nan(1,length(cells_wild));
LS_field_ratio_dir = nan(2,length(cells_wild));
for ii_cell = 1:length(cells_wild)
    cell = cells_wild(ii_cell);
    % pooled stats - check at least one direction is signif
    if any([cell.signif.TF])
        LS_field_ratio_all_wild(ii_cell) = cell.stats.all.field_ratio_LS;
    end
    % per dir stats - check signif per direction
    for ii_dir = 1:2
        if cell.signif(ii_dir).TF 
            LS_field_ratio_dir(ii_dir,ii_cell) = cell.stats.dir(ii_dir).field_ratio_LS;
        end
    end
end
LS_field_ratio_all_wild(isnan(LS_field_ratio_all_wild)) = [];

nBinEdges = 9;
edges = logspace(0,log10(25),nBinEdges);
clear h

h(1) = histogram(LS_field_ratio_all_wild);
h(1).BinEdges = edges;
h(1).FaceColor = 0.5*[1 1 1];
h(1).Normalization = hist_normalization_option;

h(2) = histogram(LS_field_ratio_all_lab);
h(2).BinEdges = edges;
h(2).FaceColor = color_Lab;
h(2).FaceAlpha = Lab_alpha;
h(2).Normalization = hist_normalization_option;
ha=gca;
ha.YScale = hist_Yscale_option;
ha.XScale = 'log';
ha.XLim = [0 27];
% % ha.YLim = [0 10];
% % ha.YLim = [7e-1 120];
% ha.YLim = [0.7 max(h(1).Values)*1.1];
% % ha.XTick = [0:5:35];
ha.XTick = [1 2 5 10 20];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
switch hist_normalization_option
    case 'count'
        ha.YLim = [0.7 max(h(1).Values)*1.1];
        ha.YTick = [1 10 100];
        ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
        ylabel('No. of cells', 'Units','normalized','Position',[-0.24 0.5])
    case 'probability'
        if strcmp(hist_Yscale_option,'log')
            ylim([0.002 0.5])
            ha.YTick = [0.01 0.1];
            ha.YTickLabel = {'10^{-2}';'10^{-1}'};
        else
            ylim([0 0.3])
            ha.YTick = [0 0.3];
        end
        ylabel('Fraction of cells', 'Units','normalized','Position',[-0.33 0.5])
end
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.17]);

%% panel E Inset
axes(panel_E_inset);
cla
hold on
clear h

h(1) = histogram(LS_field_ratio_all_wild);
h(1).BinEdges = edges;
h(1).FaceColor = 0.5*[1 1 1];
h(1).Normalization = hist_normalization_option;

h(2) = histogram(LS_field_ratio_all_lab);
h(2).BinEdges = edges;
h(2).FaceColor = color_Lab;
h(2).FaceAlpha = Lab_alpha;
h(2).Normalization = hist_normalization_option;
ha=gca;
ha.YScale = 'linear';
ha.XScale = 'log';
ha.XLim = [0 27];
ha.XTick = [1 2 5 20];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ylim([0 0.3])
ha.YTick = [0 0.15 0.3];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.001;

%% panel H - field ratio bars (largest/smallest)
axes(panel_H);
cla
hold on
text(-0.55,1.15, 'H', 'Units','normalized','FontWeight','bold');

LS_field_ratio_lab_noNans = LS_field_ratio_all_lab(~isnan(LS_field_ratio_all_lab));
LS_field_ratio_wild_noNans = LS_field_ratio_all_wild(~isnan(LS_field_ratio_all_wild));

% % % bar version:
% % x = [1,2];
% % y = [mean(LS_field_ratio_lab_noNans),mean(LS_field_ratio_wild_noNans)];
% % hb = bar(x,y);
% % hb.FaceColor = 'flat';
% % hb.CData = [color_Lab; 0.5*[1 1 1]];
% % sem_lab = std(LS_field_ratio_lab_noNans)/sqrt(length(LS_field_ratio_lab_noNans));
% % sem_wild = std(LS_field_ratio_wild_noNans)/sqrt(length(LS_field_ratio_wild_noNans));
% % 
% % he=errorbar(x,y,[sem_lab sem_wild]);
% % he.CapSize = 2;
% % he.LineStyle = 'none';
% % he.Color = 'k';
% % h=gca;
% % h.XTick = [1 2];
% % h.XTickLabels = {'Lab','Wild'};
% % h.TickLength = [0.03 0.03];
% % h.YLim = [0,6];

% boxplot version:
plot_boxplot(LS_field_ratio_lab_noNans,LS_field_ratio_wild_noNans,gca,[0 0.75 0],[0.5 0.5 0.5])
h = gca;
h.XTickLabels = {'Lab','Wild'};
h.TickLength = [0.03 0.03];

ylabel({'Field size ratio';'largest/smallest'}, 'Units','normalized','Position',[-0.17 0.5])

% [~,P_KS] = kstest2(LS_field_ratio_all_lab, LS_field_ratio_all_wild);
% text(1.5,0.8, sprintf('P_{KS} = %.02f',P_KS),'Units','normalized','FontSize',7,'HorizontalAlignment','right');
% LS_field_ratio_lab_noNans = LS_field_ratio_all_lab(~isnan(LS_field_ratio_all_lab));
% LS_field_ratio_wild_noNans = LS_field_ratio_all_wild(~isnan(LS_field_ratio_all_wild));
% text(1,0.8, sprintf('N = %d',length(LS_field_ratio_lab_noNans)),'Units','normalized','FontSize',7,'HorizontalAlignment','right','color',color_Lab);
% text(1,0.7, sprintf('N = %d',length(LS_field_ratio_wild_noNans)),'Units','normalized','FontSize',7,'HorizontalAlignment','right','color',0.5*[1 1 1]);
% [~,P_ttest] = ttest2(LS_field_ratio_all_lab, LS_field_ratio_all_wild);
% text(1.2,1, sprintf('P_{T} = %.02f',P_ttest),'Units','normalized','FontSize',7,'HorizontalAlignment','right');
[P_wilcoxon,~,wilcoxon_stats] = ranksum(LS_field_ratio_all_lab, LS_field_ratio_all_wild);
text(1.4,0.85, sprintf('P = %.02f',P_wilcoxon),'Units','normalized','FontSize',7,'HorizontalAlignment','right');

h.XRuler.TickLabelGapMultiplier = -0.2;
h.YRuler.TickLabelGapMultiplier = -0.05;

%% panel I - plot TT positions (only with signif cells!)
axes(panel_I);
cla
hold on
text(-0.3,1.18, 'I', 'Units','normalized','FontWeight','bold');

pop_signif_wild = cat(1,cells_wild.signif);
pop_signif_TF_wild = arrayfun(@(x)(x.TF), pop_signif_wild);
signif_wild_at_least_one_dir_IX = any(pop_signif_TF_wild,2);
pop_details_wild = cat(1,cells_wild.details);
pop_details_wild(~signif_wild_at_least_one_dir_IX) = [];
[TT_ID_wild, IX] = unique({pop_details_wild(:).TT_ID});
TT_pos_PD_wild = [pop_details_wild(IX).TT_pos_proximodistal_prc];
TT_pos_Long_wild = [pop_details_wild(IX).TT_pos_longitudinal_prc];
scatter(TT_pos_PD_wild*100, TT_pos_Long_wild*100, 5, [0.5 0.5 0.5]);

pop_signif_lab = cat(1,cells_lab.signif);
pop_signif_TF_lab = arrayfun(@(x)(x.TF), pop_signif_lab);
signif_lab_at_least_one_dir_IX = any(pop_signif_TF_lab,2);
pop_details_lab = cat(1,cells_lab.details);
pop_details_lab(~signif_lab_at_least_one_dir_IX) = [];
[TT_ID_lab, IX] = unique({pop_details_lab(:).TT_name});
TT_pos_PD_lab = [pop_details_lab(IX).TT_pos_proximodistal_prc];
TT_pos_Long_lab = [pop_details_lab(IX).TT_pos_longitudinal_prc];
% move 2 TTs which are exactely on top of a wild TT:
TT_pos_PD_lab(1) = TT_pos_PD_lab(1)-0.01;
TT_pos_Long_lab(1) = TT_pos_Long_lab(1)+0.001;
TT_pos_PD_lab(4) = TT_pos_PD_lab(4)-0.01;
TT_pos_Long_lab(4) = TT_pos_Long_lab(4)+0.001;

scatter(TT_pos_PD_lab*100, TT_pos_Long_lab*100, 5, color_Lab);

xlim([0 100])
ylim([16 22]);
set(gca,'ytick',16:2:22);
xlabel('Proximo-distal axis (%)');
ylabel('Longitudinal axis (%)');

%% legend panel
axes(panel_I_legend);
cla
hold on
prm = PARAMS_GetAll();
scatter(0,2,10,[0.5 0.5 0.5]);
scatter(0,1,10,color_Lab);
text(0+0.2, 2, 'Wild', 'FontSize',7,'HorizontalAlignment','Left');
text(0+0.2, 1, 'Lab', 'FontSize',7,'HorizontalAlignment','Left');
xlim([0 1]);
ylim([0 3]);
set(gca,'Visible','off');



%% print/save the figure
% fig_name_out = fullfile(res_dir, sprintf('%s__corr_%s_%d_paramset_%d',fig_name_str,corr_type,field_speed_opt,prm.parmaset));
% fig_name_out = fullfile(res_dir, sprintf('%s__example%d',fig_name_str,examples_option));
% if match_TT_pos
%     fig_name_out = fullfile(res_dir, sprintf('%s_revised_nan_corrected',[fig_name_str '_match_TT_pos']));
% else
%     fig_name_out = fullfile(res_dir, sprintf('%s_revised_',fig_name_str));
% end
fig_name_out = fullfile(res_dir, fig_name_str);

print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

