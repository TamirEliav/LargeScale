%% Large Scale - Fig. S6 - many FR maps examples

%%
% clear 
% clc
close all
clearvars -except grp data

%% params

% order_feature = 'largest';
% order_feature = 'mean';
% order_feature = 'median';
order_feature = 'ratio';

% grp_type = 'mod';
grp_type = 'div';

% grp = -1;

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S6b';
fig_caption_str = 'many FR map examples';
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
disp([fig_name_str ': ' fig_caption_str]);
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
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');
pause(0.2); % workaround to solve matlab automatically changing the axes positions...

% create panels
% create panels
panel_size_raster = [5 1];
panel_size_FR_map = [5 1];
panel_pos = [2 12];
clear panels
for ii=1:3
    for jj=1:3
        offset_x = (ii-1)*6;
        offset_y = (jj-1)*4.3;
        offset = panel_pos + [offset_x offset_y];
        panels(ii,jj,1) = axes('position', [offset+[0 2.25] panel_size_FR_map]);
        panels(ii,jj,2) = axes('position', [offset+[0 1   ] panel_size_raster]);
        panels(ii,jj,3) = axes('position', [offset+[0 0   ] panel_size_raster]);
    end
end
panels = panels(:,3:-1:1,:);
panels = reshape(panels,[9 3]);
n_examples = 9;

%% load data
if ~exist('data','var')
    data = load('L:\processed_data_structs\cells_with_FE.mat');
end

%% arrange data
cells = data.cells;
signif = cat(1,cells.signif);
signif = arrayfun(@(x)(x.TF),signif);
signif = any(signif,2);
cells = cells(signif);

for ii_cell = 1:length(cells)
    cells(ii_cell).cell_num = cells(ii_cell).details.cell_num;
    cells(ii_cell).largest = cells(ii_cell).stats.all.field_largest;
    cells(ii_cell).smallest = cells(ii_cell).stats.all.field_smallest;
    cells(ii_cell).mean = cells(ii_cell).stats.all.field_size_mean;
    cells(ii_cell).median = cells(ii_cell).stats.all.field_size_mean;
    cells(ii_cell).ratio = cells(ii_cell).stats.all.field_ratio_LS;
    cells(ii_cell).nfields = cells(ii_cell).stats.all.field_num;
    if cells(ii_cell).nfields == 1
        cells(ii_cell).largest = cells(ii_cell).mean;
        cells(ii_cell).smallest = cells(ii_cell).mean;
        cells(ii_cell).ratio = 1;
    end
end

switch order_feature
    case 'mean'
        order_feature_values = [cells.mean];
    case 'median'
        order_feature_values = [cells.median];
    case 'largest'
        order_feature_values = [cells.largest];
    case 'ratio'
        order_feature_values = [cells.ratio];
end
% sort
[order_feature_values,sort_IX] = sort(order_feature_values);
cells = cells(sort_IX);
order_feature_IX = 1:length(cells);
n_grps = ceil(length(cells) / n_examples);
switch grp_type
    case 'mod'
        grps = mod(order_feature_IX,n_grps);
    case 'div'
        grps = ceil(order_feature_IX./n_examples);
end
[cells.grp] = disperse(grps);
[cells.order_feature_IX] = disperse(order_feature_IX);
[cells.order_feature_value] = disperse(order_feature_values);
[cells.order_feature_prc] = disperse(100 * order_feature_IX / length(cells));

%% choose 
Tamir_options = [478 52 85 527 578 514 623 684 66 645 298 178 113 486 447 285 462 274 67 231 694 582 72 682 631 658 188 176 474 323];
Tamir_options = sort(Tamir_options);
switch grp
    case -1
        IX = ismember([cells.cell_num], Tamir_options(1:9) );
        [cells(IX).grp] = disperse(repelem(grp,sum(IX)));
    case -2
        IX = ismember([cells.cell_num], Tamir_options(10:18) );
        [cells(IX).grp] = disperse(repelem(grp,sum(IX)));
    case -3
        IX = ismember([cells.cell_num], Tamir_options(19:27) );
        [cells(IX).grp] = disperse(repelem(grp,sum(IX)));
    case -4
        IX = ismember([cells.cell_num], Tamir_options(28:30) );
        [cells(IX).grp] = disperse(repelem(grp,sum(IX)));
end

cell_examples = cells([cells.grp] == grp);

%% make sure none of the chosen examples are from fig 2
fig_2_cell_examples = [433;  56;  51; 609; 419; 477;  57; 628; 337];
if grp<0
    fig_supp_cell_examples = [cell_examples.cell_num];
    is_duplicate = ismember(fig_supp_cell_examples, fig_2_cell_examples);
    if any(is_duplicate)
        duplicate_cell_num = fig_supp_cell_examples(is_duplicate)
        error('Same neuron as in Fig 2!!!!');
    end
    fig_supp_cell_examples
end

%% plot
prm = PARAMS_GetAll();
for ii_cell = 1:length(cell_examples)
    %%
    cell = cell_examples(ii_cell);
    c = prm.graphics.colors.flight_directions;
    
    % map
    axes(panels(ii_cell, 1));
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
    
    % fields (same axis)
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
    text(0,1.5,sprintf('cell %d',cell.details.cell_num),'FontSize',6,'Units','normalized','HorizontalAlignment','left');
    text(0,1.2,cell.details.cell_ID,'FontSize',5,'Units','normalized','HorizontalAlignment','left','Interpreter','none','Interpreter','none');
    text(1,1.5,sprintf('max/min/mean/median/ratio'),'FontSize',6,'Units','normalized','HorizontalAlignment','Right');
    text(1,1.2,sprintf('%.1f / %.1f / %.1f / %.1f / %.1f',...
                        cell.largest,...
                        cell.smallest,...
                        cell.mean,...
                        cell.median,...
                        cell.ratio),'FontSize',6,'Units','normalized','HorizontalAlignment','Right');
    text(1,1.0,sprintf('%d/%d=%.0f%%',cell.order_feature_IX, length(cells), cell.order_feature_prc),'FontSize',6,'Units','normalized','HorizontalAlignment','Right');
    
    % rasters
    FEs = [cell.FE];
    for ii_dir=1:2
        axes(panels(ii_cell, ii_dir+1));
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
    axes(panels(ii, 3));
    xlabel('Position (m)', 'Units','normalized','Position',[0.5 -0.35]);
end
for ii = [1 4 7]
    axes(panels(ii, 3));
%     ylabel('Time (min)',   'Units','normalized','Position',[-0.1 1]);
    ylabel('Flight no.',   'Units','normalized','Position',[-0.1 1]);
    axes(panels(ii, 1));
    ylabel({'Firing rate';'(Hz)'},   'Units','normalized','Position',[-0.07 0.42]);
end

%% add direction arrows
arrow_x = 0.1 +[0 0.05];
arrow_y = repelem(0.93,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.008,  'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y      ,  'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);

%%













%%
fig_details_str = {
    sprintf('order_feature: %s',order_feature)
    sprintf('grp_type: %s',grp_type)
    sprintf('grp: %d',grp)
    };
annotation('textbox', [0.2 0.4 0 0], 'String',fig_details_str, 'HorizontalAlignment','left','Interpreter','none', 'FitBoxToText','on');

fig_S6_cells_by_largest = [405 323 176 178 67 113 684 298 462 478 486 723 578 67 658 635 645 623];
fig_2_cell_examples_str = ...
    {'fig 2 cells:';...
    sprintf('%d, ',fig_2_cell_examples);...
    '';...
    'fig S6 chosen by largest:';...
    sprintf('%d, ',fig_S6_cells_by_largest(1:9));...
    sprintf('%d, ',fig_S6_cells_by_largest(10:18))};
annotation('textbox', [0.55 0.4 0 0], 'String',fig_2_cell_examples_str, 'HorizontalAlignment','left','Interpreter','none', 'FitBoxToText','on');

%% print/save the figure
fig_name_out = fullfile(res_dir, sprintf('%s_order_by_%s_grp_type_%s_grp_%d',...
    fig_name_str,...
    order_feature,...
    grp_type,...
    grp));
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

