%% Large Scale - Fig. S6 - many FR maps examples

%%
% clear 
% clc
close all
clearvars -except grp data

%% params

% order_feature = 'largest';
order_feature = 'mean';
% order_feature = 'median';

% grp_type = 'mod';
grp_type = 'div';

grp = 1;

lw = 0.8;

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S6';
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
switch 2
    case 2
        n_cols = 2;
        n_rows = 9;
        panel_size = [7.5 1.3];
        x_positions = linspace(0,9,n_cols)+2.7;
        y_positions = linspace(0,16,n_rows);
        panels_AB_y_offset = [8 2.5];
    case 3
        n_cols = 3;
        n_rows = 6;
        panel_size = [5 1.3];
        x_positions = linspace(0,12,n_cols)+2;
        y_positions = linspace(0,9,n_rows);
        panels_AB_y_offset = [14.5 2.5];
end
y_positions = flip(y_positions);
n_examples = n_cols * n_rows;
n_examples = n_examples / 2; % add resters
clear panels 
for ii_dataset = 1
   for r = 1:n_rows
       for c = 1:n_cols
            x = x_positions(c);
            y = y_positions(r);
            y = y + panels_AB_y_offset(ii_dataset);
            panels(ii_dataset,r,c) = axes('position', [x y panel_size]);
            title(sprintf('r=%d c=%d',r,c))
       end
    end
end

% panels = panels(:,end:-1:1,:);
% panels = permute(panels, [1 3 2]);
% panels = reshape(panels, 1,[])
% for ii = 1:length(panels)
%     axes(panels(1,ii))
%     title(ii)
% end

%% load data
if ~exist('data','var')
%     cells_wild = load('L:\processed_data_structs\cells_bat_200m.mat');
    cells_wild = load('L:\processed_data_structs\cells_with_FE.mat')
    % cells_lab = load('L:\processed_data_structs\cells_lab.mat');
    cells_wild = cells_wild.cells;
    % cells_lab = cells_lab.cells;
    % data = {cells_wild;cells_lab};
    data = {cells_wild};
end

%% arrange data
cells_all = {};
for ii_dataset = 1:length(data)
    cells = data{ii_dataset};
    cells1 = cells;
    cells2 = cells;
    [cells1.dir] = disperse(repelem(1,length(cells)));
    [cells2.dir] = disperse(repelem(2,length(cells)));
    cells12=[cells1;cells2]';
    signif = cat(1,cells.signif);
    signif = arrayfun(@(x)(x.TF),signif);
    cells = cells12(signif);
    for ii_cell = 1:length(cells)
        dir = cells(ii_cell).dir;
        cells(ii_cell).stats.dir = cells(ii_cell).stats.dir(dir);
        cells(ii_cell).signif = cells(ii_cell).signif(dir);
        cells(ii_cell).FR_map = cells(ii_cell).FR_map(dir);
        cells(ii_cell).fields = cells(ii_cell).fields{dir};
        cells(ii_cell).FE = cells(ii_cell).FE{dir};
        fields_valid_speed = [cells(ii_cell).fields];
        fields_valid_speed(fields_valid_speed.in_low_speed_area) = [];
        cells(ii_cell).largest = max([fields_valid_speed.width_prc]);
        cells(ii_cell).smallest = min([fields_valid_speed.width_prc]);
        cells(ii_cell).mean = mean([fields_valid_speed.width_prc]);
        cells(ii_cell).median = median([fields_valid_speed.width_prc]);
    end
    switch order_feature
        case 'largest'
            [~,sort_IX] = sort([cells.largest]);
        case 'mean'
            [~,sort_IX] = sort([cells.mean]);
        case 'median'
            [~,sort_IX] = sort([cells.median]);
    end
    cells = cells(sort_IX);
    IDs = 1:length(cells);
    n_grps = ceil(length(cells) / n_examples);
    switch grp_type
        case 'mod'
            grps = mod(IDs,n_grps);
        case 'div'
            grps = ceil(IDs./n_examples);
    end
    [cells.ID] = disperse(IDs);
    [cells.grp] = disperse(grps);

    cells_all{ii_dataset} = cells;
    clear cells1 cells2 cells12 signif
end

%% choose 
switch grp
    case -1
        % grouped by largrst field
        % final choice
        % order by largest field
        wild_cells_IX = [327 305 301 272 256 243 210 193 182 165 143 118 103 82 66 56 25 19]; % sorted index 
        [cells_all{ii_dataset}(wild_cells_IX).grp] = disperse(repelem(-1,length(wild_cells_IX)));
    case -2
end

%% make sure none of the chosen examples are from fig 2
fig_2_cell_examples = [433;  56;  51; 609; 419; 477;  57; 628; 337];
if grp<0
    fig_supp_cell_examples = arrayfun(@(x)(x.details.cell_num), cells_all{ii_dataset}(wild_cells_IX));
    is_duplicate = ismember(fig_supp_cell_examples,fig_2_cell_examples );
    if any(is_duplicate)
        duplicate_cell_num = fig_supp_cell_examples(is_duplicate)
        error('Same neuron as in Fig 2!!!!');
    end
    fig_supp_cell_examples
end

%% plot 
x_ticks = [0:50:200];
for ii_dataset = 1:length(data)
    cells = cells_all{ii_dataset};
    cells_IX = find([cells.grp] == grp);
%     cells_IX = flip(cells_IX);
    for ii_panel = 1:length(cells_IX)
        ii_cell = cells_IX(ii_panel);
        %% choose axis
        axes(panels(ii_dataset,ii_panel));
        cla
        hold on
        %% arrange data
        cell = cells(ii_cell);
        %% plot
        plot(cell.FR_map.all.bin_centers,cell.FR_map.all.PSTH,'k','LineWidth',lw);
        text(0,1.2,sprintf('cell %d    (%s)',cell.details.cell_num,cell.details.cell_ID),'FontSize',6,'Units','normalized','HorizontalAlignment','left','Interpreter','none');
%         text(.5,1.0,strrep(sprintf('%s',cell.details.cell_ID),'_',' '),'FontSize',6,'Units','normalized','HorizontalAlignment','center');
        text(1,1.30,sprintf('max/min/mean/median/ratio'),'FontSize',6,'Units','normalized','HorizontalAlignment','Right');
        text(1,1.15,sprintf('%.1f / %.1f / %.1f / %.1f / %.1f',cell.largest,cell.smallest,cell.mean,cell.median, cell.stats.all.field_ratio_LS),'FontSize',6,'Units','normalized','HorizontalAlignment','Right');
        text(1,1.0,sprintf('%d/%d=%.0f%%',cell.ID,length(cells),100*cell.ID/length(cells)),'FontSize',6,'Units','normalized','HorizontalAlignment','Right');
%         text(0,1.1,sprintf('%.1fm (%.0f%%)',cell.largest,100*cell.ID/length(cells)),'FontSize',6,'Units','normalized','HorizontalAlignment','Left');
%         text(0.01,1.1,sprintf('%.1f',cell.stats.dir.field_ratio_LS),'FontSize',6,'Units','normalized','HorizontalAlignment','left');
%         text(0.01,1.1,sprintf('%.2fm',cell.largest),'FontSize',6,'Units','normalized','HorizontalAlignment','left');
%         text(0.99,1.1,sprintf('%d / %d (%.0f%%)',cell.ID,length(cells),100*cell.ID/length(cells)),'FontSize',6,'Units','normalized','HorizontalAlignment','right');
%         title(ii_panel)
        hax=gca;
        hax.XTick = [x_ticks];
        hax.XTickLabel = [];
        hax.XAxis.TickLength(1) = 0.02;
%         hax.TickLength = [0.01 0];
    end
end

%% labels
% axes(panels(1,1));
% text(-0.2,1.2, 'A', 'Units','normalized','FontWeight','bold');
% axes(panels(2,1));
% text(-0.2,1.2, 'B', 'Units','normalized','FontWeight','bold');

for ii_dataset = 1:size(panels,1)
    for c = 1:size(panels,2)
        axes(panels(ii_dataset,c,end));
        hax=gca;
        hax.XTick = x_ticks;
        hax.XTickLabel = x_ticks;
        hax.XRuler.TickLabelGapOffset = -1;
        xlabel('Position (m)', 'Units','normalized','Position',[0.5 -0.35]);
    end
end

switch n_cols
    case 2
        axes(panels(1,1,5));
        ylabel('Firing rate (Hz)', 'Units','normalized','Position',[-0.10 0.5]);
        % axes(panels(2,1,5));
        % ylabel('Firing rate (Hz)', 'Units','normalized','Position',[-0.12 0.5]);
    case 3
        axes(panels(1,1,4));
        ylabel('Firing rate (Hz)', 'Units','normalized','Position',[-0.12 1.3]);
        % axes(panels(2,1,4));
        % ylabel('Firing rate (Hz)', 'Units','normalized','Position',[-0.12 1.3]);
end
% for ii_dataset = 1:size(panels,1)
%     for r = 1:size(panels,2)
%         axes(panels(ii_dataset,r,1));
%         hax=gca;
%         hax.YRuler.TickLabelGapOffset = -1;
%         ylabel({'Firing rate';'(Hz)'},   'Units','normalized','Position',[-0.07 0.42]);
%     end
% end

%%
fig_details_str = {
    sprintf('order_feature: %s',order_feature)
    sprintf('grp_type: %s',grp_type)
    sprintf('grp: %d',grp)
    };
annotation('textbox', [0.2 0.25 0 0], 'String',fig_details_str, 'HorizontalAlignment','left','Interpreter','none', 'FitBoxToText','on');

fig_S6_cells_by_largest = [405 323 176 178 67 113 684 298 462 478 486 723 578 67 658 635 645 623];
fig_2_cell_examples_str = ...
    {'fig 2 cells:';...
    sprintf('%d, ',fig_2_cell_examples);...
    '';...
    'fig S6 chosen by largest:';...
    sprintf('%d, ',fig_S6_cells_by_largest(1:9));...
    sprintf('%d, ',fig_S6_cells_by_largest(10:18))};
annotation('textbox', [0.55 0.25 0 0], 'String',fig_2_cell_examples_str, 'HorizontalAlignment','left','Interpreter','none', 'FitBoxToText','on');

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

