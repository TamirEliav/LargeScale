%% Large Scale - Fig. S6 - many FR maps examples

%%
clear 
clc

%% params
grp = 0;

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
n_cols = 3;
n_rows = 6;
n_examples = n_cols * n_rows;
panel_size = [5 1.3];
x_positions = linspace(0,12,n_cols)+2;
y_positions = linspace(0,9,n_rows);
panels_AB_y_offset = [14.5 2.5];
clear panels 
for ii_dataset = 1:2
   for r = 1:n_rows
       for c = 1:n_cols
            x = x_positions(c);
            y = y_positions(r);
            y = y + panels_AB_y_offset(ii_dataset);
            panels(ii_dataset,r,c) = axes('position', [x y panel_size]);
       end
    end
end
panels = panels(:,end:-1:1,:)

%% load data
cells_wild = load('L:\processed_data_structs\cells_bat_200m.mat');
cells_lab = load('L:\processed_data_structs\cells_lab.mat');
cells_wild = cells_wild.cells;
cells_lab = cells_lab.cells;
data = {cells_wild;cells_lab};

%% arrange data
cells_all = {};
for ii_dataset = 1:2
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
%         cells(ii_cell).FE = cells(ii_cell).FE{dir};
        cells(ii_cell).largest = max([cells(ii_cell).fields.width_prc]);
    end
%     stats=[cells.stats];
%     stats_dir = [stats.dir];
    [~,sort_IX] = sort([cells.largest]);
    cells = cells(sort_IX);
%     hold on
%     histogram([cells.largest],'Normalization','pdf','NumBins',25);
%     plot([cells.largest],'o-');
    IDs = 1:length(cells);
    n_grps = ceil(length(cells) / n_examples);
    grps = mod(IDs,n_grps);
    [cells.ID] = disperse(IDs);
    [cells.grp] = disperse(grps);

    cells_all{ii_dataset} = cells;
    clear cells1 cells2 cells12 signif
end

%% plot 
x_ticks = [0:50:200];
for ii_dataset = 1:2
    cells = cells_all{ii_dataset};
    cells_IX = find([cells.grp] == grp);
    for ii_panel = 1:length(cells_IX)
        ii_cell = cells_IX(ii_panel);
        %% choose axis
        axes(panels(ii_dataset,ii_panel));
        cla
        hold on
        %% arrange data
        cell = cells(ii_cell);
        %% plot
        plot(cell.FR_map.all.bin_centers,cell.FR_map.all.PSTH,'k','LineWidth',1.2);
        text(0.01,1.1,sprintf('%.2fm',cell.largest),'FontSize',6,'Units','normalized','HorizontalAlignment','left');
        text(0.99,1.1,sprintf('%d / %d (%.0f%%)',cell.ID,length(cells),100*cell.ID/length(cells)),'FontSize',6,'Units','normalized','HorizontalAlignment','right');
%         title(ii_panel)
        hax=gca;
        hax.XTick = [x_ticks];
        hax.XTickLabel = [];
        hax.XAxis.TickLength(1) = 0.02;
%         hax.TickLength = [0.01 0];
    end
end

% labels
axes(panels(1,1));
text(-0.2,1.2, 'A', 'Units','normalized','FontWeight','bold');
axes(panels(2,1));
text(-0.2,1.2, 'B', 'Units','normalized','FontWeight','bold');

for ii_dataset = 1:size(panels,1)
    for c = 1:size(panels,3)
        axes(panels(ii_dataset,end,c));
        hax=gca;
        hax.XTick = x_ticks;
        hax.XTickLabel = x_ticks;
        hax.XRuler.TickLabelGapOffset = -1;
        xlabel('Position (m)', 'Units','normalized','Position',[0.5 -0.35]);
    end
end

axes(panels(1,4,1));
ylabel('Firing rate (Hz)', 'Units','normalized','Position',[-0.12 1.3]);
axes(panels(2,4,1));
ylabel('Firing rate (Hz)', 'Units','normalized','Position',[-0.12 1.3]);

% for ii_dataset = 1:size(panels,1)
%     for r = 1:size(panels,2)
%         axes(panels(ii_dataset,r,1));
%         hax=gca;
%         hax.YRuler.TickLabelGapOffset = -1;
%         ylabel({'Firing rate';'(Hz)'},   'Units','normalized','Position',[-0.07 0.42]);
%     end
% end

%%





%% print/save the figure
fig_name_out = fullfile(res_dir, sprintf('%s_grp_%d',fig_name_str,grp));
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

