%% Large Scale - Fig. SXXX - compare main results between bats

%%
clear 
clc

%% plotting options
% plot_style = 'errorbars';
% plot_style = 'boxplot';
plot_style = 'violin';

plot_yscale = 'linear';
% plot_yscale = 'log';

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'Fig_SXXX_results_per_bat';
fig_caption_str = 'compare main results between bats';
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
panels_size = [12 4.5];
panels_y_pos = flip(linspace(5,18,3));
panel_A = axes('position', [ 3 panels_y_pos(1) panels_size]);
panel_B = axes('position', [ 3 panels_y_pos(2) panels_size]);
panel_C = axes('position', [ 3 panels_y_pos(3) panels_size]);

%% load data
% first, decide which cells to load data for
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
clear cells stats cells_details cells_t
% now load the relevant data for these cells
cells = cellfun(@(c)(cell_load_data(c,'details','stats','signif','fields')), cells_ID, 'UniformOutput',0);
cells = [cells{:}];
cells_details = [cells.details];
cells_bat_num = [cells_details.bat];
bats_colors = arrayfun(@(bat)(prm.graphics.colors.bats(bat)), bats, 'UniformOutput',0);
whos cells

%% 1) panel A - No. of fields
axes(panel_A);
cla
hold on
text(-0.1,1.1, 'A', 'Units','normalized','FontWeight','bold');

% arrange data
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
end
values = nFields(:);
grps = repmat(cells_bat_num,2,1);
grps = grps(:);

% plot
switch plot_style 
    case 'errorbars'
        x = 1:length(bats);
        y = accumarray(grps, values, [], @nanmean, nan);
        y = y(bats);
        err = accumarray(grps, values, [], @nansem, nan);
        err = err(bats);
        h=bar(x,y);
        h.FaceColor = 0.5*[1 1 1];
        h=errorbar(x,y,err);
        h.LineStyle='none';
        h.Color = 'k';
        set(gca,'XTick',1:length(bats), 'XTickLabel', bats);
    case 'boxplot'
        boxplot(values, grps);
    case 'violin'
        rng(0);
        h=violinplot(values, grps);
        [h.ViolinColor] = disperse( bats_colors );
end

% set axis properties
ha = gca;
% ha.YScale = plot_yscale;
ha.TickDir='out';
ha.TickLength = [0.01 0.01];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
ylabel('No. of fields per direction', 'Units','normalized','Position',[-0.05 0.5]);
xlabel('Bat no.', 'Units','normalized','Position',[0.5 -0.13]);

%% 2) panel B - Field size
axes(panel_B);
cla
hold on
text(-0.1,1.1, 'B', 'Units','normalized','FontWeight','bold');

% arrange data
fields_size = [];
fields_bat_num = [];
for ii_dir = 1:2
    for ii_cell = 1:length(cells)
        cell = cells(ii_cell);
        if ~cell.signif(ii_dir).TF % check signif per direction
            continue;
        end
        fields = cell.fields{ii_dir};
        fields([fields.in_low_speed_area]) = []; % remove fields in low speed area
        fields_size = [fields_size fields.width_prc];
        fields_bat_num = [fields_bat_num repelem(cell.details.bat,1,length(fields))];
    end
end
values = fields_size(:);
grps = fields_bat_num(:);

% plot
switch plot_style
    case 'errorbars'
        x = 1:length(bats);
        y = accumarray(grps, values, [], @nanmean, nan);
        y = y(bats);
        err = accumarray(grps, values, [], @nansem, nan);
        err = err(bats);
        h=bar(x,y);
        h.FaceColor = 0.5*[1 1 1];
        h=errorbar(x,y,err);
        h.LineStyle='none';
        h.Color = 'k';
        set(gca,'XTick',1:length(bats), 'XTickLabel', bats);
    case 'boxplot'
        boxplot(values, grps);
    case 'violin'
        rng(0);
        h=violinplot(values, grps);
        [h.ViolinColor] = disperse( bats_colors );
end

% set axis properties
ha = gca;
% ha.YScale = plot_yscale;
ha.TickDir='out';
ha.TickLength = [0.01 0.01];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
ylabel('Field size', 'Units','normalized','Position',[-0.05 0.5]);
xlabel('Bat no.', 'Units','normalized','Position',[0.5 -0.13]);

%% 3) panels C - LS ratio
axes(panel_C);
cla
hold on
text(-0.1,1.1, 'C', 'Units','normalized','FontWeight','bold');

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
values = LS_field_ratio_all(:);
grps = cells_bat_num(:);

% plot
switch plot_style 
    case 'errorbars'
        x = 1:length(bats);
        y = accumarray(grps, values, [], @nanmean, nan);
        y = y(bats);
        err = accumarray(grps, values, [], @nansem, nan);
        err = err(bats);
        h=bar(x,y);
        h.FaceColor = 0.5*[1 1 1];
        h=errorbar(x,y,err);
        h.LineStyle='none';
        h.Color = 'k';
        set(gca,'XTick',1:length(bats), 'XTickLabel', bats);
    case 'boxplot'
        boxplot(values, grps);
    case 'violin'
        rng(0);
        h=violinplot(values, grps);
        [h.ViolinColor] = disperse( bats_colors );
end

% set axis properties
ha = gca;
ha.YScale = plot_yscale;
ha.TickDir='out';
ha.TickLength = [0.01 0.01];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
ylabel({'Field size ratio';'largest/smallest'}, 'Units','normalized','Position',[-0.05 0.5]);
xlabel('Bat no.', 'Units','normalized','Position',[0.5 -0.13]);


%% print/save the figure
fig_name_out = fullfile(res_dir, sprintf('%s__%s__%s',fig_name_str, plot_style, plot_yscale));
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

