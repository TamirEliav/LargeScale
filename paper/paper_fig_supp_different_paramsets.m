%% Large Scale - Fig. S7 - robustness of results with different paramsets

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'Fig_S7';
fig_caption_str = 'robustness of results with different paramsets';
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
panels_size = [3.5 2];
paramsets = 0:7;
paramsets(6) = 8;
x_positions =      linspace(2,16,4);
x_positions(end)=[];
y_positions = flip(linspace(3,23,length(paramsets)));
for ii_paramset = 1:length(paramsets)
    y = y_positions(ii_paramset);
    for ii_x = 1:length(x_positions)
        x = x_positions(ii_x);
        panel_AB(ii_paramset,ii_x) = axes('position', [ x y panels_size]);
    end
end

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

pop_data = struct();
for ii_paramset = 1:length(paramsets)
    paramset = paramsets(ii_paramset);
    fprintf('loading paramset %d results...\n',paramset)
    PARAMS_SetParamset(paramset);
    pop_data(ii_paramset).paramset = paramset;
    cells = cellfun(@(c)(cell_load_data(c,'details','stats','signif','fields')), cells_ID, 'UniformOutput',0);
    cells = [cells{:}];
    pop_data(ii_paramset).cells = cells;
end
% go back to default paramset (0)
disp('Go back to default paramset (0)')
PARAMS_SetParamset(0);

%% 1) panels A+B - No. of fields
nFields = nan(length(paramsets),2,length(cells));
for ii_paramset = 1:length(paramsets)
    %% choose axis
    axes(panel_AB(ii_paramset,1));
    cla
    hold on
    %% arrange data
    cells = pop_data(ii_paramset).cells;
    for ii_dir = 1:2
        for ii_cell = 1:length(cells)
            cell = cells(ii_cell);
            % check signif per direction
            if ~cell.signif(ii_dir).TF
                continue;
            end
            nFields(ii_paramset,ii_dir,ii_cell) = cell.stats.dir(ii_dir).field_num;
        end
    end
    %% plot
    nFields2plot = nFields(ii_paramset,:,:);
    h = histogram(nFields2plot(:));
    h.FaceColor = 0.5*[1 1 1];
    nBinEdges = 14;
    h.BinEdges = linspace(0,35,nBinEdges);
    ha = gca;
    ha.YScale = 'log';
    ha.YLim = [0.7 max(h.Values)*1.05];
    ha.XLim = [0 33];
    ha.XTick = [0:10:30];
    ha.YTick = [1 10 100];
    ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
    ha.TickDir='out';
    ha.TickLength = [0.03 0.03];
    ha.XRuler.TickLabelGapMultiplier = -0.3;
    ha.YRuler.TickLabelGapMultiplier = 0.001;
    ylabel('No. of cells')
    switch ii_paramset
        case {length(paramsets)}
            xlabel('No. of fields per direction', 'Units','normalized','Position',[0.5 -0.18]);
    end
end


%% 2) panels A+B - Field size
fields_size_all = {};
for ii_paramset = 1:length(paramsets)
    %% choose axis
    axes(panel_AB(ii_paramset,2));
    cla
    hold on
    %% arrange data
    cells = pop_data(ii_paramset).cells;
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
    fields_size_all{ii_paramset} = fields_size;

    %%  plot
    h = histogram(fields_size);
    h.FaceColor = 0.5*[1 1 1];
    h.BinWidth = 1;
    h.BinLimits=[0 40];
    ha=gca;
    ha.YScale = 'log';
    ha = gca;
    ha.YLim = [0.8 max(h.Values)*1.05];
    ha.XLim = [-0.5 40];
    ha.YTick = [1 10 100];
    ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
    ha.TickDir='out';
    ha.TickLength = [0.03 0.03];
    ha.XRuler.TickLabelGapMultiplier = -0.3;
    ha.YRuler.TickLabelGapMultiplier = 0.001;
    ylabel('Counts','Units','normalized','Position',[-0.18 0.5])
    switch ii_paramset
        case {length(paramsets)}
            xlabel('Field size (m)','Units','normalized','Position',[0.5 -0.18]);
    end
end


%% 3) panels A+B - LS ratio
LS_field_ratio_all = nan( length(paramsets),    length(cells) );
LS_field_ratio_dir = nan( length(paramsets), 2, length(cells) );
for ii_paramset = 1:length(paramsets)
    %% choose axis
    axes(panel_AB(ii_paramset,3));
    cla
    hold on
    %% arrange data
    cells = pop_data(ii_paramset).cells;
    for ii_cell = 1:length(cells)
        cell = cells(ii_cell);
        % pooled stats - check at least one direction is signif
        if any([cell.signif.TF])
            LS_field_ratio_all(ii_paramset,ii_cell) = cell.stats.all.field_ratio_LS;
        end
        % per dir stats - check signif per direction
        for ii_dir = 1:2
            if cell.signif(ii_dir).TF 
                LS_field_ratio_dir(ii_paramset,ii_dir,ii_cell) = cell.stats.dir(ii_dir).field_ratio_LS;
            end
        end
    end

    %%  plot
    nBinEdges = 9;
    edges = logspace(0,log10(25),nBinEdges);
    LS_field_ratio_all_2plot = LS_field_ratio_all(ii_paramset,:);
    h = histogram(LS_field_ratio_all_2plot);
    h.BinEdges = edges;
    h.FaceColor = 0.5*[1 1 1];
    ha=gca;
    ha.YScale = 'log';
    ha.XScale = 'log';
    ha.XLim = [0 27];
    ha.YLim = [7e-1 260];
    ha.YTick = [1 10 100];
    ha.XTick = [1 2 5 10 20];
    ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
    ha.TickDir='out';
    ha.TickLength = [0.03 0.03];
    ha.XRuler.TickLabelGapMultiplier = -0.35;
    ha.YRuler.TickLabelGapMultiplier = 0.001;
    ylabel('No. of cells','Units','normalized','Position',[-0.18 0.5])
    switch ii_paramset
        case {length(paramsets)}
            xlabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.17]);
    end
end


%% 4) panels A+B - LS ratio correlation to paramset 0
if 0
for ii_paramset = 1:length(paramsets)
    %% choose axis
    axes(panel_AB(ii_paramset,4));
    cla
    hold on
    axis equal
    axis square
    plot(LS_field_ratio_all(paramsets==0,:),...
         LS_field_ratio_all(ii_paramset,:), '.k')
    h=refline(1,0); h.Color = 'k';
    ha=gca;
    ha.YScale = 'log';
    ha.XScale = 'log';
    ha.XLim = [0 27];
    ha.YLim = [0 27];
    ha.XTick = [1 2 5 10 20];
    ha.YTick = [1 2 5 10 20];
    ha.TickDir='out';
    ha.TickLength = [0.03 0.03];
    ha.XRuler.TickLabelGapMultiplier = -0.35;
    ha.YRuler.TickLabelGapMultiplier = 0.001;
    ylabel(sprintf('paramset %d',paramsets(ii_paramset)), 'Units','normalized','Position',[-0.24 0.5]);
    switch ii_paramset
        case {length(paramsets)}
            xlabel(sprintf('paramset %d',0), 'Units','normalized','Position',[0.5 -0.17]);
        otherwise
            xlabel('')
    end
end
end

%% Add Paramsets labels (letters A-H)
paramsets_letters = char('A'+[0:length(paramsets)-1]);
for ii_paramset = 1:length(paramsets)
    axes(panel_AB(ii_paramset,1));
    text(-0.4,1, paramsets_letters(ii_paramset), 'Units','normalized','FontWeight','bold','VerticalAlignment','middle');
end

%% print/save the figure
fig_name_out = fullfile(res_dir, [fig_name_str sprintf('_%d',paramsets)]);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

