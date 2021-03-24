%% Large Scale - fig 5 - maps properties vs exposure days + dynamic examples

%%
clear 
clc

%% plotting options
field_speed_opt = 1;
% corr_type = 'pearson';
corr_type = 'spearman';

dynamic_examples_option = 5;
rng(0);

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'Fig_5';
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
panel_A_size_raster = [3.75 1];
panel_A_size_FR_map = [3.75 1];
panel_A_pos = [2 20];
panel_A = [];
for ii=1:4    
    offset_x = (ii-1)*4.9;
    offset_y = 0;
    offset = panel_A_pos + [offset_x offset_y];
    panel_A(ii,1) = axes('units', 'centimeters', 'position', [offset+[0 2.25] panel_A_size_FR_map]);
    panel_A(ii,2) = axes('units', 'centimeters','position', [offset+[0 1   ] panel_A_size_raster]);
    panel_A(ii,3) = axes('units', 'centimeters','position', [offset+[0 0   ] panel_A_size_raster]);
end
% % panel_A = permute(panel_A,[2 1 3]);
% panel_A = panel_A(:,3:-1:1,:);
% panel_A = reshape(panel_A,[9 3]);

panels_size = [4 2.5];
panel_B    = axes('position', [ 2.0  15  panels_size           ]);
panel_C    = axes('position', [ 8.37  15  panels_size           ]);
panel_D = axes('position', [ 15.1  15  panels_size              ]);
panels_size = [0.9 2.5];
panel_B_cmp    = axes('position', [ 6.1  15  panels_size           ]);
panel_C_cmp    = axes('position', [ 12.47  15  panels_size           ]);
panel_D_cmp = axes('position', [ 19.2  15  panels_size              ]);

panel_E_size_raster = [3.75 1.5];
panel_E_pos = [2 9.5];
panel_E = [];
for ii=1:4    
    offset_x = (ii-1)*4.9;
    offset_y = 0;
    offset = panel_E_pos + [offset_x offset_y];
    panel_E(ii,1) = axes('position', [offset+[0 1.5   ] panel_E_size_raster]);
    panel_E(ii,2) = axes('position', [offset+[0 0   ] panel_E_size_raster]);
end

panel_F = axes('position', [ 2.0    3.7     5     2.5 ]);
panel_F_legend = axes('position', [ 7.0    5.4     2     1 ]);
% panel_G = axes('position', [ 8.5    4.3     2     2.5 ]);

%%
% prm = PARAMS_GetAll();
prm = PARAMS_GetAll__shir_version();

%% load data
addpath('L:\processed_data_structs');
%load 6m data:
cells_6m = load('cells_bat_6m.mat');
cells_6m = cells_6m.cells;
% remove non-signif cells:
signif = cat(1,cells_6m.signif_6m);
signif = arrayfun(@(x)(x.TF), signif);
signif = any(signif,2);
cells_6m(~signif) = [];

%load 200m data:
cells_200m = load('cells_bat_200m.mat');
cells_200m = cells_200m.cells;

% %%
% prm = PARAMS_GetAll();
% cells_t = DS_get_cells_summary();
% bats = [2311,2382];
% cells_t(~ismember(cells_t.bat, bats ),:) = [];
% cells_t([cells_t.remove]==1,:) = []; %cells chosen to remove because of repititions between days
% cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
% cells = [cells{:}];
% cells = [cells.details];
% cells(~contains({cells.brain_area}, 'CA1')) = [];
% cells(ismember([cells.DayNum], [0])) = []; %take only from day1 and on data
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
% cells = cellfun(@(c)(cell_load_data(c,'details','stats','meanFR','stats','inclusion','signif','fields','fields_per_win','FR_map','FE')), cells_ID, 'UniformOutput',0);
% cells = [cells{:}];
% for ii_cell = 1:length(cells)
%     exp = exp_load_data(cells(ii_cell).details.exp_ID,'balls');
%     cells(ii_cell).balls = exp.balls;
% end
% save('Cells_130m.mat','cells');

%%
cells_130m = load('cells_bat_130m.mat');
cells_130m = cells_130m.cells;
cells_130m_IDs = arrayfun(@(x)(x.cell_ID), cat(1,cells_130m.details), 'UniformOutput', false);



%% FR map + rasters - 4 examples

cell_examples = {
    'b2311_d191218_TT3_SS02';
    'b2382_d190624_TT13_SS03';
    'b2382_d190627_TT13_SS02';
    'b2311_d191224_TT7_SS01'};

if 1
for ii_cell = 1:length(cell_examples)
    %%
    example_IX = find(strcmp(cells_130m_IDs,cell_examples{ii_cell}));
    cell = cells_130m(example_IX);
    
%     cell_ID = cell_examples{ii_cell};
%     cell = cell_load_data(cell_ID,'details','FR_map','fields','stats','FE');
%     exp = exp_load_data(cell.details.exp_ID,'balls');
    
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
    valid_speed_pos = prm.fields.valid_speed_pos+cell.balls;
    area([cell.balls(1) valid_speed_pos(1)]    , repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
    area([  valid_speed_pos(2) max(cell.balls(2),130)], repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
    
%     valid_speed_pos = prm.fields.valid_speed_pos+exp.balls;
% %     area([4 valid_speed_pos(1)]    , repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
% %     area([  valid_speed_pos(2) 135], repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
%     area([exp.balls(1) valid_speed_pos(1)]    , repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
%     area([  valid_speed_pos(2) max(exp.balls(2),130)], repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
%     
    h=plot(x,y);
    [h.Color] = disperse(c);
    box off
    h=gca;
    h.TickDir = 'out';
    h.XTick = [];
    h.YTick = [0 m];
    h.YLim = ylimits;
    h.XLim = [0 130];
    h.YRuler.TickLabelGapOffset = 1;
    h.Clipping = 'off';
    
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
    cell_num_str_pos_x   = [0.50 0.50 0.50 0.50];
    cell_num_str_pos_y   = [1.05 1.05 1.05 1.05];
    cell_stats_str_pos_x = [0.75 0.93 0.95 0.90];
    cell_stats_str_pos_y = [1.00 1.05 1.00 1.00]+0.05;
    text(cell_num_str_pos_x(ii_cell), cell_num_str_pos_y(ii_cell), "Cell "+ ii_cell+" - Day "+cell.details.DayNum,...
        'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
%     h = text(cell_num_str_pos_x(ii_cell), cell_num_str_pos_y(ii_cell), cell_ID+" - Day "+cell.details.DayNum,...
%         'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
%     h.Interpreter='none';
    switch ii_cell
        case {1,2,3,4}
            cell_stats_str = {  sprintf('max=%.1fm', cell.stats.all.field_largest);...
                                sprintf('min=%.1fm', cell.stats.all.field_smallest);...
                                sprintf('ratio=%.1f', cell.stats.all.field_ratio_LS);...
                                %sprintf('Iso dis=%.1f', cell.stats.all.IsoDist);...
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
        h.XLim = [0 130];
        h.YLim = [0 m+1];
        h.Clipping = 'off';
        switch ii_dir
            case 1
                h.XTick = [];
                h.YTickLabel = {'',num2str(m)};
                h.TickLength = [0.0133 0.025];
                h.YRuler.TickLabelGapOffset = 1;
                h.XColor = [1 1 1];
            case 2
                h.XTick = [0,130];
                h.TickLength = [0.0133 0.025];
                h.XRuler.TickLabelGapOffset = -1;
                h.YRuler.TickLabelGapOffset = 1;
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
        text(1.07, dir_offsets(ii_dir), num2str(length(fields)), 'Units','normalized', 'Color', c{ii_dir},...
            'HorizontalAlignment','right', 'VerticalAlignment','middle','FontSize',6);
    end
end

%% add x/y labels for specific panels
for ii = [1 2 3 4]
    axes(panel_A(ii, 3));
    xlabel('Position (m)', 'Units','normalized','Position',[0.5 -0.15]);
end
for ii = [1]
    axes(panel_A(ii, 3));
%     ylabel('Time (min)',   'Units','normalized','Position',[-0.1 1]);
    ylabel('Flight no.',   'Units','normalized','Position',[-0.15 1]);
    axes(panel_A(ii, 1));
    ylabel({'Firing rate';'(Hz)'},   'Units','normalized','Position',[-0.1 0.42]);
end
axes(panel_A(1, 1));
text(-0.33,1.8, 'A', 'Units','normalized','FontWeight','bold');
text(0,1.8, 'Examples of multiscale coding in the first days of exposure to the tunnel',...
    'Units','normalized','FontWeight','bold');

%% add direction arrows
arrow_x = 0.85 +[0 0.05];
arrow_y = repelem(0.9111,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.008,  'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y      ,  'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);

end







% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% add day num data

day_num_perCell = arrayfun(@(x)(x.DayNum), cat(1,cells_130m.details));
bat_num_perCell = arrayfun(@(x)(x.bat), cat(1,cells_130m.details));
day_num_perCell_perDir = [day_num_perCell day_num_perCell];
day_num_perCell_perDir = day_num_perCell_perDir(:);

day_bin_size = 5;
bin_edges = (0:day_bin_size:max(day_num_perCell)+day_bin_size) + 0.5;
[~,~,day_bin_perCell] = histcounts(day_num_perCell,bin_edges);
[~,~,day_bin_perCell_perDir] = histcounts(day_num_perCell_perDir,bin_edges);

Daylimits = [0, max(day_num_perCell)+1]; 

%% revision edit - test directionality vs. days
signif = arrayfun(@(x)(x.TF), cat(1,cells_130m.signif));
signif = all(signif,2);
cells = cells_130m(signif);
FR_maps_all = cat(1,cells.FR_map);
FR_maps_all = reshape([FR_maps_all.all],size(FR_maps_all,1),size(FR_maps_all,2),[]);
M = cat(1,FR_maps_all.PSTH);
M = reshape(M,size(FR_maps_all,1),size(FR_maps_all,2),[]);
ccc = corr(squeeze(M(:,1,:))', squeeze(M(:,2,:))' ,'rows', 'pairwise');
ccc = diag(ccc);
x = day_num_perCell(signif);
y = ccc;
[rho,rho_pval] = corr(x,y,'type','Spearman');

fprintf('Directionality vs. days:\n')
fprintf('rho = %.2f\n',rho);
fprintf('p = %.2f\n',rho_pval);

% figure
% plot(x, y,'o');
% xlabel('Days');
% ylabel('directionality (maps corr)');
% text(0.1,0.85,sprintf('rho = %.2f',rho),'Units','normalized');
% text(0.1,0.8,sprintf('P = %.g',rho_pval),'Units','normalized');
% title('Directionality vs. days of exposure')


%% panel B - field count vs days
% figure
axes(panel_B);
cla
hold on
text(-0.3125,1.1, 'B', 'Units','normalized','FontWeight','bold');
text(0,1.35, 'Population: Stability of multiscale coding across weeks, starting from day 1 in the tunnel',...
    'Units','normalized','FontWeight','bold');

signif = arrayfun(@(x)(x.TF), cat(1,cells_130m.signif));
n_fields = arrayfun(@(x)([x.dir.field_num]), cat(1,cells_130m.stats),'UniformOutput',0);
n_fields = cat(1,n_fields{:});
n_fields(~signif) = nan;
n_fields = n_fields(:);

jitter_value_y = 0.3;
jitter_value_x = 0.5;
n_fields_jitter = n_fields + jitter_value_y*(2*rand([length(n_fields),1])-1);
day_num_perCell_perDir_jitter = day_num_perCell_perDir + jitter_value_x*(2*rand([length(day_num_perCell_perDir),1])-1);
p = plot(day_num_perCell_perDir_jitter, n_fields_jitter, '.');
p.Color = 0.65*[1 1 1];

day_bins = day_bin_perCell_perDir;
day_num_perCell_perDir_tmp = day_num_perCell_perDir;
day_num_perCell_perDir_tmp(isnan(n_fields),:) = [];
day_bins(isnan(n_fields),:) = [];
n_fields(isnan(n_fields)) = [];
[day_groups, day_bins] = findgroups(day_bins);
day_groups_stat = splitapply(@(x) [mean(x) std(x)],n_fields, day_groups);
day_num = day_bins*day_bin_size - (day_bin_size/2) + 0.5;
errorbar(day_num, day_groups_stat(:,1),day_groups_stat(:,2),'k','CapSize',2,'LineWidth',1);

ylabel({'No. of fields', 'per direction'},'Units','normalized','Position',[-0.15 0.5]);
xlabel('Day','Units','normalized','Position',[0.5 -0.18]);
ha = gca;
ha.YLim = [0,20];
ha.XLim = Daylimits;
ha.XTick = [0:10:30];
ha.YTick = [5 10 15 20];
% ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
ha.TickDir='out';
ha.TickLength = [0.02 0.02];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

[rho,p_val] = corr(n_fields,day_num_perCell_perDir_tmp,'type',corr_type);
text(0.96,1.05, ['\rho' sprintf(' = %.02f',rho)],'Units','normalized','FontSize',7,'HorizontalAlignment','right');
text(0.96,0.92, sprintf('P = %.02f',p_val),'Units','normalized','FontSize',7,'HorizontalAlignment','right');
text(0.05,0.95, sprintf('n = %d',length(n_fields)),'Units','normalized','FontSize',7,'HorizontalAlignment','left');

%% panel B - compare field count
axes(panel_B_cmp);
cla
hold on

nFields_200m = nan(2,length(cells_200m));
for ii_dir = 1:2
    for ii_cell = 1:length(cells_200m)
        cell = cells_200m(ii_cell);
        % check signif per direction
        if ~cell.signif(ii_dir).TF
            continue;
        end
        nFields_200m(ii_dir,ii_cell) = cell.stats.dir(ii_dir).field_num;
    end
end
nFields_130m = nan(2,length(cells_130m));
for ii_dir = 1:2
    for ii_cell = 1:length(cells_130m)
        cell = cells_130m(ii_cell);
        % check signif per direction
        if ~cell.signif(ii_dir).TF
            continue;
        end
        nFields_130m(ii_dir,ii_cell) = cell.stats.dir(ii_dir).field_num;
    end
end
nFields_6m = nan(2,length(cells_6m));
for ii_dir = 1:2
    for ii_cell = 1:length(cells_6m)
        cell = cells_6m(ii_cell);
        % check signif per direction
        if ~cell.signif_6m(ii_dir).TF
            continue;
        end
        nFields_6m(ii_dir,ii_cell) = cell.stats_6m.dir(ii_dir).field_num;
    end
end

y_bar = [nanmean(nFields_6m(:)), nanmean(nFields_130m(:)), nanmean(nFields_200m(:))];
sem_bar = [nanstd(nFields_6m(:)), nanstd(nFields_130m(:)), nanstd(nFields_200m(:))]./ ...
    sqrt([length(~isnan(nFields_6m(:))), length(~isnan(nFields_130m(:))), length(~isnan(nFields_200m(:)))]);
hb = bar([1,2,3],y_bar);
hb.FaceColor = 'flat';
hb.CData = [0.65 0.65 0.65];
he=errorbar([1,2,3],y_bar, sem_bar);
he.CapSize = 2;
he.LineWidth = 1;
he.LineStyle = 'none';
he.Color = 'k';
h=gca;
h.XTick = [0.85,1.85,2.9];
h.XTickLabels = {'  6m','130m','200m'};
h.XTickLabelRotation = 90;
ylimits = [0,20];
h.YLim = ylimits;
% h.YTick = ylimits;
h.XLim = [0.25,3.75];
h.TickLength = [0.01 0.01];
h.XAxis.TickLength = [0 0];
h.XRuler.TickLabelGapMultiplier = -0.001;
h.YRuler.TickLabelGapMultiplier = 0.001;
h.YTick = [];
h.YColor = 'none';


%% panel C - field size vs days
% figure
axes(panel_C);
cla
hold on
text(-0.25,1.1, 'C', 'Units','normalized','FontWeight','bold');
fields_size = [];
day_num_fields = [];
for ii_dir = 1:2
    for ii_cell = 1:length(cells_130m)
        cell = cells_130m(ii_cell);
        if ~cell.signif(ii_dir).TF % check signif per direction
            continue;
        end
        fields = cell.fields{ii_dir};
        fields([fields.in_low_speed_area]) = []; % remove fields in low speed area
        fields_size = [fields_size fields.width_prc];
        day_num_fields = [day_num_fields cell.details.DayNum*ones(1,length(fields))];
    end
end

jitter_value_x = 0.5;
day_num_fields_jitter = day_num_fields + jitter_value_x*(2*rand([1,length(day_num_fields)])-1);
p = plot(day_num_fields_jitter, fields_size, '.');
p.Color = 0.65*[1 1 1];

[~,~,day_bin_perField] = histcounts(day_num_fields,bin_edges);
[day_groups, day_bins] = findgroups(day_bin_perField);
day_groups_stat = splitapply(@(x) [mean(x) std(x)],fields_size', day_groups');
day_num = day_bins*day_bin_size - (day_bin_size/2) + 0.5;
errorbar(day_num, day_groups_stat(:,1),day_groups_stat(:,2),'k','CapSize',2,'LineWidth',1);

ylabel('Field size (m)','Units','normalized','Position',[-0.15 0.5])
xlabel('Day','Units','normalized','Position',[0.5 -0.18]);
ha = gca;
ha.XLim = Daylimits;
ha.YLim = [0 25];
ha.XTick = [0:10:30];
ha.YTick = [5 10 15 20 25];
ha.TickDir='out';
ha.TickLength = [0.02 0.02];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

[rho,p_val] = corr(fields_size',day_num_fields','type',corr_type);
text(0.96,1.05, ['\rho' sprintf(' = %.02f',rho)],'Units','normalized','FontSize',7,'HorizontalAlignment','right');
text(0.96,0.92, sprintf('P = %.02f',p_val),'Units','normalized','FontSize',7,'HorizontalAlignment','right');
text(0.05,0.95, sprintf('n = %d',length(fields_size)),'Units','normalized','FontSize',7,'HorizontalAlignment','left');

% save( fullfile(res_dir,'pop_dist_fields_size'), 'fields_size');

%% panel C - field size compare setups
axes(panel_C_cmp);
cla
hold on

fields_size_200m = [];
for ii_dir = 1:2
    for ii_cell = 1:length(cells_200m)
        cell = cells_200m(ii_cell);
        if ~cell.signif(ii_dir).TF % check signif per direction
            continue;
        end
        fields = cell.fields{ii_dir};
        fields([fields.in_low_speed_area]) = []; % remove fields in low speed area
        fields_size_200m = [fields_size_200m fields.width_prc];
    end
end

fields_size_130m = [];
for ii_dir = 1:2
    for ii_cell = 1:length(cells_130m)
        cell = cells_130m(ii_cell);
        if ~cell.signif(ii_dir).TF % check signif per direction
            continue;
        end
        fields = cell.fields{ii_dir};
        fields([fields.in_low_speed_area]) = []; % remove fields in low speed area
        fields_size_130m = [fields_size_130m fields.width_prc];
    end
end

fields_size_6m = [];
for ii_dir = 1:2
    for ii_cell = 1:length(cells_6m)
        cell = cells_6m(ii_cell);
        if ~cell.signif_6m(ii_dir).TF % check signif per direction
            continue;
        end
        fields = cell.fields_6m{ii_dir};
        fields([fields.in_low_speed_area]) = []; % remove fields in low speed area
        fields_size_6m = [fields_size_6m fields.width_prc];
    end
end

y_bar = [mean(fields_size_6m), mean(fields_size_130m), mean(fields_size_200m)];
sem_bar = [std(fields_size_6m), std(fields_size_130m), std(fields_size_200m)]./ ...
    sqrt([length(fields_size_6m), length(fields_size_130m), length(fields_size_200m)]);
hb = bar([1,2,3],y_bar);
hb.FaceColor = 'flat';
hb.CData = [0.65 0.65 0.65];
he=errorbar([1,2,3],y_bar, sem_bar);
he.CapSize = 2;
he.LineWidth = 1;
he.LineStyle = 'none';
he.Color = 'k';
h=gca;
h.XTick = [0.85,1.85,2.9];
h.XTickLabels = {'  6m','130m','200m'};
h.XTickLabelRotation = 90;
ylimits = [0,25];
h.YLim = ylimits;
% h.YTick = ylimits;
h.XLim = [0.25,3.75];
h.TickLength = [0.01 0.01];
h.XAxis.TickLength = [0 0];
h.XRuler.TickLabelGapMultiplier = -0.001;
h.YRuler.TickLabelGapMultiplier = 0.001;
h.YTick = [];
h.YColor = 'none';

%% panel D - field ratio (largest/smallest)
axes(panel_D);
cla
hold on
text(-0.31,1.1, 'D', 'Units','normalized','FontWeight','bold');
LS_field_ratio_all = nan(1,length(cells_130m));
LS_field_ratio_dir = nan(2,length(cells_130m));
for ii_cell = 1:length(cells_130m)
    cell = cells_130m(ii_cell);
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
jitter_value_x = 0.5;
day_num_perCell_jitter = day_num_perCell + jitter_value_x*(2*rand([length(day_num_perCell),1])-1);
p = plot(day_num_perCell_jitter, LS_field_ratio_all,'.');
p.Color = 0.65*[1 1 1];

day_bins = day_bin_perCell;
day_num_perCell_tmp = day_num_perCell;
day_num_perCell_tmp(isnan(LS_field_ratio_all),:) = [];
day_bins(isnan(LS_field_ratio_all),:) = [];
LS_field_ratio_all(isnan(LS_field_ratio_all)) = [];
[day_groups, day_bins] = findgroups(day_bins);
day_groups_stat = splitapply(@(x) [mean(x) std(x)],LS_field_ratio_all', day_groups);
day_num = day_bins*day_bin_size - (day_bin_size/2) + 0.5;
errorbar(day_num, day_groups_stat(:,1),day_groups_stat(:,2),'k','CapSize',2,'LineWidth',1);

ha=gca;
ha.XLim = Daylimits;
ha.YLim = [0 10];
ha.XTick = [0:10:30];
ha.YTick = [2:2:10];
ha.TickDir='out';
ha.TickLength = [0.02 0.02];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
ylabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[-0.14 0.5]);
xlabel('Day','Units','normalized','Position',[0.5 -0.18]);

[rho,p_val] = corr(LS_field_ratio_all',day_num_perCell_tmp,'type',corr_type);
text(0.96,1.05, ['\rho' sprintf(' = %.02f',rho)],'Units','normalized','FontSize',7,'HorizontalAlignment','right');
text(0.96,0.92, sprintf('P = %.02f',p_val),'Units','normalized','FontSize',7,'HorizontalAlignment','right');
text(0.05,0.95, sprintf('n = %d',length(LS_field_ratio_all)),'Units','normalized','FontSize',7,'HorizontalAlignment','left');

%% panel D - field ratio compare setups
axes(panel_D_cmp);
cla
hold on

LS_ratio_all_200m = nan(1,length(cells_200m));
LS_ratio_dir_200m = nan(2,length(cells_200m));
for ii_cell = 1:length(cells_200m)
    cell = cells_200m(ii_cell);
    % pooled stats - check at least one direction is signif
    if any([cell.signif.TF])
        LS_ratio_all_200m(ii_cell) = cell.stats.all.field_ratio_LS;
    end
    % per dir stats - check signif per direction
    for ii_dir = 1:2
        if cell.signif(ii_dir).TF 
            LS_ratio_dir_200m(ii_dir,ii_cell) = cell.stats.dir(ii_dir).field_ratio_LS;
        end
    end
end

LS_ratio_all_130m = nan(1,length(cells_130m));
LS_ratio_dir_130m = nan(2,length(cells_130m));
for ii_cell = 1:length(cells_130m)
    cell = cells_130m(ii_cell);
    % pooled stats - check at least one direction is signif
    if any([cell.signif.TF])
        LS_ratio_all_130m(ii_cell) = cell.stats.all.field_ratio_LS;
    end
    % per dir stats - check signif per direction
    for ii_dir = 1:2
        if cell.signif(ii_dir).TF 
            LS_ratio_dir_130m(ii_dir,ii_cell) = cell.stats.dir(ii_dir).field_ratio_LS;
        end
    end
end

LS_ratio_all_6m = nan(1,length(cells_6m));
LS_ratio_dir_6m = nan(2,length(cells_6m));
for ii_cell = 1:length(cells_6m)
    cell = cells_6m(ii_cell);
    % pooled stats - check at least one direction is signif
    if any([cell.signif_6m.TF])
        LS_ratio_all_6m(ii_cell) = cell.stats_6m.all.field_ratio_LS;
    end
    % per dir stats - check signif per direction
    for ii_dir = 1:2
        if cell.signif_6m(ii_dir).TF 
            LS_ratio_dir_6m(ii_dir,ii_cell) = cell.stats_6m.dir(ii_dir).field_ratio_LS;
        end
    end
end

y_bar = [nanmean(LS_ratio_all_6m), nanmean(LS_ratio_all_130m), nanmean(LS_ratio_all_200m)];
sem_bar = [nanstd(LS_ratio_all_6m), nanstd(LS_ratio_all_130m), nanstd(LS_ratio_all_200m)]./ ...
    sqrt([length(~isnan(LS_ratio_all_6m)), length(~isnan(LS_ratio_all_130m)), length(~isnan(LS_ratio_all_200m))]);
hb = bar([1,2,3],y_bar);
hb.FaceColor = 'flat';
hb.CData = [0.65 0.65 0.65];
he=errorbar([1,2,3],y_bar, sem_bar);
he.CapSize = 2;
he.LineWidth = 1;
he.LineStyle = 'none';
he.Color = 'k';
h=gca;
h.XTick = [0.85,1.85,2.9];
h.XTickLabels = {'  6m','130m','200m'};
h.XTickLabelRotation = 90;
ylimits = [0,10];
h.YLim = ylimits;
% h.YTick = ylimits;
h.XLim = [0.25,3.75];
h.TickLength = [0.01 0.01];
h.XAxis.TickLength = [0 0];
h.XRuler.TickLabelGapMultiplier = -0.001;
h.YRuler.TickLabelGapMultiplier = 0.001;
h.YTick = [];
h.YColor = 'none';

%% Panel E: FR map + rasters - 3 examples of dynamic during sessions

cell_examples = {
    'b2382_d190623_TT13_SS05';
    'b2311_d191220_TT3_SS01';
    'b2382_d190808_TT15_SS01';
    'b2382_d190808_TT15_SS04'};

for ii_cell = 1:length(cell_examples)
    %%
    example_IX = find(strcmp(cells_130m_IDs,cell_examples{ii_cell}));
    cell = cells_130m(example_IX);
    
%     cell_ID = cell_examples{ii_cell};
%     cell = cell_load_data(cell_ID,'details','FR_map','fields','stats','FE');
%     exp = exp_load_data(cell.details.exp_ID,'balls');

    c = prm.graphics.colors.flight_directions;
    
    % rasters
    FEs = [cell.FE];
    for ii_dir=1:2
        axes(panel_E(ii_cell, ii_dir));
        cla
        FE = FEs{ii_dir};
        x = [FE.spikes_pos];
        pos_start = arrayfun(@(FE)(FE.pos(1)),FE,'UniformOutput',0);
        pos_end = arrayfun(@(FE)(FE.pos(end)),FE,'UniformOutput',0);
        pos_start_end = [[pos_start{:}]; [pos_end{:}]];
%         pos = [FE.pos];
%         [FE.number2] = disperse(1:length(FE));
%         y = arrayfun(@(FE)(FE.number2*ones(1,length(FE.pos))),FE,'UniformOutput',0);
%         y = [y{:}];
        y = [(1:length(FE)); (1:length(FE))];
%         plot(pos,y,'.','Color',0.93*[1 1 1],'MarkerSize',0.005)
        plot(pos_start_end,y,'-','Color',0.93*[1 1 1],'linewidth',0.5)
        hold on;
        
        [FE.number2] = disperse(1:length(FE));
        y = arrayfun(@(FE)(FE.number2*ones(1,FE.num_spikes)),FE,'UniformOutput',0);
        y = [y{:}];
        plot(x,y,'.','Color',c{ii_dir},'MarkerSize',0.05);
        box off
        h=gca;
        m = length(FE);
        h.YTick = [1 m];
        h.XLim = [0 130];
        h.YLim = [0 m+1];
        h.Clipping = 'off';
        switch ii_dir
            case 1
                h.XTick = [];
                h.YTickLabel = {'',num2str(m)};
                h.TickLength = [0.0133 0.025];
                h.YRuler.TickLabelGapOffset = 1;
            case 2
                h.XTick = [0,130];
                h.TickLength = [0.0133 0.025];
                h.XRuler.TickLabelGapOffset = -1;
                h.YRuler.TickLabelGapOffset = 1;
                h.YTickLabel = {'1',num2str(m)};
                h.TickDir = 'out';
        end
    end
    % cell details
    cell_num_str_pos_x   = [0.50 0.50 0.50 0.50];
    cell_num_str_pos_y   = [2.05 2.05 2.05 2.05];
    cell_stats_str_pos_x = [0.90 0.90 0.90 0.90];
    cell_stats_str_pos_y = [1.10 1.10 1.10 1.10]+0.05;
    text(cell_num_str_pos_x(ii_cell), cell_num_str_pos_y(ii_cell), "Cell "+ (ii_cell+4) +" - Day "+cell.details.DayNum,...
        'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
%     h = text(cell_num_str_pos_x(ii_cell), cell_num_str_pos_y(ii_cell), cell_ID+" - Day "+cell.details.DayNum,...
%         'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
%     h.Interpreter='none';
    switch cell_examples{ii_cell}
        case 'b2382_d190808_TT15_SS01'
            axes(panel_E(ii_cell, 1));
            hax = gca;
            harr = annotation('arrow', [0 1],[0 1]);
            harr.Parent = hax;
            harr.LineStyle = 'none';
            harr.HeadLength = 5;
            harr.HeadWidth = 5;
            harr.X = 41 + [6 0];
            harr.Y = 12 + [9 0];
            
            axes(panel_E(ii_cell, 2));
            hax = gca;
            for ii_arrow = 1:2
                harr(ii_arrow) = annotation('arrow', [0 1],[0 1]);
                harr(ii_arrow).Parent = hax;
                harr(ii_arrow).LineStyle = 'none';
                harr(ii_arrow).HeadLength = 5;
                harr(ii_arrow).HeadWidth = 5;
            end
            harr(1).X = 47 + [0 6];
            harr(1).Y = 6 + [9 0];
            harr(2).X = 78 + [6 0];
            harr(2).Y = 8 + [9 0];
            
        case 'b2382_d190808_TT15_SS04'
            axes(panel_E(ii_cell, 1));
            hax = gca;
            harr = annotation('arrow', [0 1],[0 1]);
            harr.Parent = hax;
            harr.LineStyle = 'none';
            harr.HeadLength = 5;
            harr.HeadWidth = 5;
            harr.X = 46 + [0 6];
            harr.Y = 16 + [9 0];
            harrW = annotation('arrow', [0 1],[0 1]);
            harrW.Parent = hax;
            harrW.LineStyle = 'none';
            harrW.Color = [1 1 1];
            harrW.HeadLength = 2;
            harrW.HeadWidth = 2;
            harrW.X = 44.4 + [0 6];
            harrW.Y = 16.5 + [9 0];
            
        case 'b2382_d190623_TT13_SS05'
            axes(panel_E(ii_cell, 1));
            hax = gca;
            for ii_arrow = 1:3
                harr(ii_arrow) = annotation('arrow', [0 1],[0 1]);
                harr(ii_arrow).Parent = hax;
                harr(ii_arrow).LineStyle = 'none';
                harr(ii_arrow).HeadLength = 5;
                harr(ii_arrow).HeadWidth = 5;
            end
            harr(1).X = 116 + [0 6];
            harr(1).Y = 15 + [9 0];
            harr(2).X = 94 + [0 6];
            harr(2).Y = 16 + [9 0];
            harr(3).X = 119 + [6 0];
            harr(3).Y = 5 + [9 0];
            harrW = annotation('arrow', [0 1],[0 1]);
            harrW.Parent = hax;
            harrW.LineStyle = 'none';
            harrW.Color = [1 1 1];
            harrW.HeadLength = 2;
            harrW.HeadWidth = 2;
            harrW.X = 120.7 + [6 0];
            harrW.Y = 5.4 + [9 0];
                
            axes(panel_E(ii_cell, 2));
            hax = gca;
            harr = annotation('arrow', [0 1],[0 1]);
            harr.Parent = hax;
            harr.LineStyle = 'none';
            harr.HeadLength = 5;
            harr.HeadWidth = 5;
            harr.X = 22 + [6 0];
            harr.Y = 5 + [9 0];
            
        case 'b2311_d191220_TT3_SS01'
            axes(panel_E(ii_cell, 1));
            hax = gca;
            for ii_arrow = 1:2
                harr(ii_arrow) = annotation('arrow', [0 1],[0 1]);
                harr(ii_arrow).Parent = hax;
                harr(ii_arrow).LineStyle = 'none';
                harr(ii_arrow).HeadLength = 5;
                harr(ii_arrow).HeadWidth = 5;
            end
            harr(1).X = 21 + [0 6];
            harr(1).Y = 12 + [9 0];
            harr(2).X = 20 + [0 6];
            harr(2).Y = 34.5 + [9 0];
            harrW = annotation('arrow', [0 1],[0 1]);
            harrW.Parent = hax;
            harrW.LineStyle = 'none';
            harrW.Color = [1 1 1];
            harrW.HeadLength = 2;
            harrW.HeadWidth = 2;
            harrW.X = 18.4 + [0 6];
            harrW.Y = 35 + [9 0];
            
            axes(panel_E(ii_cell, 2));
            hax = gca;
            harr = annotation('arrow', [0 1],[0 1]);
            harr.Parent = hax;
            harr.LineStyle = 'none';
            harr.HeadLength = 5;
            harr.HeadWidth = 5;
            harr.X = 28 + [0 6];
            harr.Y = 6 + [9 0];
            
    end
end

%% add x/y labels for specific panels
for ii = [1 2 3 4]
    axes(panel_E(ii, 2));
    xlabel('Position (m)', 'Units','normalized','Position',[0.5 -0.1]);
end
for ii = [1]
    axes(panel_E(ii, 2));
    ylabel('Flight no.',   'Units','normalized','Position',[-0.15 1]);
end
axes(panel_E(1, 2));
text(-0.33,2.6, 'E', 'Units','normalized','FontWeight','bold');
text(0,2.6, 'Examples of within-day dynamics in spatial tuning',...
    'Units','normalized','FontWeight','bold');


%% Panel F,G: probability for appear/disappear segment
day_th1 = 2;
day_th2 = 4;
day_th3 = 6;
% day_th4 = 8;
field_size_th = 4.3;%median(fields_size);

count_seg_change = [];
count_seg_change_part1 = [];
count_seg_change_part2 = [];
count_seg_change_part3 = [];
count_seg_change_part4 = [];
% count_seg_change_part5 = [];
count_seg_change_small_fields = [];
count_seg_change_large_fields = [];
count_start_end_small_fields = [];
count_start_end_large_fields = [];

count_map_part1 = [];
count_map_part2 = [];
count_map_part3 = [];
count_map_part4 = [];
flight_ID_app = [];
flight_ID_dis = [];
flight_ID_part1 = [];
flight_ID_part2 = [];
flight_ID_part3 = [];
flight_ID_part4 = [];

n_small_fields = 0;
n_large_fields = 0;
count_map = 0;
for ii_cell = 1:length(cells_130m)
    cell_ID = cells_130m_IDs{ii_cell};
    cell = cells_130m(ii_cell);
%     cell = cell_load_data(cell_ID, 'details','stats','fields_per_win','FE');

    for ii_dir = 1:2
        field_num = cell.stats.dir(ii_dir).field_num;
        map_signif = cell.stats.dir(ii_dir).map_signif;
        if field_num>0 && map_signif
            count_map = count_map+1;
            fields_per_win = cell.fields_per_win{ii_dir};
            FE = cell.FE{ii_dir};
            flight_ID_tmp = [];
            for ii_flight = 1:length(FE)
                flight_ID_tmp(ii_flight) = str2num([num2str(cell.details.DayNum) num2str(FE(ii_flight).start_IX)]);
            end
            [segment_appear_flight, segment_disappear_flight, ~] = map_detect_segments_appear_disappear(fields_per_win,FE);
            count_seg_change = [count_seg_change ; [segment_appear_flight' segment_disappear_flight']];
            flight_ID_app = [flight_ID_app ; flight_ID_tmp(2+(1:length(segment_appear_flight)))'];
            flight_ID_dis = [flight_ID_dis ; flight_ID_tmp((end-length(segment_appear_flight)):(end-1))'];
            
            if cell.details.DayNum <= day_th1
                count_seg_change_part1 = [count_seg_change_part1 ; [segment_appear_flight' segment_disappear_flight']];
                count_map_part1 = [count_map_part1;  count_map*ones([length(segment_appear_flight),1])];
                flight_ID_part1 = [flight_ID_part1; flight_ID_tmp(2+(1:length(segment_appear_flight)))'];
            elseif cell.details.DayNum <= day_th2
                count_seg_change_part2 = [count_seg_change_part2 ; [segment_appear_flight' segment_disappear_flight']];
                count_map_part2 = [count_map_part2;  count_map*ones([length(segment_appear_flight),1])];
                flight_ID_part2 = [flight_ID_part2; flight_ID_tmp(2+(1:length(segment_appear_flight)))'];
            elseif cell.details.DayNum <= day_th3
                count_seg_change_part3 = [count_seg_change_part3 ; [segment_appear_flight' segment_disappear_flight']];
                count_map_part3 = [count_map_part3;  count_map*ones([length(segment_appear_flight),1])];
                flight_ID_part3 = [flight_ID_part3; flight_ID_tmp(2+(1:length(segment_appear_flight)))'];
%             elseif cell.details.DayNum <= day_th4
%                 count_seg_change_part4 = [count_seg_change_part4 ; [segment_appear_flight' segment_disappear_flight']];
            else
%                 count_seg_change_part5 = [count_seg_change_part5 ; [segment_appear_flight' segment_disappear_flight']];
                count_seg_change_part4 = [count_seg_change_part4 ; [segment_appear_flight' segment_disappear_flight']];
                count_map_part4 = [count_map_part4;  count_map*ones([length(segment_appear_flight),1])];
                flight_ID_part4 = [flight_ID_part4; flight_ID_tmp(2+(1:length(segment_appear_flight)))'];
            end
            
            % sort fields by size:
            field_size_per_win = [fields_per_win(:).field_size];
            size_per_field = nanmedian(field_size_per_win);
            fields_per_win_small = fields_per_win;
            fields_per_win_large = fields_per_win;
            fields_per_win_small(size_per_field >= field_size_th) = [];
            fields_per_win_large(size_per_field < field_size_th) = [];
            n_small_fields = n_small_fields + sum(size_per_field < field_size_th);
            n_large_fields = n_large_fields + sum(size_per_field >= field_size_th);
            
            if sum(size_per_field < field_size_th)>0
                [segment_appear_flight, segment_disappear_flight, ~] = map_detect_segments_appear_disappear(fields_per_win_small,FE);
                count_seg_change_small_fields = [count_seg_change_small_fields ; [segment_appear_flight' segment_disappear_flight']];
                % only field start, end (to avoid field size biasses):
                [segment_appear_flight, segment_disappear_flight, ~] = map_detect_segments_appear_disappear(fields_per_win_small,FE,{'none'});
                count_start_end_small_fields = [count_start_end_small_fields ; [segment_appear_flight' segment_disappear_flight']];
            end
            if sum(size_per_field >= field_size_th)>0
                [segment_appear_flight, segment_disappear_flight, ~] = map_detect_segments_appear_disappear(fields_per_win_large,FE);
                count_seg_change_large_fields = [count_seg_change_large_fields ; [segment_appear_flight' segment_disappear_flight']];
                % only field start, end (to avoid field size biasses):
                [segment_appear_flight, segment_disappear_flight, ~] = map_detect_segments_appear_disappear(fields_per_win_large,FE,{'none'});
                count_start_end_large_fields = [count_start_end_large_fields ; [segment_appear_flight' segment_disappear_flight']];
            end
        end
    end
end

%% probability - devided for days:
axes(panel_F);
cla
hold on
text(-0.25,1.67, 'F', 'Units','normalized','FontWeight','bold');
text(0,1.58, {'Population: Within-day dynamics is most prominent',... 
    'during the first 2 days, but occurs also later'},...
    'Units','normalized','FontWeight','bold');

P0_appear1 = mean(count_seg_change_part1(:,1)==0);
P0_disappear1 = mean(count_seg_change_part1(:,2)==0);
se_appear1 = sqrt(P0_appear1*(1-P0_appear1)/length(count_seg_change_part1(:,1)));
se_disappear1 = sqrt(P0_disappear1*(1-P0_disappear1)/length(count_seg_change_part1(:,1)));

P0_appear2 = mean(count_seg_change_part2(:,1)==0);
P0_disappear2 = mean(count_seg_change_part2(:,2)==0);
se_appear2 = sqrt(P0_appear2*(1-P0_appear2)/length(count_seg_change_part2(:,1)));
se_disappear2 = sqrt(P0_disappear2*(1-P0_disappear2)/length(count_seg_change_part2(:,1)));

P0_appear3 = mean(count_seg_change_part3(:,1)==0);
P0_disappear3 = mean(count_seg_change_part3(:,2)==0);
se_appear3 = sqrt(P0_appear3*(1-P0_appear3)/length(count_seg_change_part3(:,1)));
se_disappear3 = sqrt(P0_disappear3*(1-P0_disappear3)/length(count_seg_change_part3(:,1)));

P0_appear4 = mean(count_seg_change_part4(:,1)==0);
P0_disappear4 = mean(count_seg_change_part4(:,2)==0);
se_appear4 = sqrt(P0_appear4*(1-P0_appear4)/length(count_seg_change_part4(:,1)));
se_disappear4 = sqrt(P0_disappear4*(1-P0_disappear4)/length(count_seg_change_part4(:,1)));

% P0_appear5 = mean(count_seg_change_part5(:,1)==0);
% P0_disappear5 = mean(count_seg_change_part5(:,2)==0);
% se_appear5 = sqrt(P0_appear5*(1-P0_appear5)/length(count_seg_change_part5(:,1)));
% se_disappear5 = sqrt(P0_disappear5*(1-P0_disappear5)/length(count_seg_change_part5(:,1)));

x = [1,2, 4,5, 7,8 10,11];%, 13,14];
y = [1-P0_appear1, 1-P0_disappear1, ...
    1-P0_appear2, 1-P0_disappear2, ...
    1-P0_appear3, 1-P0_disappear3 ...
    1-P0_appear4, 1-P0_disappear4];%, ...
%     1-P0_appear5, 1-P0_disappear5];
hb = bar(x,y);
hb.FaceColor = 'flat';
hb.CData = [[0 0 0]; [1 1 1]; [0 0 0]; [1 1 1]; [0 0 0]; [1 1 1];...
    [0 0 0]; [1 1 1]];%; [0 0 0]; [1 1 1]];
he=errorbar(x,y,[se_appear1 se_disappear1 ...
    se_appear2 se_disappear2 ...
    se_appear3 se_disappear3 ...
    se_appear4 se_disappear4]);% ...
%     se_appear5 se_disappear5]);
he.CapSize = 2;
he.LineWidth = 1;
he.LineStyle = 'none';
he.Color = 'k';
h=gca;
h.XTick = [1.5, 4.5, 7.5 10.5];
% h.XTickLabels = {'Appear','Disappear','Appear','Disappear','Appear','Disappear',...
%     'Appear','Disappear','Appear','Disappear','Appear','Disappear'};
h.XTickLabels = {'Day 1-2','Days 3-4','Days 5-6','Days \geq 7'};
h.XTickLabelRotation = 45;
ylimits = [0,0.2];
h.YLim = ylimits;
h.YTick = ylimits;
h.TickLength = [0.01 0.01];
h.XAxis.TickLength = [0 0];
h.XRuler.TickLabelGapMultiplier = -0.1;
h.YRuler.TickLabelGapMultiplier = 0.001;
ylabel('Probability',  'Units','normalized','Position',[-0.07 0.5]);

% test significance:
% n = length(count_seg_change_part1);
% p_pooled_1 = ((1-P0_appear1) + (1-P0_disappear1))*n/(2*n);
% se_pooled = sqrt(p_pooled_1*(1-p_pooled_1)*(2/n));
% z_score = ((1-P0_appear1) - (1-P0_disappear1))/se_pooled;

% parametric significance test:
n1 = length(count_seg_change_part1);
n2 = length(count_seg_change_part2);
p_pooled = ((1-P0_appear1)*n1 + (1-P0_appear2)*n2)/(n1 + n2);
se_pooled = sqrt(p_pooled*(1-p_pooled)*((1/n1)+(1/n2)));
z_score_1_2_app = ((1-P0_appear1) - (1-P0_appear2))/se_pooled;
p_value_1_2_app = normcdf(-1*z_score_1_2_app)*2;
p_pooled = ((1-P0_disappear1)*n1 + (1-P0_disappear2)*n2)/(n1 + n2);
se_pooled = sqrt(p_pooled*(1-p_pooled)*((1/n1)+(1/n2)));
z_score_1_2_dis = ((1-P0_disappear1) - (1-P0_disappear2))/se_pooled;
p_value_1_2_dis = normcdf(-1*z_score_1_2_dis)*2;
p_value_1_2 = max(p_value_1_2_app, p_value_1_2_dis)

n3 = length(count_seg_change_part3);
p_pooled = ((1-P0_appear1)*n1 + (1-P0_appear3)*n3)/(n1 + n3);
se_pooled = sqrt(p_pooled*(1-p_pooled)*((1/n1)+(1/n3)));
z_score_1_3_app = ((1-P0_appear1) - (1-P0_appear3))/se_pooled;
p_value_1_3_app = normcdf(-1*z_score_1_3_app)*2;
p_pooled = ((1-P0_disappear1)*n1 + (1-P0_disappear3)*n3)/(n1 + n3);
se_pooled = sqrt(p_pooled*(1-p_pooled)*((1/n1)+(1/n3)));
z_score_1_3_dis = ((1-P0_disappear1) - (1-P0_disappear3))/se_pooled;
p_value_1_3_dis = normcdf(-1*z_score_1_3_dis)*2;
p_value_1_3 = max(p_value_1_3_app, p_value_1_3_dis)

n4 = length(count_seg_change_part4);
p_pooled = ((1-P0_appear1)*n1 + (1-P0_appear4)*n4)/(n1 + n4);
se_pooled = sqrt(p_pooled*(1-p_pooled)*((1/n1)+(1/n4)));
z_score_1_4_app = ((1-P0_appear1) - (1-P0_appear4))/se_pooled;
p_value_1_4_app = normcdf(-1*z_score_1_4_app)*2;
p_pooled = ((1-P0_disappear1)*n1 + (1-P0_disappear4)*n4)/(n1 + n4);
se_pooled = sqrt(p_pooled*(1-p_pooled)*((1/n1)+(1/n4)));
z_score_1_4_dis = ((1-P0_disappear1) - (1-P0_disappear4))/se_pooled;
p_value_1_4_dis = normcdf(-1*z_score_1_4_dis)*2;
p_value_1_4 = max(p_value_1_4_app, p_value_1_4_dis)

% n5 = length(count_seg_change_part5);
% p_pooled = ((1-P0_appear1)*n1 + (1-P0_appear5)*n5)/(n1 + n5);
% se_pooled = sqrt(p_pooled*(1-p_pooled)*((1/n1)+(1/n5)));
% z_score_1_5_app = ((1-P0_appear1) - (1-P0_appear5))/se_pooled;
% p_value_1_5_app = normcdf(-1*z_score_1_5_app)*2;
% p_pooled = ((1-P0_disappear1)*n1 + (1-P0_disappear5)*n5)/(n1 + n5);
% se_pooled = sqrt(p_pooled*(1-p_pooled)*((1/n1)+(1/n5)));
% z_score_1_5_dis = ((1-P0_disappear1) - (1-P0_disappear5))/se_pooled;
% p_value_1_5_dis = normcdf(-1*z_score_1_5_dis)*2;
% p_value_1_5 = max(p_value_1_5_app, p_value_1_5_dis)

% non-parametric significance test:
% [p_value_app_1_2, p_value_disapp_1_2] = non_parametric_compare(count_seg_change_part1, count_map_part1, count_seg_change_part2, count_map_part2)
% [p_value_app_1_3, p_value_disapp_1_3] = non_parametric_compare(count_seg_change_part1, count_map_part1, count_seg_change_part3, count_map_part3)
% [p_value_app_1_4, p_value_disapp_1_4] = non_parametric_compare(count_seg_change_part1, count_map_part1, count_seg_change_part4, count_map_part4)
% [p_value_app_1_2, p_value_disapp_1_2] = non_parametric_compare(count_seg_change_part1, flight_ID_part1, count_seg_change_part2, flight_ID_part2)
% [p_value_app_1_3, p_value_disapp_1_3] = non_parametric_compare(count_seg_change_part1, flight_ID_part1, count_seg_change_part3, flight_ID_part3)
% [p_value_app_1_4, p_value_disapp_1_4] = non_parametric_compare(count_seg_change_part1, flight_ID_part1, count_seg_change_part4, flight_ID_part4)


set(gca,'clipping','off')
line_w = 0.8;
line([x(1), x(2)], ylimits(2)*0.95*[1,1], 'color', 'k', 'linewidth',line_w);
line([x(3), x(4)], ylimits(2)*0.55*[1,1], 'color', 'k', 'linewidth',line_w);
line([x(5), x(6)], ylimits(2)*0.55*[1,1], 'color', 'k', 'linewidth',line_w);
line([x(7), x(8)], ylimits(2)*0.55*[1,1], 'color', 'k', 'linewidth',line_w);

line([mean([x(1),x(2)]), mean([x(3),x(4)])], ylimits(2)*1*[1,1], 'color', 'k', 'linewidth',line_w);
line(mean([x(3),x(4)])*[1,1], ylimits(2)*[0.55,1], 'color', 'k', 'linewidth',line_w);
text(mean([x(1),x(4)]), ylimits(2)*1.02, '***', 'HorizontalAlignment','center','FontSize',12);

line([mean([x(1),x(2)]), mean([x(5),x(6)])], ylimits(2)*1.12*[1,1], 'color', 'k', 'linewidth',line_w);
line(mean([x(5),x(6)])*[1,1], ylimits(2)*[0.55,1.12], 'color', 'k', 'linewidth',line_w);
text(mean([x(1),x(6)]), ylimits(2)*1.14, '****', 'HorizontalAlignment','center','FontSize',12);

line([mean([x(1),x(2)]), mean([x(7),x(8)])], ylimits(2)*1.24*[1,1], 'color', 'k', 'linewidth',line_w);
line(mean([x(7),x(8)])*[1,1], ylimits(2)*[0.55,1.24], 'color', 'k', 'linewidth',line_w);
text(mean([x(1),x(8)]), ylimits(2)*1.26, '*****', 'HorizontalAlignment','center','FontSize',12);

line(mean([x(1),x(2)])*[1,1], ylimits(2)*[0.95,1.24], 'color', 'k', 'linewidth',line_w);

%% legend:
axes(panel_F_legend);
cla
box off
hax = gca;
h_app = annotation('rectangle',[0 0.7 0.25 0.25],'FaceColor','k');
h_app.Parent = hax;
h_dis = annotation('rectangle',[0 0.3 0.25 0.25],'Color','k');
h_dis.Parent = hax;
text(0.35,0.85,'Appearance', 'HorizontalAlignment','left','FontSize',8);
text(0.35,0.45,'Disappearance', 'HorizontalAlignment','left','FontSize',8);
hax.Visible = 'off';

%% P-value for appeared vs disappeared poolled for all days:
count_seg_change_all_parts = [count_seg_change_part1; count_seg_change_part2;...
    count_seg_change_part3; count_seg_change_part4];
count_map_all_parts = [count_map_part1; count_map_part2;...
    count_map_part3; count_map_part4];
P0_appear_all_parts = mean(count_seg_change_all_parts(:,1)==0);
P0_disappear_all_parts = mean(count_seg_change_all_parts(:,2)==0);
% se_appear_all_parts = sqrt(P0_appear_all_parts*(1-P0_appear_all_parts)/length(count_seg_change_all_parts(:,1)));
% se_disappear_all_parts = sqrt(P0_disappear_all_parts*(1-P0_disappear_all_parts)/length(count_seg_change_all_parts(:,1)));

% non-parametric significance test:
% p_value_shuffle = shuffle_appear_disappear(count_seg_change_all_parts, count_map_all_parts, 10000)
% p_value_shuffle_flight = shuffle_appear_disappear(count_seg_change_all_parts, flight_ID_app, 10000)

% parametric (difference in proportions z-test) significance test:
n = length(count_seg_change_all_parts(:,1));
p_pooled = ((1-P0_appear_all_parts)*n + (1-P0_disappear_all_parts)*n)/(n + n);
se_pooled = sqrt(p_pooled*(1-p_pooled)*((1/n)+(1/n)));
z_score_app_vs_dis = ((1-P0_appear_all_parts) - (1-P0_disappear_all_parts))/se_pooled;
if z_score_app_vs_dis>=0
    p_value_app_vs_dis = normcdf(-1*z_score_app_vs_dis)*2
else
    p_value_app_vs_dis = normcdf(z_score_app_vs_dis)*2
end


%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
% fig_name_out = fullfile(res_dir, sprintf('%s__bin%d_example%d',fig_name_str,day_bin_size,dynamic_examples_option));
% fig_name_out = fullfile(res_dir, sprintf('%s__corr_%s_%d_paramset_%d',fig_name_str,corr_type,field_speed_opt,prm.parmaset));
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');


%%
function p_value = shuffle_appear_disappear(count_seg_change, count_map, n_shuffle)

rng(0);
p0_app_shuffle = zeros([1,n_shuffle]);
p0_dis_shuffle = zeros([1,n_shuffle]);

unique_ID = unique(count_map);
for ii_shuffle = 1:n_shuffle
    flip_rand = rand([1,length(unique_ID)]);
    ID_flip = unique_ID(flip_rand>0.5);
    count_seg_change_new = count_seg_change;
    
    for ii_flip = 1:length(ID_flip)
        count_seg_change_new(count_map == ID_flip(ii_flip),1) = ...
            count_seg_change(count_map == ID_flip(ii_flip),2);
        count_seg_change_new(count_map == ID_flip(ii_flip),2) = ...
            count_seg_change(count_map == ID_flip(ii_flip),1);
    end
    p0_app_shuffle(ii_shuffle) = mean(count_seg_change_new(:,1)==0);
    p0_dis_shuffle(ii_shuffle) = mean(count_seg_change_new(:,2)==0);
end

P0_appear = mean(count_seg_change(:,1)==0);
P0_disappear = mean(count_seg_change(:,2)==0);
diff_value = (1-P0_appear) - (1-P0_disappear);
diff_value_shuffle = (1-p0_app_shuffle) - (1-p0_dis_shuffle);
if diff_value<=0
    p_small_than = sum(diff_value_shuffle <= diff_value)/n_shuffle;
    p_large_than = sum(diff_value_shuffle >= (-1*diff_value))/n_shuffle;
else
    p_small_than = sum(diff_value_shuffle <= (-1*diff_value))/n_shuffle;
    p_large_than = sum(diff_value_shuffle >= diff_value)/n_shuffle;
end
p_value = p_small_than+p_large_than;

end

%% 
function [p_value_app, p_value_disapp] = non_parametric_compare(count_seg_change1, count_map1, count_seg_change2, count_map2)
rng(0);
count_seg_change_pooled = [count_seg_change1; count_seg_change2];
count_map_pooled = [count_map1; count_map2];
count_map_pooled_unique = unique(count_map_pooled);
n_maps = length(count_map_pooled_unique);
chunk1 = length(unique(count_map1));
chunk2 = length(unique(count_map2));

p_chunk1 = [];
p_chunk2 = [];
n_shuffle = 10000;
for ii_rand = 1:n_shuffle
    rand_maps_IX = count_map_pooled_unique(randperm(n_maps));
    rand_maps_chunk1 = [];
    rand_maps_chunk2 = [];
    for ii_map = 1:chunk1
        IX = find(count_map_pooled == rand_maps_IX(ii_map));
        rand_maps_chunk1 = [rand_maps_chunk1; count_seg_change_pooled(IX,:)];
    end
    for ii_map = chunk1 + (1:chunk2)
        IX = find(count_map_pooled == rand_maps_IX(ii_map));
        rand_maps_chunk2 = [rand_maps_chunk2; count_seg_change_pooled(IX,:)];
    end
    
    p_chunk1(ii_rand,1) = 1 - mean(rand_maps_chunk1(:,1)==0);
    p_chunk1(ii_rand,2) = 1 - mean(rand_maps_chunk1(:,2)==0);
    p_chunk2(ii_rand,1) = 1 - mean(rand_maps_chunk2(:,1)==0);
    p_chunk2(ii_rand,2) = 1 - mean(rand_maps_chunk2(:,2)==0);
end

p_app_diff_real = (1 - mean(count_seg_change1(:,1)==0)) - (1 - mean(count_seg_change2(:,1)==0));
p_disapp_diff_real = (1 - mean(count_seg_change1(:,2)==0)) - (1 - mean(count_seg_change2(:,2)==0));

p_app_diff_shuffle = p_chunk1(:,1) - p_chunk2(:,1);
p_disapp_diff_shuffle = p_chunk1(:,2) - p_chunk2(:,2);

if p_app_diff_real<=0
    p_small_than = sum(p_app_diff_shuffle <= p_app_diff_real)/n_shuffle;
    p_large_than = sum(p_app_diff_shuffle >= (-1*p_app_diff_real))/n_shuffle;
else
    p_small_than = sum(p_app_diff_shuffle <= (-1*p_app_diff_real))/n_shuffle;
    p_large_than = sum(p_app_diff_shuffle >= p_app_diff_real)/n_shuffle;
end
p_value_app = p_small_than+p_large_than;

if p_disapp_diff_real<=0
    p_small_than = sum(p_disapp_diff_shuffle <= p_disapp_diff_real)/n_shuffle;
    p_large_than = sum(p_disapp_diff_shuffle >= (-1*p_disapp_diff_real))/n_shuffle;
else
    p_small_than = sum(p_disapp_diff_shuffle <= (-1*p_disapp_diff_real))/n_shuffle;
    p_large_than = sum(p_disapp_diff_shuffle >= p_disapp_diff_real)/n_shuffle;
end
p_value_disapp = p_small_than+p_large_than;

end

