%% Large Scale - supp fig - comparing flight directions

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_supp_compare_directions';
fig_caption_str = 'compare firing patterns between flight directions';
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
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none');
pause(0.2); % workaround to solve matlab automatically changing the axes positions...

% create panels
panels_size = [4 4];
panel_A = axes('position', [ 2 20 panels_size]);
panel_B = axes('position', [ 7 20 panels_size]);
panel_C = axes('position', [12 20 panels_size]);
panel_D = axes('position', [ 2 14 panels_size]);
panel_E = axes('position', [ 7 14 panels_size]);
panel_F = axes('position', [12 14 panels_size]);
panel_G = axes('position', [ 2  8 panels_size]);
panel_H = axes('position', [ 7  8 panels_size]);
panel_I = axes('position', [12  8 panels_size]);
panel_legend_bat_colors = axes('position', [17 20 3 4]);

%% load population data
% =========================================================================
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
whos cells

%% get population stats
pop_stats = cat(1,cells.stats);
pop_stats_dir = cat(1,pop_stats.dir);
pop_signif = cat(1,cells.signif);
pop_signif_TF = arrayfun(@(x)(x.TF), pop_signif);
signif_both_dir_IX = all(pop_signif_TF,2);
signif_at_least_one_dir_IX = any(pop_signif_TF,2);
pop_details = cat(1,cells.details);
pop_bat_number = [pop_details.bat];
pop_bat_color = arrayfun(@(x)(prm.graphics.colors.bats(x)),pop_bat_number,'UniformOutput',0);

%%
% % % pop_data_field_count = nan(2,length(cells));
% % % pop_data_field_size_mean = nan(2,length(cells));
% % % pop_data_field_size_median = nan(2,length(cells));
% % % pop_data_field_size_max = nan(2,length(cells));
% % % pop_data_field_size_min = nan(2,length(cells));
% % % pop_data_ratio_LS = nan(2,length(cells));
% % % pop_data_SI = nan(2,length(cells));
% % % pop_data_sparsity = nan(2,length(cells));
% % % % arrange data
% % % for ii_dir = 1:2
% % %     for ii_cell = 1:length(cells)
% % %         cell = cells(ii_cell);
% % %         if ~cell.signif(ii_dir).TF
% % %             continue;
% % %         end
% % %         fields = cell.fields{ii_dir};
% % %         fields( [fields.in_low_speed_area] ) = []; % TODO: decide if we want those fields near the landing balls
% % %         pop_data_field_count(ii_dir,ii_cell) = length(fields);
% % %         pop_data_field_size_mean(ii_dir,ii_cell) = mean([fields.width_prc]);
% % %         pop_data_field_size_median(ii_dir,ii_cell) = median([fields.width_prc]);
% % %         pop_data_field_size_max(ii_dir,ii_cell) = max([fields.width_prc]);
% % %         pop_data_field_size_min(ii_dir,ii_cell) = min([fields.width_prc]);
% % %         pop_data_field_size_min(ii_dir,ii_cell) = min([fields.width_prc]);
% % %     end
% % %     
% % % end

%% bat color legend panel
axes(panel_legend_bat_colors);
cla
hold on
cmap = prm.graphics.colors.bats;
bats_num = cmap.keys;
bats_num = [bats_num{:}];
bats_colors = cmap.values;
for ii_bat = 1:length(bats_num)
    plot(1, ii_bat, '.', 'Color', bats_colors{ii_bat}, 'MarkerSize', 15);
    text(2, ii_bat, "bat "+bats_num(ii_bat))
end
xlim([0.5 10])
set(gca,'Visible','off')

%% panel A - number of cells signif
axes(panel_A);
cla
hold on

% histogram(randn(1,100))
% histogram(categorical(pop_bat_number));
histogram(categorical(pop_bat_number(signif_at_least_one_dir_IX)));
histogram(categorical(pop_bat_number(signif_both_dir_IX)));

text(-0.1,1.1, 'A', 'Units','normalized','FontWeight','bold');
xlabel('Bat')
ylabel('No. cells')
legend({'signif one dir';'signif both dir'},'Location','northoutside','Position',[0.15 0.91 0.03 0.02])


%% panel B - Fields count
axes(panel_B);
cla
hold on
text(-0.1,1.1, 'B', 'Units','normalized','FontWeight','bold');

x = [pop_stats_dir(signif_both_dir_IX,1).field_num];
y = [pop_stats_dir(signif_both_dir_IX,2).field_num];
c = cat(1,pop_bat_color{signif_both_dir_IX});
plot_jitter_sigma = 0.18;
x = x + plot_jitter_sigma*randn(size(x));
y = y + plot_jitter_sigma*randn(size(y));
h=scatter(x,y,5,c,'filled');

xlabel('dir1')
ylabel('dir2')
title('Fields count')
axis equal
xlim([0 20])
ylim([0 20])

%% panel C - Fields Size (mean)
axes(panel_C);
cla
hold on
text(-0.1,1.1, 'C', 'Units','normalized','FontWeight','bold');

x = [pop_stats_dir(signif_both_dir_IX,1).field_size_mean];
y = [pop_stats_dir(signif_both_dir_IX,2).field_size_mean];
c = cat(1,pop_bat_color{signif_both_dir_IX});
h=scatter(x,y,5,c,'filled');

xlabel('dir1')
ylabel('dir2')
title('Fields size (mean)')
axis equal
xlim([0 max([x y])+1])
ylim([0 max([x y])+1])

%% panel D - Fields Size (median)
axes(panel_D);
cla
hold on
text(-0.1,1.1, 'D', 'Units','normalized','FontWeight','bold');

x = [pop_stats_dir(signif_both_dir_IX,1).field_size_median];
y = [pop_stats_dir(signif_both_dir_IX,2).field_size_median];
c = cat(1,pop_bat_color{signif_both_dir_IX});
h=scatter(x,y,5,c,'filled');

xlabel('dir1')
ylabel('dir2')
title('Fields size (median)')
axis equal
xlim([0 max([x y])+1])
ylim([0 max([x y])+1])

%% panel E - Fields Size (Largest)
axes(panel_E);
cla
hold on
text(-0.1,1.1, 'D', 'Units','normalized','FontWeight','bold');

x = [pop_stats_dir(signif_both_dir_IX,1).field_largest];
y = [pop_stats_dir(signif_both_dir_IX,2).field_largest];
c = cat(1,pop_bat_color{signif_both_dir_IX});
h=scatter(x,y,5,c,'filled');

xlabel('dir1')
ylabel('dir2')
title('Fields size (largest)')
axis equal
xlim([0 max([x y])+1])
ylim([0 max([x y])+1])

%% panel F - Fields Size (mean)
axes(panel_F);
cla
hold on
text(-0.1,1.1, 'F', 'Units','normalized','FontWeight','bold');

x = [pop_stats_dir(signif_both_dir_IX,1).field_smallest];
y = [pop_stats_dir(signif_both_dir_IX,2).field_smallest];
c = cat(1,pop_bat_color{signif_both_dir_IX});
h=scatter(x,y,5,c,'filled');

xlabel('dir1')
ylabel('dir2')
title('Fields size (smallest)')
axis equal
xlim([0 max([x y])+1])
ylim([0 max([x y])+1])

%% panel G - ratio largest/smallest
axes(panel_G);
cla
hold on
text(-0.1,1.1, 'G', 'Units','normalized','FontWeight','bold');

x = [pop_stats_dir(signif_both_dir_IX,1).field_ratio_LS];
y = [pop_stats_dir(signif_both_dir_IX,2).field_ratio_LS];
c = cat(1,pop_bat_color{signif_both_dir_IX});
h=scatter(x,y,5,c,'filled');

xlabel('dir1')
ylabel('dir2')
title('ratio L/S')
axis equal
xlim([1 max([x y])+1])
ylim([1 max([x y])+1])

%% panel H - spatial info
axes(panel_H);
cla
hold on
text(-0.1,1.1, 'H', 'Units','normalized','FontWeight','bold');

x = [pop_stats_dir(signif_both_dir_IX,1).SI_bits_spike];
y = [pop_stats_dir(signif_both_dir_IX,2).SI_bits_spike];
c = cat(1,pop_bat_color{signif_both_dir_IX});
h=scatter(x,y,5,c,'filled');

xlabel('dir1')
ylabel('dir2')
title('spatial info')
axis equal
xlim([0 max([x y])+1])
ylim([0 max([x y])+1])

%% panel I - sparsity
axes(panel_I);
cla
hold on
text(-0.1,1.1, 'I', 'Units','normalized','FontWeight','bold');

x = [pop_stats_dir(signif_both_dir_IX,1).sparsity];
y = [pop_stats_dir(signif_both_dir_IX,2).sparsity];
c = cat(1,pop_bat_color{signif_both_dir_IX});
h=scatter(x,y,5,c,'filled');

xlabel('dir1')
ylabel('dir2')
title('sparsity')
axis equal
xlim([0 1])
ylim([0 1])

%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');








