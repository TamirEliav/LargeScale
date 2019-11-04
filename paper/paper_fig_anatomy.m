%% Large Scale - fig. supp - functional anatomy

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_supp_anatomy';
fig_caption_str = 'Functional anatomy';
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
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');
pause(0.2); % workaround to solve matlab automatically changing the axes positions...

% create panels
panel_A_size = [4 5];
panel_BCD_size = [6 4];
panel_A(1) = axes('position', [ 3 20.5 panel_A_size]);
panel_A(2) = axes('position', [ 8 21.5 4 3]);
panel_B(1) = axes('position', [ 3 15 panel_BCD_size]);
panel_B(2) = axes('position', [11 15 panel_BCD_size]);
panel_C(1) = axes('position', [ 3 10 panel_BCD_size]);
panel_C(2) = axes('position', [11 10 panel_BCD_size]);
panel_D(1) = axes('position', [ 3  5 panel_BCD_size]);
panel_D(2) = axes('position', [11  5 panel_BCD_size]);
panel_legend = axes('position', [16 21.5 1 3]);

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

load('L:\TTs_position\atlas_coordinates_Tamir_sep2019.mat');
load('L:\TTs_position\Prox2dist_curve.mat');

%% get population stats
pop_stats = cat(1,cells.stats);
pop_stats_dir = cat(1,pop_stats.dir);
pop_signif = cat(1,cells.signif);
pop_signif_TF = arrayfun(@(x)(x.TF), pop_signif);
signif_both_dir_IX = all(pop_signif_TF,2);
signif_at_least_one_dir_IX = any(pop_signif_TF,2);
pop_details = cat(1,cells.details);
pop_bat_number = [pop_details.bat];
color_by_bat = 1;
if color_by_bat 
    pop_bat_color = arrayfun(@(x)(prm.graphics.colors.bats(x)),pop_bat_number,'UniformOutput',0);
else
    pop_bat_color = arrayfun(@(x)([0 0 0]),pop_bat_number,'UniformOutput',0);
end

%% legend panel
axes(panel_legend);
cla
hold on
prm = PARAMS_GetAll();
bats_colors = prm.graphics.colors.bats;
for ii_bat = 1:length(bats)
    bat_num = bats(ii_bat);
    c = bats_colors(bat_num);
    scatter(1,ii_bat,9,c);
    text(0.6, ii_bat, sprintf('bat %d',bat_num), 'FontSize',7,'HorizontalAlignment','right');
end
xlim([0 1]);
ylim([1 length(bats)]);
set(gca,'Visible','off');

%% panel A - plot TT positions (only with signif cells!)
axes(panel_A(1));
cla
hold on
text(-0.3,1.1, 'A', 'Units','normalized','FontWeight','bold');

axis ij
plot(Prox2dist_curve(1,:),Prox2dist_curve(2,:), 'k');
% plot([0 Prox2dist_curve(1,1)], Prox2dist_curve(2,[1 1]),'k')
% plot([0 Prox2dist_curve(1,end)], Prox2dist_curve(2,[end end]),'k')

cells_signif = pop_details(signif_at_least_one_dir_IX);
cells_signif_color = {pop_bat_color{signif_at_least_one_dir_IX}};
[C,IA,IC] = unique({cells_signif.TT_ID});
TT_pos = cat(1,cells_signif(IA).TT_pos);
TT_pos_color = cat(1,cells_signif_color{IA});
scatter(TT_pos(:,1),TT_pos(:,2), 5, TT_pos_color)

xlim([0,2000])
ylim([0,4000]);
xlabel('Proximo-distal axis (\mum)');
ylabel('Longitudinal axis (\mum)');

%% histology slice example
axes(panel_A(2));
cla
histology_slice_example_file = 'L:\resources\Histology_bat_148_TT1.jpg';
image = imread(histology_slice_example_file);
imshow(image);

%% panels BCD labels
TT_pos_labels = {'Proximo-distal (\mum)';'Longitudinal (\mum)'};
TT_pos_limits = [0 1500; 2000 3000];

%% panel B - plot L/S ratio vs. TT positions (proximo-distal or Longitudinal)
for ii_pos_TT_opt = 1:2
    axes(panel_B(ii_pos_TT_opt));
    cla
    hold on

    cells_signif = pop_details(signif_at_least_one_dir_IX);
    cells_signif_color = {pop_bat_color{signif_at_least_one_dir_IX}};
    stats_signif = [pop_stats(signif_at_least_one_dir_IX).all];
    x = cat(1,cells_signif.TT_pos);
    y = [stats_signif.field_ratio_LS];
    TT_pos_color = cat(1,cells_signif_color{:});
    rng(0);
    h=scatter(x(:,ii_pos_TT_opt)+10*randn(size(x,1),1), y, 5, TT_pos_color);
    % h.Marker = 'x';
    % h.SizeData = 20;
    
    [r,pval_r]     = corr(x(:,ii_pos_TT_opt),y','rows','pairwise','type','Pearson');
    [rho,pval_rho] = corr(x(:,ii_pos_TT_opt),y','rows','pairwise','type','Spearman');
    text(1,1, {sprintf('r=%.2f',r);sprintf('P=%.2f',pval_r)}, ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    text(1,0.8, {sprintf('%s=%.2f','\rho',rho);sprintf('P=%.2f',pval_rho)}, ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    
    xlim(TT_pos_limits(ii_pos_TT_opt,:))
    ylim([0 15]);
    xlabel(TT_pos_labels{ii_pos_TT_opt});
    ylabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[-0.1 0.5]);
end
axes(panel_B(1));
text(-0.2,1.1, 'B', 'Units','normalized','FontWeight','bold');

%% panel C - plot mean Field size vs. TT positions (proximo-distal or Longitudinal)
for ii_pos_TT_opt = 1:2
    axes(panel_C(ii_pos_TT_opt));
    cla
    hold on

    cells_signif = pop_details(signif_at_least_one_dir_IX);
    cells_signif_color = {pop_bat_color{signif_at_least_one_dir_IX}};
    stats_signif = [pop_stats(signif_at_least_one_dir_IX).all];
    x = cat(1,cells_signif.TT_pos);
    y = [stats_signif.field_size_mean];
    TT_pos_color = cat(1,cells_signif_color{:});
    rng(0);
    h=scatter(x(:,ii_pos_TT_opt)+10*randn(size(x,1),1), y, 5, TT_pos_color);
    % h.Marker = 'x';
    % h.SizeData = 20;

    [r,pval_r]     = corr(x(:,ii_pos_TT_opt),y','rows','pairwise','type','Pearson');
    [rho,pval_rho] = corr(x(:,ii_pos_TT_opt),y','rows','pairwise','type','Spearman');
    text(1,1, {sprintf('r=%.2f',r);sprintf('P=%.2f',pval_r)}, ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    text(1,0.8, {sprintf('%s=%.2f','\rho',rho);sprintf('P=%.2f',pval_rho)}, ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    
    xlim(TT_pos_limits(ii_pos_TT_opt,:))
    ylim([0 20]);
    xlabel(TT_pos_labels{ii_pos_TT_opt});
    ylabel('Averaged field size (m)','Units','normalized','Position',[-0.1 0.5]);
end
axes(panel_C(1));
text(-0.2,1.1, 'C', 'Units','normalized','FontWeight','bold');


%% panel D - plot number of fields vs. TT positions (proximo-distal or Longitudinal)
for ii_pos_TT_opt = 1:2
    axes(panel_D(ii_pos_TT_opt));
    cla
    hold on

    cells_signif = pop_details(signif_at_least_one_dir_IX);
    cells_signif_color = {pop_bat_color{signif_at_least_one_dir_IX}};
    stats_signif = [pop_stats(signif_at_least_one_dir_IX).all];
    x = cat(1,cells_signif.TT_pos);
    y = [stats_signif.field_num];
    TT_pos_color = cat(1,cells_signif_color{:});
    rng(0);
    h=scatter(x(:,ii_pos_TT_opt)+10*randn(size(x,1),1), y, 5, TT_pos_color);
    % h.Marker = 'x';
    % h.SizeData = 20;

    [r,pval_r]     = corr(x(:,ii_pos_TT_opt),y','rows','pairwise','type','Pearson');
    [rho,pval_rho] = corr(x(:,ii_pos_TT_opt),y','rows','pairwise','type','Spearman');
    text(1,1, {sprintf('r=%.2f',r);sprintf('P=%.2f',pval_r)}, ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    text(1,0.8, {sprintf('%s=%.2f','\rho',rho);sprintf('P=%.2f',pval_rho)}, ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    
    xlim(TT_pos_limits(ii_pos_TT_opt,:))
    ylim([0 40]);
    xlabel(TT_pos_labels{ii_pos_TT_opt});
    ylabel('No. of fields','Units','normalized','Position',[-0.1 0.5]);
end
axes(panel_D(1));
text(-0.2,1.1, 'D', 'Units','normalized','FontWeight','bold');


%%






%%


%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');


