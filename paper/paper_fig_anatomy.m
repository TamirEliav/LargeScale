%% Large Scale - fig. supp - functional anatomy

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S1';
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
panel_CDE_size = [4.5 3];
panel_A(1) = axes('position', [ 0.6 12.5 8 16]);
panel_B(1) = axes('position', [ 8.6 19 4 5]);
panel_B(2) = axes('position', [13.6 21 6 5]);
panel_B(3) = axes('position', [13.6 16.5 6 5]);
panel_C(1) = axes('position', [3 12   panel_CDE_size]);
% panel_C(2) = axes('position', [9 12   panel_CDE_size]);
% panel_D(1) = axes('position', [3  7.5 panel_CDE_size]);
% panel_D(2) = axes('position', [9  7.5 panel_CDE_size]);
% panel_E(1) = axes('position', [3  3   panel_CDE_size]);
% panel_E(2) = axes('position', [9  3   panel_CDE_size]);
panel_legend = axes('position', [8.5 23.8 1 1]);


%% bat ID-num map
bat_ID_num_map = containers.Map([34 79 148 2289 9861],...
                                 1:5);
                             
%% prepare 3D data for panel A
% load('L:\TTs_position\plot_3D_surf_CA1\atlas_coordinates_Tamir_Nov2019.mat')
load('L:\TTs_position\plot_3D_surf_CA1\atlas_coordinates_Tamir_Sep2020.mat')
load('L:\TTs_position\plot_3D_surf_CA1\CA1_3D_struct.mat')

X = CA1_3D_struct.x;
Y = CA1_3D_struct.y;
Z = CA1_3D_struct.z;

TT_x = [atlas_coordinates_Tamir(:).rel_AP_position]/1000; %mm
TT_y = [atlas_coordinates_Tamir(:).ML_position]/1000;
TT_z = -[atlas_coordinates_Tamir(:).DV_position]/1000;

TT_xyz = [TT_x' TT_y' TT_z'];
ALL_fitted_points = [X(:) Y(:) Z(:)];

TT_xyz_tensor = repmat(TT_xyz,1,1,length(ALL_fitted_points(:,1)));
TT_xyz_tensor = shiftdim(TT_xyz_tensor,1); %3xNxN_tt

fitted_tensor = repmat(ALL_fitted_points',1,1,length(TT_x)); %3xNxN_tt

%find closest point on fitted_tensor surface to each of the TTs positions:
norma_TT_to_curve_XY = shiftdim(sqrt(sum((fitted_tensor - TT_xyz_tensor).^2,1))); %NxN_tt
[~, ind] = min(norma_TT_to_curve_XY,[],1);
TT_z_new = ALL_fitted_points(ind,3)+0.050;%so the points will be above the surface

%% panel A - 3D plot of hippocampus CA1 region + TT tracks
axes(panel_A);
cla
hold on
text(0.15,0.955, 'A', 'Units','normalized','FontWeight','bold');

long_axis = linspace(100,0,size(X,1))';
s = surf(X,Y,Z, repmat(long_axis,1,size(X,2)));
% s = surf(X/1000,Y/1000,Z/1000); %DV colormap
set(gca, 'xtick',0:1:4,'ytick', 0:1:6, 'ztick',-10:1:-4);
s.FaceAlpha = 0.6;
s.EdgeColor = 'none';
light(gca,'position',[-0.5,-0.05,0.8],'style','infinite');
% s.FaceLighting = 'gouraud';
xlabel('AP (\mum)'); ylabel('ML (\mum)'); zlabel('DV (\mum)');
hold on;
axis equal
% TT tracks
% scatter3(TT_x, TT_y, TT_z_new, 10, 'k', 'filled');
plot3(repmat(TT_x,2,1), repmat(TT_y,2,1), [TT_z_new'; TT_z_new'+5], '-k', 'LineWidth',1);

grid on
grid off
axis off
view([-60 16])
ha = gca;
ha.XLim = [-0.5 4.5];
ha.YLim = [0.5 7.25];
ha.ZLim = [-11 -2];

% add x/y/z scale bar
scale_bar_center = [2,3,-7];
scale_bar_length = [-1 1];
line(scale_bar_center(1)+scale_bar_length,scale_bar_center(2)+[0 0],scale_bar_center(3)+[0 0],'Color','k','LineWidth',1.5);
line(scale_bar_center(1)+[0 0],scale_bar_center(2)+scale_bar_length,scale_bar_center(3)+[0 0],'Color','k','LineWidth',1.5);
line(scale_bar_center(1)+[0 0],scale_bar_center(2)+[0 0],scale_bar_center(3)+scale_bar_length,'Color','k','LineWidth',1.5);

text([0.7 3.5],[3 3.05],[-7 -7], {'A';'P'}, 'FontWeight','bold','FontSize',10,'HorizontalAlignment','center', 'VerticalAlignment','middle');
text([2 2],[1.7 4.25],[-7 -7], {'M';'L'}, 'FontWeight','bold','FontSize',10,'HorizontalAlignment','center', 'VerticalAlignment','middle');
text([2 2],[3 3],[-8.23 -5.7], {'V';'D'}, 'FontWeight','bold','FontSize',10,'HorizontalAlignment','center', 'VerticalAlignment','middle');

text(1,1.5,-5,    {'Septal';'pole'}, 'FontSize',8,'HorizontalAlignment','left');
text(1,2.3,-10.3, {'Temporal';'pole'}, 'FontSize',8,'HorizontalAlignment','left');
text(1,4,-1.4,       {'Recording';'locations'}, 'FontSize',8,'HorizontalAlignment','left');

% add colorbar
POSf = ds2nfu(panel_A(1), [0 0 0 0]);
POSf([1 2]) = POSf([1 2]) + [0.05 0.1];
POSf([3 4]) = [0.01 0.10];
% colorbar_ax_pos = [panel_A(1).Position(1:2)+[panel_A(1).Position(3) 1] 0.7 3];
c=colorbar('location','EastOutside', 'position', POSf);
c.Ticks = [0,100];
c.TickLabels = {'100 ','0'};
c.Label.String = 'Longitudinal axis (%)';
c.Label.FontSize = 8;

colormap cool


%% load population data
% =========================================================================
prm = PARAMS_GetAll();
cells_t = DS_get_cells_summary();
bats = [79,148,2289,9861];
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
cells = cellfun(@(c)(cell_load_data(c,'details','stats','meanFR','stats','inclusion','signif','fields','FR_map')), cells_ID, 'UniformOutput',0);
cells = [cells{:}];
whos cells

% load('L:\TTs_position\atlas_coordinates_Tamir_sep2019.mat');
% load('L:\TTs_position\Prox2dist_curve.mat');
% load('L:\TTs_position\atlas_coordinates_Tamir_Nov2019.mat');
load('L:\TTs_position\atlas_coordinates_Tamir_Sep2020.mat');
atlas_coordinates = struct2table(atlas_coordinates_Tamir);

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

%% panel B - plot TT positions (only with signif cells!)
axes(panel_B(1));
cla
hold on
text(-0.25,1.24, 'B', 'Units','normalized','FontWeight','bold');

axis ij
axis xy

cells_signif = pop_details(signif_at_least_one_dir_IX);
cells_signif_color = {pop_bat_color{signif_at_least_one_dir_IX}};
[C,IA,IC] = unique({cells_signif.TT_ID});
TT_pos_PD = cat(1,cells_signif(IA).TT_pos_proximodistal_prc);
TT_pos_Long = cat(1,cells_signif(IA).TT_pos_longitudinal_prc);
TT_pos_color = cat(1,cells_signif_color{IA});
scatter(TT_pos_PD*100, TT_pos_Long*100, 5, TT_pos_color);

xlim([0 100])
ylim([15 22]);
xlabel('Proximo-distal axis (%)');
ylabel('Longitudinal axis (%)');

%% draw lines from points to histology slides
slide_examples_TT = {
'bat_0148_TT1';
'bat_9861_TT3';
};
slide_examples_TT_pos = [];
for ii_TT = 1:length(slide_examples_TT)
    TT_IX = find(contains({atlas_coordinates_Tamir.name}, slide_examples_TT{ii_TT}));
    TT_pos_details = atlas_coordinates_Tamir(TT_IX);
    slide_examples_TT_pos(ii_TT,:) = [TT_pos_details.proximo_distal_TT_precentage TT_pos_details.longitudinal_TT_precentage];
end
slide_examples_TT_pos = slide_examples_TT_pos .* 100; % change to percentage

xa = slide_examples_TT_pos(:,1);
ya = slide_examples_TT_pos(:,2);
axes(panel_B(1));
[xaf,yaf] = ds2nfu(xa,ya);
xaf = xaf + 0.003;
yaf = yaf + [1; -1].*0.0015;
% annotation('line', xaf, yaf)
% plot(xa,ya)
annotation('line', [xaf(1) 0.64], [yaf(1) 0.90]);
annotation('line', [xaf(2) 0.64], [yaf(2) 0.74]);

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
    text(1+0.25, ii_bat, sprintf('Bat %d',bat_ID_num_map(bat_num)), 'FontSize',7,'HorizontalAlignment','Left');
end
xlim([0 1]);
ylim([1 length(bats)]);
hax=gca;
hax.YDir = 'reverse';
set(gca,'Visible','off');


%% histology slice example
axes(panel_B(2));
cla
hold on
histology_slice_example_file = 'L:\resources\Histology\processed\Tamir_bat148_Sec16a_X4_after_WB.tif';
image = imread(histology_slice_example_file);
imshow(image);
% add arrow to indicate the track/lesion
harr = annotation('arrow', [0 1],[0 1]);
hax = gca;
harr.Parent = hax;
harr.X = 870 + [0 6];
harr.Y = 350 + [9 0];
harr.LineStyle = 'none';
harr.HeadLength = 7;
harr.HeadWidth = 7;
harr.Color = [1 0 0];
xlim([324 1795]);
ylim([ 67 1116]);
hax.DataAspectRatio(:) = 1;
% add CA1 border lines
plot([1588 1484],[347 539], 'k-', 'LineWidth',0.7); % border correction by nachum
% plot([686 763],[480 550], 'k-', 'LineWidth',0.7);
plot([590 760],[525 570], 'k-', 'LineWidth',0.7); % border correction by nachum
text(1000,222,'CA1','FontSize',8,'FontWeight','bold')
% add scale bar
scale_bar_mm = 1;
pixel2mm = 2.325e-3;
scale_bar_pixels = scale_bar_mm / pixel2mm ;
xa = hax.XLim(1) + [0 scale_bar_pixels];
ya = hax.YLim([1 1]);
[xaf,yaf] = ds2nfu(xa,ya);
% yaf = yaf + 0.016;
hl=annotation('line',xaf,yaf);
hl.LineWidth = 2;
hl.Color = 'k';
% h=annotation('textbox',[xaf(2)+0.001 yaf(1) 0 0], 'String',sprintf('%dmm',scale_bar_mm), 'FitBoxToText','on', 'LineStyle','none');
text(550, 1193, sprintf('%dmm',scale_bar_mm), 'HorizontalAlignment','center', 'VerticalAlignment','middle','FontSize',8);

axes(panel_B(3));
cla
hold on
histology_slice_example_file = 'L:\resources\Histology\processed\Tamir_bat9861_Sec15a_X4_after_WB.tif';
image = imread(histology_slice_example_file);
imshow(image);
% add arrow to indicate the track/lesion
harr = annotation('arrow', [0 1],[0 1]);
hax = gca;
harr.Parent = hax;
harr.X = 920 + [0 11];
harr.Y = 260 + [9 0];
harr.LineStyle = 'none';
harr.HeadLength = 7;
harr.HeadWidth = 7;
harr.Color = [1 0 0];
xlim([80 1934])
ylim([30 1188])
hax.DataAspectRatio(:) = 1;
% add CA1 border lines
% plot([1640 1543],[234 430], 'k-', 'LineWidth',0.7);
plot([1635 1505],[192 437], 'k-', 'LineWidth',0.7); % border correction by nachum
% plot([476 557],[457 580], 'k-', 'LineWidth',0.7);
% plot([370 525],[535 592], 'k-', 'LineWidth',0.7); % border correction by nachum
plot([336 550],[524 602], 'k-', 'LineWidth',0.7); % border correction by nachum - another correction
text(1000,130,'CA1','FontSize',8,'FontWeight','bold')
% add scale bar
scale_bar_mm = 1;
pixel2mm = 2.325e-3;
scale_bar_pixels = scale_bar_mm / pixel2mm ;
xa = hax.XLim(1) + [0 scale_bar_pixels];
ya = hax.YLim([1 1]);
[xaf,yaf] = ds2nfu(xa,ya);
% yaf = yaf + 0.025;
hl=annotation('line',xaf,yaf);
hl.LineWidth = 2;
hl.Color = 'k';
% annotation('textbox',[xaf(2)+0.001 yaf(1) 0.1 0.1], 'String',sprintf('%dmm',scale_bar_mm), 'FitBoxToText','on', 'LineStyle','none');
text(300, 1290, sprintf('%dmm',scale_bar_mm), 'HorizontalAlignment','center', 'VerticalAlignment','middle','FontSize',8);

%%
TT_pos_labels = {'Proximo-distal axis (%)';'Longitudinal axis (%)'};
TT_pos_limits = [0 100; 10 20];
scatter_jitter_std = [1 0.10];
%%
ii_pos_TT_opt=1;
axes(panel_C(ii_pos_TT_opt));
cla('reset')
hold on
text(-0.3,1.1, 'C', 'Units','normalized','FontWeight','bold');
cells_signif = pop_details(signif_at_least_one_dir_IX);
cells_signif_color = {pop_bat_color{signif_at_least_one_dir_IX}};
stats_signif = [pop_stats(signif_at_least_one_dir_IX).all];
x = cat(1,cells_signif.TT_pos_prc) .* 100;
y = [stats_signif.field_ratio_LS];
TT_pos_color = cat(1,cells_signif_color{:});
rng(0);
h=scatter(x(:,ii_pos_TT_opt)+scatter_jitter_std(ii_pos_TT_opt)*randn(size(x,1),1), y, 5, 'k');

[r,pval_r]     = corr(x(:,ii_pos_TT_opt),y','rows','pairwise','type','Pearson');
[rho,pval_rho] = corr(x(:,ii_pos_TT_opt),y','rows','pairwise','type','Spearman');
fprintf('field size ratio vs proximo-distal correlation (pearson): r=%.2f p=%.2f df=%d \n',r,pval_r,sum(~isnan(y))-2);
fprintf('field size ratio vs proximo-distal correlation (spearman): rho=%.2f p=%.2f df=%d \n',rho,pval_rho,sum(~isnan(y))-2);
% text(0.05,0.95, {sprintf('r = %.2f',r);sprintf('P = %.2f',pval_r)}, ...
%     'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);
text(0.05,0.95, {sprintf('%s = %.2f','\rho',rho);sprintf('P = %.2f',pval_rho)}, ...
    'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);

hax=gca;
hax.YScale = 'log';
hax.YTick = [1 2 3 5 10 15 20];
xlim(TT_pos_limits(ii_pos_TT_opt,:))
% ylim([0 15]);
xlabel(TT_pos_labels{ii_pos_TT_opt});
ylabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[-0.1 0.5]);


%%
if 0
    
%% panels BCD labels
TT_pos_labels = {'Proximo-distal (%)';'Longitudinal (%)'};
TT_pos_limits = [0 100; 10 20];
scatter_jitter_std = [1 0.10];

%% panel C - plot L/S ratio vs. TT positions (proximo-distal or Longitudinal)
for ii_pos_TT_opt = 1:2
    axes(panel_C(ii_pos_TT_opt));
    cla
    hold on

    cells_signif = pop_details(signif_at_least_one_dir_IX);
    cells_signif_color = {pop_bat_color{signif_at_least_one_dir_IX}};
    stats_signif = [pop_stats(signif_at_least_one_dir_IX).all];
    x = cat(1,cells_signif.TT_pos_prc) .* 100;
    y = [stats_signif.field_ratio_LS];
    TT_pos_color = cat(1,cells_signif_color{:});
    rng(0);
    h=scatter(x(:,ii_pos_TT_opt)+scatter_jitter_std(ii_pos_TT_opt)*randn(size(x,1),1), y, 5, TT_pos_color);
    % h.Marker = 'x';
    % h.SizeData = 20;
    
    [r,pval_r]     = corr(x(:,ii_pos_TT_opt),y','rows','pairwise','type','Pearson');
    [rho,pval_rho] = corr(x(:,ii_pos_TT_opt),y','rows','pairwise','type','Spearman');
    text(0.05,0.95, {sprintf('r = %.2f',r);sprintf('P = %.2f',pval_r)}, ...
        'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);
%     text(1,0.8, {sprintf('%s=%.2f','\rho',rho);sprintf('P=%.2f',pval_rho)}, ...
%         'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    
    xlim(TT_pos_limits(ii_pos_TT_opt,:))
    ylim([0 15]);
    xlabel(TT_pos_labels{ii_pos_TT_opt});
    ylabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[-0.1 0.5]);
end
axes(panel_C(1));
text(-0.25,1.1, 'C', 'Units','normalized','FontWeight','bold');

%% panel D - plot mean Field size vs. TT positions (proximo-distal or Longitudinal)
for ii_pos_TT_opt = 1:2
    axes(panel_D(ii_pos_TT_opt));
    cla
    hold on

    cells_signif = pop_details(signif_at_least_one_dir_IX);
    cells_signif_color = {pop_bat_color{signif_at_least_one_dir_IX}};
    stats_signif = [pop_stats(signif_at_least_one_dir_IX).all];
    x = cat(1,cells_signif.TT_pos_prc) .* 100;
    y = [stats_signif.field_size_mean];
    TT_pos_color = cat(1,cells_signif_color{:});
    rng(0);
    h=scatter(x(:,ii_pos_TT_opt)+scatter_jitter_std(ii_pos_TT_opt)*randn(size(x,1),1), y, 5, TT_pos_color);
    % h.Marker = 'x';
    % h.SizeData = 20;

    [r,pval_r]     = corr(x(:,ii_pos_TT_opt),y','rows','pairwise','type','Pearson');
    [rho,pval_rho] = corr(x(:,ii_pos_TT_opt),y','rows','pairwise','type','Spearman');
    text(0.05,1, {sprintf('r = %.2f',r);sprintf('P = %.2f',pval_r)}, ...
        'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);
%     text(1,0.8, {sprintf('%s=%.2f','\rho',rho);sprintf('P=%.2f',pval_rho)}, ...
%         'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    
    xlim(TT_pos_limits(ii_pos_TT_opt,:))
    ylim([0 20]);
    xlabel(TT_pos_labels{ii_pos_TT_opt});
    ylabel('Average field size (m)','Units','normalized','Position',[-0.1 0.5]);
end
axes(panel_D(1));
text(-0.25,1.1, 'D', 'Units','normalized','FontWeight','bold');


%% panel E - plot number of fields vs. TT positions (proximo-distal or Longitudinal)
for ii_pos_TT_opt = 1:2
    axes(panel_E(ii_pos_TT_opt));
    cla
    hold on

    cells_signif = pop_details(signif_at_least_one_dir_IX);
    cells_signif_color = {pop_bat_color{signif_at_least_one_dir_IX}};
    stats_signif = [pop_stats(signif_at_least_one_dir_IX).all];
    x = cat(1,cells_signif.TT_pos_prc) .* 100;
    y = [stats_signif.field_num];
    TT_pos_color = cat(1,cells_signif_color{:});
    rng(0);
    h=scatter(x(:,ii_pos_TT_opt)+scatter_jitter_std(ii_pos_TT_opt)*randn(size(x,1),1), y, 5, TT_pos_color);
    % h.Marker = 'x';
    % h.SizeData = 20;

    [r,pval_r]     = corr(x(:,ii_pos_TT_opt),y','rows','pairwise','type','Pearson');
    [rho,pval_rho] = corr(x(:,ii_pos_TT_opt),y','rows','pairwise','type','Spearman');
    text(0.05,1, {sprintf('r = %.2f',r);sprintf('P = %.2f',pval_r)}, ...
        'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);
%     text(1,0.8, {sprintf('%s=%.2f','\rho',rho);sprintf('P=%.2f',pval_rho)}, ...
%         'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    
    xlim(TT_pos_limits(ii_pos_TT_opt,:))
    ylim([0 40]);
    xlabel(TT_pos_labels{ii_pos_TT_opt});
    ylabel('No. of fields','Units','normalized','Position',[-0.1 0.5]);
end
axes(panel_E(1));
text(-0.25,1.1, 'E', 'Units','normalized','FontWeight','bold');

end
%%






%%


%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');


