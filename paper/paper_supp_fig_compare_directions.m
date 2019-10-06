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

color_by_bat = 0;

% create panels
panels_size = [4 4];
panel_A = axes('position', [ 2   20 panels_size]);
panel_B = axes('position', [ 7.5 20 panels_size]);
panel_C = axes('position', [13   20 panels_size]);
if color_by_bat
    panel_legend_bat_colors = axes('position', [17 20 3 4]);
end

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
if color_by_bat 
    pop_bat_color = arrayfun(@(x)(prm.graphics.colors.bats(x)),pop_bat_number,'UniformOutput',0);
else
    pop_bat_color = arrayfun(@(x)([0 0 0]),pop_bat_number,'UniformOutput',0);
end

%% bat color legend panel
if color_by_bat 
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
end 

%% panel A - number of fields
axes(panel_A);
cla
hold on
text(-0.2,1.13, 'A', 'Units','normalized','FontWeight','bold');

x = [pop_stats_dir(signif_both_dir_IX,1).field_num];
y = [pop_stats_dir(signif_both_dir_IX,2).field_num];
c = cat(1,pop_bat_color{signif_both_dir_IX});
plot_jitter_sigma = 0.16;
rng(0);
h=scatter(  x+plot_jitter_sigma*randn(size(x)),...
            y+plot_jitter_sigma*randn(size(y)),5,c,'filled');

[r,rpval] = corr(x',y','type','Pearson');
% [rho,rhopval] = corr(x',y','type','Spearman');
% signtest_pval = signtest(x,y);
% stats_str = {   sprintf('r=%.2f,P=%g',r,rpval);
%                 sprintf('rho=%.2f,P=%g',rho,rhopval);
%                 sprintf('sign test: P=%g',signtest_pval)};
stats_str = {sprintf('r=%.2g',r); sprintf('P=%.2g',rpval)};
text(0.1,0.95, stats_str, 'units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);

xlabel('Direction 1','Color',prm.graphics.colors.flight_directions{1});
ylabel('Direction 2','Color',prm.graphics.colors.flight_directions{2});
title('No. of fields','units','normalized','Position',[0.5 1.05])
axis equal
xlim([0 20])
ylim([0 20])
h=refline(1,0);
h.Color = 0.5*[1 1 1];

%% panel B - Fields Size (median)
axes(panel_B);
cla
hold on
text(-0.2,1.13, 'B', 'Units','normalized','FontWeight','bold');

x = [pop_stats_dir(signif_both_dir_IX,1).field_size_median];
y = [pop_stats_dir(signif_both_dir_IX,2).field_size_median];
c = cat(1,pop_bat_color{signif_both_dir_IX});
h=scatter(x,y,5,c,'filled');

[r,rpval] = corr(x',y','type','Pearson','rows','complete');
% [rho,rhopval] = corr(x',y','type','Spearman','rows','complete');
% signtest_pval = signtest(x,y);
% stats_str = {   sprintf('r=%.2f,P=%g',r,rpval);
%                 sprintf('rho=%.2f,P=%g',rho,rhopval);
%                 sprintf('sign test: P=%g',signtest_pval)};
stats_str = {sprintf('r=%.2g',r); sprintf('P=%.2g',rpval)};
text(0.1,0.95, stats_str, 'units',...
    'normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);


xlabel('Direction 1','Color',prm.graphics.colors.flight_directions{1});
ylabel('Direction 2','Color',prm.graphics.colors.flight_directions{2});
title('Field size (median per neuron)','units','normalized','Position',[0.5 1.05])
axis equal
xlim([0 max([x y])+1])
ylim([0 max([x y])+1])
h=refline(1,0);
h.Color = 0.5*[1 1 1];

%% panel C - spatial info
axes(panel_C);
cla
hold on
text(-0.2,1.13, 'C', 'Units','normalized','FontWeight','bold');

x = [pop_stats_dir(signif_both_dir_IX,1).SI_bits_spike];
y = [pop_stats_dir(signif_both_dir_IX,2).SI_bits_spike];
c = cat(1,pop_bat_color{signif_both_dir_IX});
h=scatter(x,y,5,c,'filled');

[r,rpval] = corr(x',y','type','Pearson');
% [rho,rhopval] = corr(x',y','type','Spearman');
% signtest_pval = signtest(x,y);
% stats_str = {   sprintf('r=%.2f,P=%g',r,rpval);
%                 sprintf('rho=%.2f,P=%g',rho,rhopval);
%                 sprintf('sign test: P=%g',signtest_pval)};
stats_str = {sprintf('r=%.2g',r); sprintf('P=%.2g',rpval)};
text(0.1,0.95, stats_str, 'units',...
    'normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);

xlabel('Direction 1','Color',prm.graphics.colors.flight_directions{1});
ylabel('Direction 2','Color',prm.graphics.colors.flight_directions{2});
title('Spatial information','units','normalized','Position',[0.5 1.05])
axis equal
xlim([0 max([x y])+1])
ylim([0 max([x y])+1])
h=refline(1,0);
h.Color = 0.5*[1 1 1];


%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');








