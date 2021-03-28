%% Large Scale - supp fig - comparing flight directions

%%
clear 
clc

%% params
color_by_bat = 0;
panel_C_avg = 'median';
% panel_C_avg = 'mean';
% panel_D_scale = 'linear';
panel_D_scale = 'log';

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S10';
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
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');
pause(0.2); % workaround to solve matlab automatically changing the axes positions...

% create panels
panels_size = [4 4];
panel_A = axes('position', [ 5  20 panels_size]);
panel_B = axes('position', [11  20 panels_size]);
panel_C = axes('position', [ 5  14 panels_size]);
panel_D = axes('position', [11  14 panels_size]);
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
cells = cellfun(@(c)(cell_load_data(c,'details','meanFR')), {cells.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];
cells_details = [cells.details];
cells_ID = {cells_details.cell_ID};
meanFR = [cells.meanFR];
cells_ID([meanFR.all]>prm.inclusion.interneuron_FR_thr)=[];
clear cells stats cells_details cells_t meanFR
cells = cellfun(@(c)(cell_load_data(c,'details','stats','meanFR','stats','signif','fields')), cells_ID, 'UniformOutput',0);
% cells = cellfun(@(c)(cell_load_data(c,'details')), cells_ID, 'UniformOutput',0);
cells = [cells{:}];

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


%% panel A - spatial info
axes(panel_A);
cla
hold on
text(-0.23,1.13, 'A', 'Units','normalized','FontWeight','bold');

x = [pop_stats_dir(signif_both_dir_IX,1).SI_bits_spike];
y = [pop_stats_dir(signif_both_dir_IX,2).SI_bits_spike];
c = cat(1,pop_bat_color{signif_both_dir_IX});
h=scatter(x,y,5,c,'filled');

[r,rpval] = corr(x',y','type','Pearson');
[rho,rhopval] = corr(x',y','type','Spearman');
if rhopval==0
    rhopval = realmin;
end
[~,P_ttest] = ttest(x,y);
[P_ranksum,~,STATS_ranksum] = ranksum(x,y);
[P_signtest,~,STATS_signtest] = signtest(x,y);
stats_str = {[sprintf('%s = %.2g','\rho',rho), ',  ' sprintf('P_{%s} = %.2g','\rho',rhopval)];...
              sprintf('P_{wilc} = %.3g',P_ranksum);};
text(0.06,1, stats_str, 'units',...
    'normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);
text(0.95,0.2, "n = "+sum(~any(isnan([x;y]))), 'units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
xlabel('Direction 1','Color',prm.graphics.colors.flight_directions{1},'units','normalized','Position',[0.5 -0.14]);
ylabel('Direction 2','Color',prm.graphics.colors.flight_directions{2},'units','normalized','Position',[-0.14 0.5]);
title('Spatial information (bits/spike)','units','normalized','Position',[0.5 1.05])
axis equal
hax=gca;
hax.XLim = ([0 max([x y])+1]);
hax.YLim = ([0 max([x y])+1]);
hax.XTick = 1:6;
hax.YTick = hax.XTick;
hax.TickLength = [0.02 0.02];
h=refline(1,0);
h.Color = 0.5*[1 1 1];

%% panel B - number of fields
axes(panel_B);
cla
hold on
text(-0.23,1.13, 'B', 'Units','normalized','FontWeight','bold');

x = [pop_stats_dir(signif_both_dir_IX,1).field_num];
y = [pop_stats_dir(signif_both_dir_IX,2).field_num];
c = cat(1,pop_bat_color{signif_both_dir_IX});
plot_jitter_sigma = 0.16;
rng(0);
h=scatter(  x+plot_jitter_sigma*randn(size(x)),...
            y+plot_jitter_sigma*randn(size(y)),5,c,'filled');

[r,rpval] = corr(x',y','type','Pearson');
[rho,rhopval] = corr(x',y','type','Spearman');
if rhopval==0
    rhopval = realmin;
end
[~,P_ttest] = ttest(x,y);
[P_ranksum,~,STATS_ranksum] = ranksum(x,y);
[P_signtest,~,STATS_signtest] = signtest(x,y);
stats_str = {[sprintf('%s = %.2g','\rho',rho), ',  ' sprintf('P_{%s} = %.2g','\rho',rhopval)];...
              sprintf('P_{wilc} = %.2g',P_ranksum);};
text(0.06,1, stats_str, 'units',...
    'normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);
text(0.95,0.2, "n = "+sum(~any(isnan([x;y]))), 'units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);

xlabel('Direction 1','Color',prm.graphics.colors.flight_directions{1},'units','normalized','Position',[0.5 -0.14]);
ylabel('Direction 2','Color',prm.graphics.colors.flight_directions{2},'units','normalized','Position',[-0.14 0.5]);
title('No. of fields','units','normalized','Position',[0.5 1.05])
axis equal
xlim([0 20])
ylim([0 20])
hax=gca;
hax.XScale = 'linear';
hax.YScale = 'linear';
% hax.XScale = 'log';
% hax.YScale = 'log';
hax.XTick = [0:5:20];
hax.YTick = [0:5:20];
hax.TickLength = [0.02 0.02];
hax.YRuler.TickLabelGapOffset = 1;

h=refline(1,0);
h.Color = 0.5*[1 1 1];

%% panel C - Fields Size (median)
axes(panel_C);
cla
hold on
text(-0.23,1.13, 'C', 'Units','normalized','FontWeight','bold');

switch panel_C_avg
    case 'mean'
        x = [pop_stats_dir(signif_both_dir_IX,1).field_size_mean];
        y = [pop_stats_dir(signif_both_dir_IX,2).field_size_mean];
        title_str = 'Mean field size (m)';
    case 'median'
        x = [pop_stats_dir(signif_both_dir_IX,1).field_size_median];
        y = [pop_stats_dir(signif_both_dir_IX,2).field_size_median];
        title_str = 'Median field size (m)';
end
c = cat(1,pop_bat_color{signif_both_dir_IX});
h=scatter(x,y,5,c,'filled');

[r,rpval] = corr(x',y','type','Pearson');
[rho,rhopval] = corr(x',y','type','Spearman');
if rhopval==0
    rhopval = realmin;
end
[~,P_ttest] = ttest(x,y);
[P_ranksum,~,STATS_ranksum] = ranksum(x,y);
[P_signtest,~,STATS_signtest] = signtest(x,y);
stats_str = {[sprintf('%s = %.2g','\rho',rho), ',  ' sprintf('P_{%s} = %.2g','\rho',rhopval)];...
              sprintf('P_{wilc} = %.3g',P_ranksum);};
text(0.06,1, stats_str, 'units',...
    'normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);
text(0.95,0.2, "n = "+sum(~any(isnan([x;y]))), 'units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);

xlabel('Direction 1','Color',prm.graphics.colors.flight_directions{1},'units','normalized','Position',[0.5 -0.14]);
ylabel('Direction 2','Color',prm.graphics.colors.flight_directions{2},'units','normalized','Position',[-0.14 0.5]);
title(title_str,'units','normalized','Position',[0.5 1.05])
axis equal
xlim([0 max([x y])+1])
ylim([0 max([x y])+1])
h=refline(1,0);
h.Color = 0.5*[1 1 1];

%% panel D - Field size ratio (per direction)
axes(panel_D);
cla
hold on
text(-0.23,1.13, 'D', 'Units','normalized','FontWeight','bold');

x = [pop_stats_dir(signif_both_dir_IX,1).field_ratio_LS];
y = [pop_stats_dir(signif_both_dir_IX,2).field_ratio_LS];
c = cat(1,pop_bat_color{signif_both_dir_IX});
h=scatter(x,y,5,c,'filled');

[r,rpval] = corr(x',y','type','Pearson','rows','complete');
[rho,rhopval] = corr(x',y','type','Spearman','rows','complete');
if rhopval==0
    rhopval = realmin;
end
[~,P_ttest] = ttest(x,y);
[P_ranksum,~,STATS_ranksum] = ranksum(x,y);
[P_signtest,~,STATS_signtest] = signtest(x,y);
stats_str = {[sprintf('%s = %.2g','\rho',rho), ',  ' sprintf('P_{%s} = %.2g','\rho',rhopval)];...
              sprintf('P_{wilc} = %.2g',P_ranksum);};
text(0.06,1, stats_str, 'units',...
    'normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',7);
text(1.05,0.28, "n = "+sum(~any(isnan([x;y]))), 'units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);

xlabel('Direction 1','Color',prm.graphics.colors.flight_directions{1},'units','normalized','Position',[0.5 -0.14]);
ylabel('Direction 2','Color',prm.graphics.colors.flight_directions{2},'units','normalized','Position',[-0.14 0.5]);
title('Field size ratio, largest/smallest','units','normalized','Position',[0.5 1.05])
axis equal
xlim([1 max([x y])+1])
ylim([1 max([x y])+1])
hax=gca;
hax.XScale = panel_D_scale;
hax.YScale = panel_D_scale;
hax.XTick = [1 10];
hax.YTick = [1 10];
hax.TickLength = [0.02 0.02];
h=refline(1,0);
h.Color = 0.5*[1 1 1];


%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
fig_name_out = [fig_name_out '_panel_C_' panel_C_avg];
fig_name_out = [fig_name_out '_panel_D_' panel_D_scale];
if color_by_bat
    fig_name_out = [fig_name_out '_color_by_bat'];
end
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');








