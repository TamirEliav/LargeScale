%% Large Scale - supp fig - speed modualtion on AVERAGED field size

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_supp_speed_effects';
fig_caption_str = 'speed effects on AVERAGED field size';
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
panel_AB_size = [3 4];
panel_AB(1) = axes('position', [ 6 20 panel_AB_size]);
panel_AB(2) = axes('position', [10 20 panel_AB_size]);
panel_AB(3) = axes('position', [ 6 15 panel_AB_size]);
panel_AB(4) = axes('position', [10 15 panel_AB_size]);
panel_AB(5) = axes('position', [16 17.5 panel_AB_size]);

panel_stats_text = axes('position', [1 9 4 12]);

panel_AB2_size = [3 3];
panel_AB2(1) = axes('position', [ 6 20 panel_AB2_size]+[0 -12 0 0]);
panel_AB2(2) = axes('position', [10 20 panel_AB2_size]+[0 -12 0 0]);
panel_AB2(3) = axes('position', [ 6 15 panel_AB2_size]+[0 -12 0 0]);
panel_AB2(4) = axes('position', [10 15 panel_AB2_size]+[0 -12 0 0]);
panel_AB2(5) = axes('position', [16 17.5 panel_AB2_size]+[0 -12 0 0]);


%% load population data
% =========================================================================
prm = PARAMS_GetAll();
cells_t = DS_get_cells_summary();
bats = [148,34,9861,2289,79];
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
signif_both_dir = all(pop_signif_TF,2);
pop_details = cat(1,cells.details);
pop_bat_number = [pop_details.bat];

%% panels A&B
% =========================================================================
averaging_method = 'mean';
% averaging_method = 'median';
% averaging_method = 'AC_width';
stats_str = {};
for ii_bat = 1:length(bats)

    %% arrange data
    bat = bats(ii_bat);
    c = prm.graphics.colors.bats(bat);
    bat_cells = ismember(pop_bat_number, bat);
    IX = signif_both_dir & bat_cells';
    switch averaging_method
        case 'mean'
            x1 = [pop_stats_dir(IX,1).field_size_mean];
            x2 = [pop_stats_dir(IX,2).field_size_mean];
        case 'median'
            x1 = [pop_stats_dir(IX,1).field_size_median];
            x2 = [pop_stats_dir(IX,2).field_size_median];
        case 'AC_width'
            x1 = [pop_stats_dir(IX,1).AC_width];
            x2 = [pop_stats_dir(IX,2).AC_width];
    end
    
    % calc stats
    [~,ttest_pval,~,ttest_stats] = ttest(x1,x2);
    signtest_pval = signtest(x1,x2);
    signrank_pval = signrank(x1,x2);
    stats_str{end+1} = sprintf('bat %d',bat);
    stats_str{end+1} = sprintf('ttest, P=%.3f',ttest_pval);
    stats_str{end+1} = sprintf('signtest, P=%.3f',signtest_pval);
    stats_str{end+1} = sprintf('signrank, P=%.3f',signrank_pval);
    stats_str{end+1} = '';
    
    %% plot (graphical option 1)
    axes(panel_AB(ii_bat));
    cla
    hold on
    plot([x1;x2],'.-', 'Color',c);
    ha = gca;
    ha.XTick = [1 2];
    ha.XTickLabel = "Direction "+[1 2];
    ha.XLim = [0.5 2.5];
    switch averaging_method
        case 'mean'
            ylabel('Averaged field size (m)');
        case 'median'
            ylabel('Median field size (m)');
        case 'AC_width'
            ylabel('AC width (m)');
    end
    xlabel('');
    title("bat "+bat)
    
    %% plot (graphical option 2)
    axes(panel_AB2(ii_bat));
    cla
    hold on
    axis equal
    plot(x1,x2,'.', 'Color',c);
    ha = gca;
    ha.XLim = [0 max([x1,x2])+1];
    ha.YLim = [0 max([x1,x2])+1];
    xlabel('Direction 1');
    ylabel('Direction 2');
    title("bat "+bat)
    h=refline(1,0);
    h.Color = 'k';
end

axes(panel_AB(1));
text(-0.2,1.1, 'A', 'Units','normalized','FontWeight','bold');
axes(panel_AB(5));
text(-0.2,1.1, 'B', 'Units','normalized','FontWeight','bold');

%%
axes(panel_stats_text)
cla
set(gca,'Visible','off');
text(0,1, stats_str, 'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top');

%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
fig_name_out = [fig_name_out '_' averaging_method]
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');







%%

