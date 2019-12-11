%% Large Scale - supp fig - speed modualtion on AVERAGED field size

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S5';
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
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');
pause(0.2); % workaround to solve matlab automatically changing the axes positions...

% create panels
panel_B_size = [2 3.5];
panel_A       = axes('position', [ 2 18 12 4]);
panel_A_ex(1) = axes('position', [ 2 13 4 2.5]);
panel_A_ex(2) = axes('position', [ 7 13 4 2.5]);
panel_A_ex(3) = axes('position', [12 13 4 2.5]);
panel_B(1) = axes('position', [ 2 7 panel_B_size]);
panel_B(2) = axes('position', [ 5 7 panel_B_size]);
panel_B(3) = axes('position', [ 8 7 panel_B_size]);
panel_B(4) = axes('position', [11 7 panel_B_size]);
panel_B(5) = axes('position', [14 7 panel_B_size]);

% panel_stats_text = axes('position', [16 13 4 12]);

%% load population data
% =========================================================================
prm = PARAMS_GetAll();
cells_t = DS_get_cells_summary();
bats = [2289 9861 148 34 79];
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

%% panel A - population speed comparison per bat

% arrange data
signif = cat(1,cells.signif);
signif = arrayfun(@(x)(x.TF), signif);
signif = any(signif,2);
signif_cells_details = [cells(signif).details];
signif_cells_expIDs = {signif_cells_details.exp_ID};
[signif_cells_expIDs_unique,IA,IC] = unique(signif_cells_expIDs, 'stable');

% load exp data
exps = cellfun(@(c)(exp_load_data(c,'details','flight')), signif_cells_expIDs_unique, 'UniformOutput',0);
exps = [exps{:}];
exps_details = [exps.details];

%
bins_centers = exps(1).flight.speed_traj(1).bins_centers;
valid_IX = get_data_in_ti(bins_centers, prm.fields.valid_speed_pos);
speed = nan(length(exps), 2, length(valid_IX));
for ii_exp = 1:length(exps)
    for ii_dir = 1:2
        vel = exps(ii_exp).flight.speed_traj(ii_dir).vel_median;
        speed(ii_exp,ii_dir,:) = abs(vel(valid_IX));
    end
end
exps_bat_IX = interp1(bats,1:length(bats), [exps_details.batNum]);
n_sessions_per_bat = accumarray(exps_bat_IX',ones(size(exps_bat_IX))');

speed_mean = [];
speed_sem = [];
speed_pval = [];
for ii_bat = 1:length(bats)
    bat_IX = exps_bat_IX == ii_bat;
    for ii_dir = 1:2
        M = squeeze(speed(bat_IX, ii_dir, :));
        M = mean(M,2);
        bat_speed{ii_bat,ii_dir} = M;
        speed_mean(ii_bat,ii_dir) = nanmean(M(:));
        speed_sem(ii_bat,ii_dir) = nansem(M(:));
    end
    % calc stats
    x = bat_speed{ii_bat,1}(:);
    y = bat_speed{ii_bat,2}(:);
    speed_pval(ii_bat) = ranksum(x,y,'tail','both');
end
speed_pval

%%
axes(panel_A); 
cla
hold on
text(-0.1,1.2, 'A', 'Units','normalized','FontWeight','bold');

x = [1:length(bats)]' +0.15.*[-1 1];
y = speed_mean;
err = speed_sem;
hbar = bar(x,y);
[hbar.FaceColor] = disperse( prm.graphics.colors.flight_directions );
[hbar.BarWidth] = disperse( repelem(3,2) );
herr = errorbar(x,y,err);
[herr.LineStyle] = disperse(repelem('none',2));
[herr.Color] = disperse( prm.graphics.colors.flight_directions );
herr(1).XData = herr(1).XData-0.04; % workaround...
herr(2).XData = herr(2).XData+0.04; % workaround...

% draw significance comparison lines
for ii_bat = 1:length(bats)
    plot(x(ii_bat,[1 1]), [y(ii_bat,1)+0.5 9.5],'k-','LineWidth',1.1)
    plot(x(ii_bat,[2 2]), [y(ii_bat,2)+0.5 9.5],'k-','LineWidth',1.1)
    plot(x(ii_bat,[1 2]), [9.5 9.5],'k-','LineWidth',1.1)
end
significant_strs = {'n.s';'n.s';'n.s';'n.s';'*'};
ht = text( mean(x,2), repelem(11.5,length(bats)), "n="+n_sessions_per_bat,...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontSize',8);
ht = text( mean(x,2), repelem(10.5,length(bats)), "sessions",...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontSize',8);
ht = text( mean(x,2), repelem(9.6,length(bats)), significant_strs,...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom');
[ht.FontSize] = disperse([8 8 8 8 15 ]);
ht(end).Position(2) = 9;

xlim([0.5 5.5])
ylim([0 10])

ha = gca;
ha.XTick = 1:length(bats);
ha.XTickLabel = {};
ha.TickDir = 'out';
ticklabels = arrayfun(@(bat)({"Bat";num2str(bat)}),bats,'UniformOutput',0)
text(1:length(bats), repelem(-0.3,length(bats)), ticklabels, ...
    'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8);
ylabel('Speed (m/s)','Units','normalized','Position',[-0.05 0.45]);

%% panel A - speed trajectory exmaples
example_opt = [20 15 4];
example_list = {
    'b0079_d160915'; % 1
    'b0079_d160920';
    'b0079_d160921';
    'b0079_d160925';
    'b0079_d160927'; % 5
    'b0079_d160928';
    'b0079_d160929';
    'b0079_d160930';
    'b0079_d161004';
    'b0079_d161005'; % 10
    
    'b0034_d180313'; % 11
    'b0034_d180315';
    'b0079_d161003';
    'b0148_d170613';
    'b0148_d170615'; % 15
    'b0148_d170626';
    'b0148_d170802';
    'b0148_d170807';
    'b2289_d180520';
    'b2289_d180615'; % 20
};
% ylimits = [9 10];
for ii_ex = 1:3
    axes(panel_A_ex(ii_ex));
    cla
    hold on
    exp = exp_load_data(example_list{example_opt(ii_ex)},'flight');
    FE=exp.flight.FE;
    FE([exp.flight.FE.distance]<prm.flight.full_min_distance) = [];
    directions = [1 -1];
    for ii_dir = 1:2
        dir_IX = [FE.direction]==directions(ii_dir);
        FE_dir = [FE(dir_IX)];
        x = [FE_dir.pos];
        y = [FE_dir.vel];
        y = abs(y);
        plot(x,y,'.','Color', prm.graphics.colors.flight_directions{ii_dir},'MarkerSize',1);
    end
    set(gca,'xtick',0:50:200,'ytick',[-10 0 8],'xlim',[0 200])
    set(gca,'tickdir','out','TickLength',repelem(0.01,2));
    xlim([0 200]);
%     ylim([0 ylimits(ii_ex)]);
    ylim([0 10]);
    xlabel('Position (m)','Units','normalized','Position',[0.5 -0.2]);
    ylabel('Speed (m/s)','Units','normalized','Position',[-0.05 0.4]);
end

%% add lines going from bat to example
annotation('arrow',[0.150 0.18], [0.64 0.6],'LineWidth',0.8,'HeadWidth',6,'HeadLength',6);
annotation('arrow',[0.375 0.41], [0.64 0.6],'LineWidth',0.8,'HeadWidth',6,'HeadLength',6);
annotation('arrow',[0.595 0.64], [0.64 0.6],'LineWidth',0.8,'HeadWidth',6,'HeadLength',6);

%% panel B - compare field size (dir1 vs. dir 2)
% =========================================================================
% averaging_method = 'mean';
averaging_method = 'median';
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
    
    %% plot
    axes(panel_B(ii_bat));
    cla
    hold on
    switch 2
        case 1
            plot([x1;x2],'.-', 'Color',c);
        case 2
            plot([x1;x2],'-', 'Color',.5*[1 1 1]);
            plot(1,x1,'.','Color',prm.graphics.colors.flight_directions{1});
            plot(2,x2,'.','Color',prm.graphics.colors.flight_directions{2});
    end
    ha = gca;
    ha.XTick = [1 2];
    ha.XTickLabel = {};
    ha.TickDir = 'out';
%     ha.XTickLabel = "Direction "+[1 2];
%     ha.XTickLabel = {'{\rightarrow}';'{\leftarrow}'};
%     ha.XTickLabelRotation = 0;
    ha.XLim = [0.5 2.5];
    ha.YLim(1) = 0;
    ht = text( [1 2], repelem(ha.YLim(1)-0.16*diff(ha.YLim),2), {'{\rightarrow}';'{\leftarrow}'},...
          'HorizontalAlignment','center', 'VerticalAlignment','bottom');
    [ht.FontSize] = disperse([10 10]);
    [ht.Color] = disperse( prm.graphics.colors.flight_directions );
    [ht.FontWeight] = disperse(repelem('bold',2));
    switch averaging_method
        case 'mean'
            ylabel('Averaged field size (m)');
        case 'median'
            ylabel('Median field size (m)');
        case 'AC_width'
            ylabel('AC width (m)');
    end
    xlabel('');
    title("Bat "+bat)

end

axes(panel_B(1));
text(-0.6,1.1, 'B', 'Units','normalized','FontWeight','bold');

%% stats panel
% axes(panel_stats_text)
% cla
% set(gca,'Visible','off');
% text(0,1, stats_str, 'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top');

%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
fig_name_out = [fig_name_out '_' averaging_method];
fig_name_out = [fig_name_out '_opt_' char(join(""+example_opt,"_"))]
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');







%%

