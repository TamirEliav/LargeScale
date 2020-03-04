%% Large Scale - Fig. S3 - 1-D behavior

%%
clear 
clc

%% data options
panel_A_opt = [4 6 7];

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S3';
fig_caption_str = '1-Dim behavior';
log_name_str = [fig_name_str '_log_file' '.txt'];
log_name_str = strrep(log_name_str , ':', '-');
log_name_str = strrep(log_name_str , ' ', '_');
log_name_out = fullfile(res_dir, log_name_str);

%%
% fig_name_str = sprintf('%s_opt_%d_%d_%d', fig_name_str , panel_A_opt)

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
panel_A(1) = axes('position', [ 2 9    5 15]);
panel_A(2) = axes('position', [ 8 9    5 15]);
panel_A(3) = axes('position', [14 9    5 15]);
panel_B    = axes('position', [ 4 4.6  6  2]);
panel_C    = axes('position', [12 4.6  2  2]);


%% load population data
% load params
prm = PARAMS_GetAll();
% get list of significant cells (at least in one direction)
% ---------------------------------------------------------
cells_t = DS_get_cells_summary();
% choose bats
bats = [79,148,34,9861,2289];
cells_t(~ismember(cells_t.bat, bats ),:) = [];
% choose good clusters from CA1
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.details];
cells(~contains({cells.brain_area}, 'CA1')) = [];
cells(~ismember([cells.ClusterQuality], [2])) = [];
% choose pyramidal cells only
cells = cellfun(@(c)(cell_load_data(c,'details','meanFR')), {cells.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];
meanFR = [cells.meanFR];
cells([meanFR.all]>prm.inclusion.interneuron_FR_thr)=[];
% choose only signif cells (in any direction)
cells_details = [cells.details];
cells = cellfun(@(c)(cell_load_data(c,'details','signif')), {cells_details.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];
signif = cat(1,cells.signif);
signif = arrayfun(@(x)(x.TF), signif);
signif = any(signif,2);
cells(~signif)=[];

% load the final list of cells
cells_details = [cells.details];
clear cells stats signif cells_t
cells = cellfun(@(c)(cell_load_data(c,'details','stats','meanFR','stats','inclusion','signif','fields','FR_map')), {cells_details.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];

% load exp data
cells_details = [cells.details];
exps_IDs = {cells_details.exp_ID};
exp_list = unique(exps_IDs);
exps = cellfun(@(x)(exp_load_data(x,'details','flight')), exp_list);
exps_details = [exps.details];
exps_flight = [exps.flight];
pos_y_dev_all = cat(1,exps_flight.pos_y_std);
speed_traj_all = cat(1,exps_flight.speed_traj);
ystd_median_all = arrayfun(@(x)(x.ystd_median), pos_y_dev_all);

%% panel A - Behavioral trajectories (time vs pos) - examples
exp_ID_options = {
    'b0034_d180310';
    'b0034_d180312';
    'b0034_d180313';
    'b0079_d160915';
    'b0079_d160921';
    'b0079_d160928';
    'b0079_d161004';
    'b0079_d161005';
    'b0148_d170613';
    'b0148_d170802';
    'b0148_d170807';
    'b9861_d180524';
    'b9861_d180525';
    'b9861_d180527';
    'b9861_d180529';
    'b9861_d180603';
    'b9861_d180615';
};
for ii_exp = 1:length(panel_A_opt)
    % get exp data
    exp_ID = exp_ID_options{panel_A_opt(ii_exp)};
    exp = exp_load_data(exp_ID,'flight','pos');

    x = exp.pos.proc_1D.pos;
    y = exp.pos.proc_1D.ts;
    session_ts = exp_get_sessions_ti(exp_ID, 'Behave');
    t0 = session_ts(1);
    
    % down-sample (too many dots!!)
    ds = 10;
    x = x(1:ds:end);
    y = y(1:ds:end);
    
    % plot!
    axes(panel_A(ii_exp));
    cla('reset')
    hold on
%     plot(x,y,'.k','MarkerSize',4)
    plot(x,y,'-k','LineWidth',1)
    
    % labels
    xlabel('Position (m)','Units','normalized','Position',[0.5 -0.03]);
    ylabel('Time (min)','Units','normalized','Position',[-0.1 0.5]);
    ha = gca;
    ha.XLim = [0 200];
    ha.YLim = session_ts;
%     ha.YLim = [-1.5 1.5];
    ha.XTick = [0:50:200];
%     ha.YTick = [-2:1:2];
    ha.TickDir='out';
    ha.TickLength = [0.007 0.01];
    ha.XRuler.TickLabelGapOffset = 0.5;
    ha.YRuler.TickLabelGapOffset = 1.5;
    
    rescale_plot_data('y',[1e-6/60,t0]);
    
end
axes(panel_A(1));
text(-0.15,1.03, 'A', 'Units','normalized','FontWeight','bold');

%% panel B - behavioral trajectory is 1D (small y deviations) - example
axes(panel_B);
cla
hold on
text(-0.14,1.15, 'B', 'Units','normalized','FontWeight','bold');
exp_ID = 'b0034_d180413';
exp = exp_load_data(exp_ID,'flight');
for ii_dir = 1:2
    c = prm.graphics.colors.flight_directions{ii_dir};
    ydev = exp.flight.pos_y_std(ii_dir);
    x = ydev.xy(:,1);
    y = ydev.xy(:,2);
    ymean = interp1(ydev.bin_centers, ydev.ymean, x);
    y = y-ymean;
    
%     % the data is not in FE structure, so using a line plot creates jumps.
%     % WORKAROUND: add nans in the big jumps
%     x(find(abs(diff(x)) > 10)) = nan;
%     plot(x, y, '-', 'Color',c, 'LineWidth',0.0001);

    plot(x, y, '.', 'Color',c, 'MarkerSize',.0001);
end
xlabel('Position (m)','Units','normalized','Position',[0.5 -0.2]);
ylabel('Y (m)','Units','normalized','Position',[-0.055 0.5]);
ha = gca;
ha.XLim = [0 200];
ha.YLim = [-1.5 1.5];
ha.XTick = [0:50:200];
ha.YTick = [-2:1:2];
ha.TickDir='out';
ha.TickLength = [0.01 0.01];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

% add direction arrows
arrow_x = [0 0.05] + 0.38;
arrow_y = repelem(0.24,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.01, 'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y     , 'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);

%% panel C - behavioral trajectory is 1D (small y deviations) - Population
axes(panel_C);
cla
text(-0.4,1.15, 'C', 'Units','normalized','FontWeight','bold');
hold on
% arrange data
data = {};
for ii_bat = 1:length(bats)
    bat = bats(ii_bat);
    c = prm.graphics.colors.bats(bat);
    bat_exp_IX = ismember([exps_details.batNum],bat);
    bat_exp_ydev = pos_y_dev_all(bat_exp_IX,:);
    valid_pos_IX = get_data_in_ti(bat_exp_ydev(1).bin_centers, prm.fields.valid_speed_pos);
    ystd_bat_all = cat(1,bat_exp_ydev.ystd); % pooling directions together!
    ystd_bat_all = ystd_bat_all(:,valid_pos_IX); % take only valid positions
    ystd_bat_all = ystd_bat_all(:)';% pooling over all positions!
    ystd_bat_all = ystd_bat_all.*100; % convert to cm
    data{ii_bat} = ystd_bat_all;
end
% t2 = table2cell(t);
x = [data{:}];
grps = repelem(1:length(bats), cellfun(@length,data));

% plot!
cmap = prm.graphics.colors.bats;
c = arrayfun(@(x)(cmap(x)), bats,'UniformOutput',0);
yvar_pop_plot = 'ksdensity_pooled_bats';
% yvar_pop_plot = 'violin';
% yvar_pop_plot = 'boxplot';
% yvar_pop_plot = 'median_std';
% yvar_pop_plot = 'none';
switch yvar_pop_plot 
    case 'ksdensity_pooled_bats'
        for ii_dir = 1:2
            ystd_by_dir = cat(1,pos_y_dev_all(:,ii_dir).ystd); % pool over sessions
            ystd_by_dir = ystd_by_dir(:,valid_pos_IX); % take only valid position (high-speed)
            ystd_by_dir = ystd_by_dir(:); % pool over positions
            ystd_by_dir = ystd_by_dir .* 100; % convert to cm
            xi = linspace(0,50,100);
            ystd_density = ksdensity(ystd_by_dir,xi);
            plot(xi,ystd_density,'Color',prm.graphics.colors.flight_directions{ii_dir},'LineWidth',1.5);
        end
        xlabel('Y s.d. (cm)');
        ylabel('Probability','Units','normalized','Position',[-0.1 0.5]);
        ha=gca;
        ha.XLim = [0 40];
        ha.XTick = 0:10:40;
        ha.YLim = [0 0.1];
        ha.YTick = [0 0.1];
        ha.TickLength = [0.03 0.03];
        ha.XRuler.TickLabelGapMultiplier = -0.3;
        ha.YRuler.TickLabelGapMultiplier = 0.001;
    case 'violin'
        hv=violinplot(x,grps);
        [hv.BoxWidth] = disperse(repelem(0.02,length(hv)));
        [hv.EdgeColor ] = disperse(repelem({'none'},length(hv)));
        hs = [hv.ScatterPlot];
        [hs.SizeData] = disperse(repelem(5,length(hs)));
        hm = [hv.MedianPlot];
        [hm.SizeData] = disperse(repelem(10,length(hm)));
        ha=gca;
        ha.YLim = [0 90];
        ha.YTick = [0:30:90];
    case 'boxplot'
        boxplot(x,grps)
        box off
        ha=gca;
        ha.YLim = [0 90];
        ha.YTick = [0:30:90];
    case 'median_std'
        medians = cellfun(@median,data);
        stds    = cellfun(@std,data);
        errorbar(medians,stds,'-o','MarkerSize',3);
        ha=gca;
        ha.YLim = [0 40];
        ha.YTick = [0:20:40];
    case 'none'
        ha=gca;
end
switch yvar_pop_plot
    case {'ksdensity_pooled_bats'}
    otherwise
        ha.XTick = 1:length(bats);
        ha.XTickLabel = bats;
        ha.XLim = [0.5 length(bats)+.5];
        ha.TickDir='out';
        ha.TickLength = [0.01 0.01];
        ha.XRuler.TickLabelGapMultiplier = -0.3;
        ha.YRuler.TickLabelGapMultiplier = 0.001;
        xlabel('Bats');
        ylabel('Y s.d. (cm)');
end

% add direction arrows
arrow_x = [0 0.03] + 0.61;
arrow_y = repelem(0.24,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.01, 'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y     , 'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);



%% print/save the figure
fig_name_out = fullfile(res_dir, sprintf('%s_ds=%d',fig_name_str,ds));
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');


%%
