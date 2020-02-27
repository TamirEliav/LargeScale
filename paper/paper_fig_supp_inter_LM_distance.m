%% Large Scale - supp fig - Inter LM distance

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S9';
fig_caption_str = 'field properties vs. distance between landmrks';
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
panel_A_size = [3 3];
panel_B_size = [7 3];
panel_A(1) = axes('position', [ 4    21    panel_A_size]);
panel_A(2) = axes('position', [ 4    15.8  panel_A_size]);
panel_B(1) = axes('position', [10.6  21    panel_B_size]);
panel_B(2) = axes('position', [10.6  15.8  panel_B_size]);


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
clear cells stats cells_details cells_t
cells = cellfun(@(c)(cell_load_data(c,'details','stats','meanFR','stats','inclusion','signif','fields','FR_map')), cells_ID, 'UniformOutput',0);
cells = [cells{:}];
whos cells

% load LM data
exp_ID = 'b2289_d180615';
exp = exp_load_data(exp_ID,'LM');
LM = exp.LM;
LM( contains({LM.name},{'ball','enter'}) ) = [];

%% panel A - field size vs. inter-LM-distanance
% =========================================================================
% figure
for ii_dir = 1:2

    %% arrange data
    signif = cat(1,cells.signif);
    cells_dir = cells([signif(:,ii_dir).TF]);
    fields = cellfun(@(x)(x{ii_dir}), {cells_dir.fields},'UniformOutput',0);
    fields =[fields{:}];
    switch ii_dir
        case 1
            start_ind = 1;
            end_ind = 2;
            interp_next = 'next';
            interp_prev = 'previous';
        case 2
            start_ind = 2;
            end_ind = 1;
            interp_next = 'previous';
            interp_prev = 'next';
    end
    edges = cat(1,fields.edges_prc);
    [fields.start] = disperse(edges(:,start_ind));
    [fields.end]   = disperse(edges(:,end_ind  ));
  
    [fields.LM_nearest_by_peak]  = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.loc]  , 'nearest'));
    [fields.LM_nearest_by_start] = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.start] , 'nearest'));
    [fields.LM_nearest_by_end]   = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.end]   , 'nearest'));
    [fields.LM_next_by_peak]  = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.loc]  , interp_next));
    [fields.LM_next_by_start] = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.start] , interp_next));
    [fields.LM_next_by_end]   = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.end]   , interp_next));
    [fields.LM_prev_by_peak]  = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.loc]  , interp_prev));
    [fields.LM_prev_by_start] = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.start] , interp_prev));
    [fields.LM_prev_by_end]   = disperse(interp1( [LM.pos_proj], [LM.pos_proj], [fields.end]   , interp_prev));
    
    % remove fields near balls
    fields([fields.in_low_speed_area])=[];
    fields(isnan([fields.LM_nearest_by_peak])) = [];
    
    %% scatter plot - field size vs. inter-LM-dist
    axes(panel_B(ii_dir));
    cla
    hold on
    c = prm.graphics.colors.flight_directions{ii_dir};
    x = abs([fields.LM_prev_by_peak]-[fields.LM_next_by_peak]);
    y = [fields.width_prc];
    rng(0);
%     plot(x, y, '.', 'Color',c); % no jitter
%     plot(x+0.5*(rand(size(x))-.5), y, '.', 'Color',c); % uniform jitter
    plot(x+0.1*randn(size(x)), y, '.', 'Color',c); % gaussian jitter
    lm = fitlm(x,y);
    [r,r_pval] = corr(x',y','type','Pearson');
    [rho,rho_pval] = corr(x',y','type','Spearman');
    text(0.95,1.07,{    ['Spearman {\rho} = '   sprintf('%.2f',rho) ', P = '      sprintf('%.2f',rho_pval) ];...
                        ['Pearson r = '         sprintf('%.2f',r)   ', P = '      sprintf('%.2f',r_pval)   ]},...
                    'Units','normalized', 'HorizontalAlignment','right', 'FontSize',7);
    xlabel('Inter-Landmark distance (m)','Units','normalized','Position',[0.5 -0.13])
    ylabel('Field size (m)','Units','normalized','Position',[-0.09 0.5]);
    xlim([0 25]);
    ylim([0 35]);
    set(gca,'xtick',[0:5:25])
    set(gca,'ytick',[0:5:35])
    ha=gca;
    ha.TickDir = 'out';
    ha.TickLength = repelem(0.01,2);
    ha.XRuler.TickLabelGapMultiplier = -.3;
    ha.YRuler.TickLabelGapMultiplier = 0.1;

    %% scatter plot - field size vs. distance to nearest LM
    axes(panel_A(ii_dir));
    cla
    hold on
    x = abs([fields.loc] - [fields.LM_nearest_by_peak]);
    y = [fields.width_prc];
    plot(x,y,'.', 'Color', prm.graphics.colors.flight_directions{ii_dir})
    [r,r_pval] = corr(x',y','type','Pearson');
    [rho,rho_pval] = corr(x',y','type','Spearman');
    text(1.15,1.07,{    ['Spearman {\rho} = '   sprintf('%.2f',rho) ', P = '      sprintf('%.2f',rho_pval) ];...
                        ['Pearson r = '         sprintf('%.2f',r)   ', P = '      sprintf('%.2f',r_pval)   ]},...
                    'Units','normalized', 'HorizontalAlignment','right', 'FontSize',7);

    % set axis properties
    ha = gca;
    ha.YLim = [0 35];
    ha.YTick = 0:5:35;
    ha.TickDir='out';
    ha.TickLength = [0.025 0.025];
    ha.XRuler.TickLabelGapMultiplier = -0.2;
    ha.YRuler.TickLabelGapMultiplier = 0.2;
    xlabel({'Distance of fields';'to nearest landmark (m)'}, 'Units','normalized','Position',[0.5 -0.11]);
    ylabel('Field size (m)', 'Units','normalized','Position',[-0.19 0.5]);
    
end
axes(panel_A(1));
text(-0.35,1.2, 'A', 'Units','normalized','FontWeight','bold');
axes(panel_B(1));
text(-0.16,1.2, 'B', 'Units','normalized','FontWeight','bold');

%% add direction arrows
arrow_x = 0.29 +[0 0.04];
arrow_y = repelem(0.87,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.008,  'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y      ,  'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);


%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');







%%

