%% Large Scale - Fig. SXXX - field-size vs. field-distance-to-nearest-LM

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'Fig_SXXX_field_size_vs_dist2LM';
fig_caption_str = 'Field size vs. field distance-to-LM';
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
panels_size = [4 4];
panel_A(1) = axes('position', [ 3 20 panels_size]);
panel_A(2) = axes('position', [ 9 20 panels_size]);

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

% load LM data
exp_ID = 'b2289_d180615';
exp = exp_load_data(exp_ID,'LM');
LM = exp.LM;
LM( contains({LM.name},{'ball','enter'}) ) = [];

%% panel A
for ii_dir = 1:2
    axes(panel_A(ii_dir));
    cla
    hold on
    
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
    
    %% plot
    x = abs([fields.loc] - [fields.LM_nearest_by_peak]);
    y = [fields.width_prc];
    plot(x,y,'.', 'Color', prm.graphics.colors.flight_directions{ii_dir})

    % set axis properties
    ha = gca;
    ha.TickDir='out';
    ha.TickLength = [0.025 0.025];
    ha.XRuler.TickLabelGapMultiplier = -0.2;
    ha.YRuler.TickLabelGapMultiplier = 0.001;
    xlabel({'Distance of fields';'to nearest landmark (m)'}, 'Units','normalized','Position',[0.5 -0.11]);
    ylabel('Field size (m)', 'Units','normalized','Position',[-0.15 0.5]);
    
end


%% add direction arrows
arrow_vec = 0.025*[-1 1];
arrow_y = 0.915*[1 1];
clear h
h(1)=annotation('arrow',      arrow_vec +0.23, arrow_y, 'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow', flip(arrow_vec)+0.51, arrow_y, 'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);



%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

