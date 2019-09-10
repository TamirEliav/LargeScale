%% Large Scale - fig 2 - Behavioral and neural recordings from bats fliying over large spatial scales.

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_3';
fig_caption_str = 'Place cells are distributed uniformly across the tunnel, with slight over-representation at particular landmarks';
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

pause(0.2); % workaround to solve matlab automatically changing the axes positions...

%% create panels
panel_A_size = [8 11];
panel_B_size = [8 11];
panel_C_size = [8 2];
panel_D_size = [2 2];
panel_E_size = [8 2];
panel_F_size = [2 2];
panel_G_size = [4 3];
panel_H_size = [4 3];
panel_A(1,1) = axes('position', [ 2  13 panel_A_size ]);
panel_A(2,1) = axes('position', [ 11 13 panel_A_size ]);
panel_A(1,2) = axes('position', [ 2  13+panel_A_size(2) panel_A_size.*[1 .2] ]);
panel_A(2,2) = axes('position', [ 11 13+panel_A_size(2) panel_A_size.*[1 .2] ]);
panel_B(1) = axes('position', [ 2   1 panel_B_size ]);
panel_B(2) = axes('position', [ 11  1 panel_B_size ]);
% panel_C(1) = axes('position', [ 2   5 panel_C_size ]);
% panel_C(2) = axes('position', [ 11  5 panel_C_size ]);
% panel_D(1,1) = axes('position', [ 2  1 panel_D_size ]);
% panel_D(1,2) = axes('position', [ 5  1 panel_D_size ]);
% panel_D(1,3) = axes('position', [ 8  1 panel_D_size ]);
% panel_D(2,1) = axes('position', [11  1 panel_D_size ]);
% panel_D(2,2) = axes('position', [14  1 panel_D_size ]);
% panel_D(2,3) = axes('position', [17  1 panel_D_size ]);
% panel_E(1) = axes('position', [ 2  1 panel_E_size ]);
% panel_E(2) = axes('position', [ 11 1 panel_E_size ]);
% panel_F(1,1) = axes('position', [ 2  6 panel_F_size ]);
% panel_F(1,2) = axes('position', [ 5  6 panel_F_size ]);
% panel_F(1,3) = axes('position', [ 8  6 panel_F_size ]);
% panel_F(2,1) = axes('position', [11  6 panel_F_size ]);
% panel_F(2,2) = axes('position', [14  6 panel_F_size ]);
% panel_F(2,3) = axes('position', [17  6 panel_F_size ]);
% panel_G = axes('position', [2 2 panel_H_size ]);
% panel_H = axes('position', [8 2 panel_H_size ]);

% error

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% load population data
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
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% plot population fields distribution - lines and edges
axes(panel_A(1,2));
cla
hold on
text(-0.12,1, 'A', 'Units','normalized','FontWeight','bold');
% xlabel('Position (m)')
% ylabel('Cell no.')

axes(panel_A(1,1));
cla
hold on
text(-0.07,1.1, 'A', 'Units','normalized','FontWeight','bold');
% xlabel('Position (m)')
ylabel('Cell no.')

vline_w = 1;
vline_h = 1*[-1 1];

color_edges_L = 'g';
color_edges_R = 'r';
color_peaks = 'm';

% color_edges_L = 'g';
% color_edges_R = [255 127 0]/255;
% color_peaks = 'm';

% color_edges_L = [102,194,165]/255;
% color_edges_R = [166,216,84]/255;
% color_peaks = 'm';

for ii_dir = 1:2
    axes(panel_A(ii_dir,1));
    cla
%     yyaxis left
    hold on
    signif = cat(1,cells.signif);
    cells_dir = cells([signif(:,ii_dir).TF]);
    stats = cat(1,cells_dir.stats);
    stats = cat(1,stats.dir);
    SI = [stats(:,ii_dir).SI_bits_spike];
    [SI,sort_IX] = sort(SI,'descend');
    cells_dir = cells_dir(sort_IX);
    fields_edges = [];
    fields_peaks = [];
    for ii_cell = 1:length(cells_dir)
        cell = cells_dir(ii_cell);
        fields = cell.fields{ii_dir};
        x = cat(1,fields.edges_prc);
        fields_edges = [fields_edges ; x];
        y = ii_cell;
        plot( x',[y y],'-','Color',prm.graphics.colors.flight_directions{ii_dir},'LineWidth',0.5);
        plot( x(:,[1 1])', y+vline_h,'-','Color',color_edges_L,'LineWidth',vline_w);
        plot( x(:,[2 2])', y+vline_h,'-','Color',color_edges_R,'LineWidth',vline_w);
        x = cat(1,fields.loc);
        fields_peaks = [fields_peaks ; x];
        plot( [x x]', y+vline_h,'-','Color',color_peaks,'LineWidth',vline_w);
    end
    ylim([0 length(cells_dir)+1])
    
    % plot peak/edges distribution
    axes(panel_A(ii_dir,2));
    cla
%     yyaxis right
    hold on
    x = 0:0.2:200;
    [f,xi] = ksdensity(fields_peaks,x,'Kernel','normal','Bandwidth',1.5);
    plot(xi,f, '-', 'LineWidth',0.01, 'Color',color_peaks);
    [f,xi] = ksdensity(fields_edges(:,1),x,'Kernel','normal','Bandwidth',1.5);
    plot(xi,f, '-', 'LineWidth',0.01, 'Color',color_edges_L);
    [f,xi] = ksdensity(fields_edges(:,2),x,'Kernel','normal','Bandwidth',1.5);
    plot(xi,f, '-', 'LineWidth',0.01, 'Color',color_edges_R);
    ha = gca;
    ha.XLim = [0 200];
    ha.XTickLabel = {};
%     hl=columnlegend(3,{' peak      ','         left ',' right '});
%     hl.Position = ha.Position.*[1 1 0.3 0.1] + ha.Position([3 4 1 2]).*[0.3 0.8 0 0];
end

%% plot FR maps ordered by field location (heatmap)
axes(panel_B(1));
hold on
text(-0.07,1.08, 'B', 'Units','normalized','FontWeight','bold');
xlabel('Position (m)')
ylabel('Field no.')
for ii_dir = 1:2
    axes(panel_B(ii_dir));
    hold on
    xlabel('Position (m)')
    
    %% arrange data
    signif = cat(1,cells.signif);
    cells_dir = cells([signif(:,ii_dir).TF]);
    pop_data = struct('field',{}, 'PSTH',{}, 'fields',{});
    ii_field = 1;
    for ii_cell = 1:length(cells_dir)
        cell = cells_dir(ii_cell);
        fields = cell.fields{ii_dir};
        for jj_field = 1:length(fields)
            pop_data(ii_field).field = fields(jj_field);
            pop_data(ii_field).PSTH = cell.FR_map(ii_dir).all.PSTH;
            pop_data(ii_field).fields = fields;
            ii_field = ii_field + 1;
        end
    end
    locs = arrayfun(@(x)(x.loc), [pop_data.field]);
    [locs, sort_IX] = sort(locs,'ascend');
    pop_data = pop_data(sort_IX);
    
    PSTH_pos_bins = cell.FR_map(ii_dir).all.bin_centers;
    M = nan(length(pop_data),length(PSTH_pos_bins));
    for ii_field = 1:length(pop_data)
        PSTH = pop_data(ii_field).PSTH;
        PSTH = PSTH ./ pop_data(ii_field).field.peak;
%         PSTH = PSTH ./ (pop_data(ii_field).field.peak*prm.fields.width_href);
        M(ii_field, :) = PSTH;
    end
    
%     exp_list = arrayfun(@(x)(x.exp_ID), [cells_dir.details],'UniformOutput',0);
%     exps = cellfun(@(exp_ID)(exp_load_data(exp_ID,'details','LM','flight')), exp_list ,'UniformOutput',0);
%     exps = [exps{:}];
%     vel_all = arrayfun(@(x)(x.speed_traj(ii_dir).vel_median), [exps.flight],'UniformOutput',0);
%     vel_all =cat(2,vel_all{:});
%     speed_all = abs(vel_all);
%     speed_mean = nanmean(abs(vel_all),2);

    %% plot
    m = 1.5;
    M(M>m) = m;
    imagesc(PSTH_pos_bins,1:size(M,1),M)
%     colorbar
    set(gca,'CLim',[0 m]);
    colormap parula
%     colormap jet
    ylim([0 size(M,1)+1])
    
% %     yyaxis right
% %     x = exps(1).flight.speed_traj.bins_centers;
% % %     plot(x, speed_all )
% %     plot(x, speed_mean, 'k','LineWidth',2)

%     exp = exp_load_data(cell.details.exp_ID,'LM');
%     plot_LM(exp.LM)
    
    
end

%%
if 0
%%
axes(panel_C(1));
hold on
text(-0.07,1.1, 'C', 'Units','normalized','FontWeight','bold');
xlabel('')
ylabel('')


%%
axes(panel_D(1,1));
hold on
text(-0.25,1.15, 'D', 'Units','normalized','FontWeight','bold');
xlabel('')
ylabel('')

%%
axes(panel_E(1));
hold on
text(-0.07,1.1, 'E', 'Units','normalized','FontWeight','bold');
xlabel('')
ylabel('')


%%
axes(panel_F(1,1));
hold on
text(-0.25,1.15, 'F', 'Units','normalized','FontWeight','bold');
xlabel('')
ylabel('')

%%
axes(panel_G);
hold on
text(-0.13,1.1, 'G', 'Units','normalized','FontWeight','bold');
xlabel('')
ylabel('')

%%
axes(panel_H);
hold on
text(-0.15,1.1, 'H', 'Units','normalized','FontWeight','bold');
xlabel('')
ylabel('')

end


%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');






%% local function to plot lines for landmarks
function plot_LM(LM)

ylimits = get(gca,'ylim');
for ii_LM=1:length(LM)
    x = LM(ii_LM).pos_proj;
    plot(repelem(x,2) , ylimits, '-', 'color', 0.9.*[1 1 1], 'LineWidth',0.5)
    text(x, ylimits(2)+0.02*diff(ylimits), LM(ii_LM).name, 'Rotation', 45, 'FontSize',8)
end

end
