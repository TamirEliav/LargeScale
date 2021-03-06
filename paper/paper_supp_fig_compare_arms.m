%% Large Scale - supp fig 6 - comparing the two arms

%%
clear 
clc

%% plotting options
% yscale_opt = 'linear';
yscale_opt = 'log';
long_arm_lw      = 2;
long_arm_lc      = 0.3*[1 1 1];
entire_tunnel_lw = 1;
entire_tunnel_lc = 0.0*[1 1 1];

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S7';
fig_caption_str = 'compare firing patterns between the two arms (long vs. short)';
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
panel_A = axes('position', [ 4    18 panels_size]);
panel_B = axes('position', [10.8  18 panels_size]);
panel_A_legend = axes('position', [ 6.5 21 0.5 0.3]);
panel_B_legend = axes('position', [13.3   21 0.5 0.3]);
% panel_A(1) = axes('position', [ 4   18 panels_size]);
% panel_A(2) = axes('position', [ 4   13 panels_size]);
% panel_A(3) = axes('position', [ 4    8 panels_size]);
% panel_B(1) = axes('position', [10.5 18 panels_size]);
% panel_B(2) = axes('position', [10.5 13 panels_size]);
% panel_B(3) = axes('position', [10.5  8 panels_size]);
% panel_A_legend(1) = axes('position', [ 6.5 20.5 0.5 0.3]);
% panel_A_legend(2) = axes('position', [ 6.5 15.5 0.5 0.3]);
% panel_A_legend(3) = axes('position', [ 6.5 10.5 0.5 0.3]);
% panel_B_legend(1) = axes('position', [12.5 20.5 0.5 0.3]);
% panel_B_legend(2) = axes('position', [12.5 15.5 0.5 0.3]);
% panel_B_legend(3) = axes('position', [12.5 10.5 0.5 0.3]);

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
cells = cellfun(@(c)(cell_load_data(c,'details','stats','meanFR','stats','inclusion','signif','fields','FR_map','Ipos')), cells_ID, 'UniformOutput',0);
cells = [cells{:}];
whos cells

%% load LM data
exp_ID = 'b2289_d180615';
exp = exp_load_data(exp_ID,'LM');
LM = exp.LM;
% LM( contains({LM.name},{'ball','enter'}) ) = [];
turn_point_LM = LM( contains({LM.name},{'turn-point'}) );

%% panel A - Arrange data - compare fields size between short/long arms
fields_size = builtin('cell',length(cells),2);
for ii_dir = 1:2
    for ii_cell = 1:length(cells)
        cell = cells(ii_cell);
        if ~cell.signif(ii_dir).TF % check signif per direction
            continue;
        end
        fields = cell.fields{ii_dir};
        fields([fields.in_low_speed_area]) = []; % remove fields in low speed area
        pos = [fields.loc];
        width = [fields.width_prc];
        pos_thr = turn_point_LM.pos_proj;
        IX1 = pos <  pos_thr;
        IX2 = pos >= pos_thr;
        fields_size{ii_cell,ii_dir,1} = width(IX1);
        fields_size{ii_cell,ii_dir,2} = width(IX2);
    end
end
nFields = cellfun(@length, fields_size);
nFields(nFields==0) = nan;

%% panel A - plot (per direction)
% % % % for ii_dir = 1:2
% % % %     axes(panel_A(ii_dir));
% % % %     cla
% % % %     hold on
% % % %     
% % % %     x1 = [fields_size{:,ii_dir,1}];
% % % %     x2 = [fields_size{:,ii_dir,2}];
% % % %     
% % % %     c = prm.graphics.colors.flight_directions{ii_dir};
% % % %     nBins = 15;
% % % %     h1 = histogram( x1 );
% % % %     h1.NumBins = nBins;
% % % %     h1.Normalization = 'pdf';
% % % %     h1.DisplayStyle = 'stairs';
% % % %     h1.EdgeColor = c;
% % % %     h1.LineWidth = 2; % long arm
% % % %     h2 = histogram( x2 );
% % % %     h2.NumBins = nBins;
% % % %     h2.Normalization = 'pdf';
% % % %     h2.DisplayStyle = 'stairs';
% % % %     h2.EdgeColor = c;
% % % %     h2.LineWidth = 1; % short arm
% % % %     [H,P_KS,KSSTAT] = kstest2(x1,x2);
% % % %     text(1,0.9,sprintf('P_{KS} = %.2f',P_KS),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
% % % % %     P_ranksum = ranksum(x1,x2);
% % % % %     text(1,0.85,sprintf('P_{Wilc}=%.2f',P_ranksum),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
% % % %     ha=gca;
% % % %     ha.YLim(1) = 6e-4;
% % % %     ha.YScale = yscale_opt;
% % % %     ha.XRuler.TickLabelGapOffset = -1;
% % % %     ha.TickLength = [0.03 0.03];
% % % % %     title("dir "+ii_dir);
% % % %     xlabel('Field size (m)', 'Units','normalized','Position',[0.5 -0.12]);
% % % %     ylabel('Probability', 'Units','normalized','Position',[-0.20 0.5]);
% % % % end
% % % % axes(panel_A(1));
% % % % text(-0.3,1.15, 'A', 'Units','normalized','FontWeight','bold');

%% add direction arrows
% arrow_vec = 0.025*[-1 1] + 0.29;
% arrow_y = 0.83*[1 1];
% clear h
% h(1)=annotation('arrow',      arrow_vec , arrow_y-0.01, 'Color', prm.graphics.colors.flight_directions{1});
% h(2)=annotation('arrow', flip(arrow_vec), arrow_y-0.2,      'Color', prm.graphics.colors.flight_directions{2});
% [h.HeadWidth] = disperse([5 5]);
% [h.HeadLength] = disperse([5 5]);

%% legend for long/short arms
% % % for ii_dir = 1:2
% % %     axes(panel_A_legend(ii_dir));
% % %     cla
% % %     hold on
% % %     plot([1 2], [2 2], 'Color', prm.graphics.colors.flight_directions{ii_dir}, 'LineWidth', 2); % long arm
% % %     plot([1 2], [1 1], 'Color', prm.graphics.colors.flight_directions{ii_dir}, 'LineWidth', 1); % shorty arm
% % %     xlim([0.5 2.5])
% % %     text(2.5, 1, 'Short arm', 'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',7);
% % %     text(2.5, 2, 'Long arm',  'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',7);
% % %     set(gca,'visible','off');
% % % end

%% calculate stats for comparing short vs. long arm (per direction)
for ii_dir = 1:2
    x1 = [fields_size{:,ii_dir,1}];
    x2 = [fields_size{:,ii_dir,2}];
    
    [H,P_KS,KSSTAT] = kstest2(x1,x2);
    fprintf('field size long vs. short: direction %d\t P_KS=%.2f D=%.2f n=%d m=%d\n',ii_dir,P_KS, KSSTAT, length(x1), length(x2));
end

%% panel A - plot (directions pooled)
axes(panel_A);
cla
hold on
text(-0.33,1.15, 'A', 'Units','normalized','FontWeight','bold');
x1 = [fields_size{:,:,1}];
x2 = [fields_size{:,:,:}];
x1(isnan(x1))=[];
x2(isnan(x2))=[];

nBins = 15;
% long arm
h1 = histogram( x1 );
h1.NumBins = nBins;
h1.Normalization = 'pdf';
h1.DisplayStyle = 'stairs';
h1.EdgeColor = long_arm_lc;
h1.LineWidth = long_arm_lw;
% entire tunnel
h2 = histogram( x2 );
h2.NumBins = nBins;
h2.Normalization = 'pdf';
h2.DisplayStyle = 'stairs';
h2.EdgeColor = entire_tunnel_lc;
h2.LineWidth = entire_tunnel_lw;
[H,P_KS,KSSTAT] = kstest2(x1,x2,'Tail','unequal');
fprintf('field size long vs. all\t\t\t\t\t P_KS=%.2f D=%.2f n=%d m=%d\n',P_KS, KSSTAT, length(x1), length(x2));
text(1,1.05,sprintf('P_{KS} = %.2f',P_KS),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
% P_ranksum = ranksum(x1,x2);
% text(1,0.85,sprintf('P_{Wilc}=%.3f',P_ranksum),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);

ha=gca;
ha.YLim = [0.6e-3 0.2];
ha.YScale = yscale_opt;
ha.XRuler.TickLabelGapOffset = -1;
ha.TickLength = [0.03 0.03];
ha.YTick = 10.^[-3 -2 -1];
ha.YTickLabel = {'10^{ -3}'; '10^{ -2}'; '10^{ -1}';};
xlabel('Field size (m)', 'Units','normalized','Position',[0.5 -0.14]);
ylabel('Fraction of fields', 'Units','normalized','Position',[-0.23 0.5]);

%% legend for long/short arms
axes(panel_A_legend);
cla
hold on
plot([1 2], [2 2], 'Color', long_arm_lc,      'LineWidth', long_arm_lw,      'Clipping', 'off'); % long arm
plot([1 2], [1 1], 'Color', entire_tunnel_lc, 'LineWidth', entire_tunnel_lw, 'Clipping', 'off'); % entire tunnel
xlim([0.5 2.5])
text(2.5, 2, 'Long arm',  'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',7);
text(2.5, 1, 'Entire tunnel', 'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',7);
set(gca,'visible','off');


%% field size ratio - from fields only in the long arm
% ------------------------------------------------------------------------
LS_field_ratio_all = nan(1,length(cells));
LS_field_ratio_all2 = nan(1,length(cells)); % just to double check we did it correctly as in cell_calc_stats
LS_field_ratio_all_long  = nan(1,length(cells));
LS_field_ratio_all_short = nan(1,length(cells));
LS_field_ratio_dir_long = nan(2,length(cells));
LS_field_ratio_dir_short = nan(2,length(cells));
pos_thr = turn_point_LM.pos_proj;
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    
    %% pooled stats - check at least one direction is signif
    if any([cell.signif.TF])
        LS_field_ratio_all(ii_cell) = cell.stats.all.field_ratio_LS;
        fields=[];
        for ii_dir = 1:2
            if cell.signif(ii_dir).TF 
                fields_to_add = cell.fields{ii_dir};
                % workaround to solve the problem that sometimes I don't have the
                % field 'overlap_edges'... maybe change that in 'cell_calc_fields'...
                if isfield(fields_to_add,'overlap_edges')
                    fields_to_add = rmfield(fields_to_add,'overlap_edges');
                end
                fields = [fields fields_to_add];
            end
        end
        fields([fields.in_low_speed_area])=[];
        IX1 = [fields.loc] <  pos_thr;
        IX2 = [fields.loc] >= pos_thr;
        long_arm_fields = fields(IX1);
        short_arm_fields = fields(IX2);
        if length(long_arm_fields) >= 2
            LS_field_ratio_all_long(ii_cell) = max([long_arm_fields.width_prc]) / min([long_arm_fields.width_prc]);
        end
        if length(short_arm_fields) >= 2
            LS_field_ratio_all_short(ii_cell) = max([short_arm_fields.width_prc]) / min([short_arm_fields.width_prc]);
        end
        if length(fields) >= 2
            LS_field_ratio_all2(ii_cell) = max([fields.width_prc]) / min([fields.width_prc]);
        end
    end
    
    %% per dir stats - check signif per direction
    for ii_dir = 1:2
        if cell.signif(ii_dir).TF 
            fields = cell.fields{ii_dir};
            fields([fields.in_low_speed_area])=[];
            IX1 = [fields.loc] <  pos_thr;
            IX2 = [fields.loc] >= pos_thr;
            long_arm_fields = fields(IX1);
            short_arm_fields = fields(IX2);
            if length(long_arm_fields) >= 2
                LS_field_ratio_dir_long(ii_dir,ii_cell) = max([long_arm_fields.width_prc]) / min([long_arm_fields.width_prc]);
            end
            if length(short_arm_fields) >= 2
                LS_field_ratio_dir_short(ii_dir,ii_cell) = max([short_arm_fields.width_prc]) / min([short_arm_fields.width_prc]);
            end
        end
    end
end

%% panel B - plot (per direction)
% % % for ii_dir = 1:2
% % %     axes(panel_B(ii_dir));
% % %     cla
% % %     hold on
% % % 
% % %     x1=LS_field_ratio_dir_long(ii_dir,:);
% % %     x2=LS_field_ratio_dir_short(ii_dir,:);
% % %     x1(isnan(x1))=[];
% % %     x2(isnan(x2))=[];
% % % 
% % %     nBinEdges = 9;
% % %     edges = logspace(0,log10(25),nBinEdges);
% % %     c = prm.graphics.colors.flight_directions{ii_dir};
% % %     h1=histogram(x1);
% % %     h1.BinEdges = edges;
% % %     h1.Normalization = 'pdf';
% % %     h1.DisplayStyle = 'stairs';
% % %     h1.EdgeColor = c;
% % %     h1.LineWidth = 2; % long arm
% % %     h2=histogram(x2);
% % %     h2.BinEdges = edges;
% % %     h2.Normalization = 'pdf';
% % %     h2.DisplayStyle = 'stairs';
% % %     h2.EdgeColor = c;
% % %     h2.LineWidth = 1; % short arm
% % % 
% % %     [H,P_KS,KSSTAT] = kstest2(x1, x2, 'Tail','unequal');
% % %     text(1,1,sprintf('P_{KS} = %.3f',P_KS),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
% % % %     P_ranksum = ranksum(x1,x2);
% % % %     text(1,0.95,sprintf('P_{Wilc}=%.3f',P_ranksum),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
% % % 
% % %     ha=gca;
% % %     ha.YScale = yscale_opt;
% % %     % ha.YScale = 'linear';
% % %     ha.XScale = 'log';
% % %     ha.XLim = [0 27];
% % %     ha.XTick = [1 2 5 10 20];
% % %     ha.YLim = [3e-4 1.2];
% % %     ha.YTick = 10.^[-3 -2 -1];
% % %     ha.YTickLabel = {'10^{ -3}'; '10^{ -2}'; '10^{ -1}';};
% % %     ha.TickDir='out';
% % %     ha.TickLength = [0.03 0.03];
% % %     ha.XRuler.TickLabelGapMultiplier = -0.35;
% % %     ha.YRuler.TickLabelGapMultiplier = 0.001;
% % %     xlabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.12]);
% % %     ylabel('Probability','Units','normalized','Position',[-0.18 0.5]);
% % % 
% % %     % add legend
% % %     axes(panel_B_legend);
% % %     cla
% % %     hold on
% % %     plot([1 2], [2 2], 'Color', 'k', 'LineWidth', 2,'Clipping','off'); % long arm
% % %     plot([1 2], [1 1], 'Color', 'k', 'LineWidth', 1,'Clipping','off'); % shorty arm
% % %     xlim([0.5 2.5])
% % %     ylim([0 1])
% % %     text(2.5, 1, 'Short arm', 'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',7);
% % %     text(2.5, 2, 'Long arm',  'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',7);
% % %     set(gca,'visible','off');
% % % end
% % % axes(panel_B(1));
% % % text(-0.3,1.15, 'B', 'Units','normalized','FontWeight','bold');

%% panel B - plot (directions pooled)
axes(panel_B);
cla
hold on
text(-0.34,1.15, 'B', 'Units','normalized','FontWeight','bold');

x1=LS_field_ratio_all_long;
x2=LS_field_ratio_all;
x1(isnan(x1))=[];
x2(isnan(x2))=[];

nBinEdges = 9;
edges = logspace(0,log10(25),nBinEdges);
 % long arm
h1=histogram(x1);
h1.BinEdges = edges;
h1.Normalization = 'pdf';
h1.DisplayStyle = 'stairs';
h1.EdgeColor = long_arm_lc;
h1.LineWidth = long_arm_lw;
% entire_tunnel
h2=histogram(x2);
h2.BinEdges = edges;
h2.Normalization = 'pdf';
h2.DisplayStyle = 'stairs';
h2.EdgeColor = entire_tunnel_lc;
h2.LineWidth = entire_tunnel_lw;

[H,P_KS,KSSTAT] = kstest2(x1, x2, 'Tail','unequal');
fprintf('field size ratio long vs. all\t\t\t P_KS=%.2f D=%.2f n=%d m=%d\n',P_KS, KSSTAT, length(x1), length(x2));
text(1,1.05,sprintf('P_{KS} = %.2f',P_KS),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
% P_ranksum = ranksum(x1,x2);
% text(1,0.95,sprintf('P_{Wilc}=%.3f',P_ranksum),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);

ha=gca;
ha.YScale = yscale_opt;
% ha.YScale = 'linear';
ha.XScale = 'log';
ha.XLim = [0 27];
ha.XTick = [1 2 5 10 20];
ha.YLim = [3e-4 1.2];
ha.YTick = 10.^[-3 -2 -1];
ha.YTickLabel = {'10^{ -3}'; '10^{ -2}'; '10^{ -1}';};
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapOffset = -1;
xlabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.12]);
ylabel('Fraction of cells','Units','normalized','Position',[-0.23 0.5]);

%% add legend
axes(panel_B_legend);
cla
hold on
plot([1 2], [2 2], 'Color', long_arm_lc,      'LineWidth', long_arm_lw,      'Clipping', 'off'); % long arm
plot([1 2], [1 1], 'Color', entire_tunnel_lc, 'LineWidth', entire_tunnel_lw, 'Clipping', 'off'); % entire tunnel
xlim([0.5 2.5])
text(2.5, 2, 'Long arm',  'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',7);
text(2.5, 1, 'Entire tunnel', 'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',7);
set(gca,'visible','off');

%% print/save the figure
fig_name_out = fullfile(res_dir, sprintf('%s_yscale_%s',fig_name_str,yscale_opt));
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');








