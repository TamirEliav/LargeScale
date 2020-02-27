%% Large Scale - Fig. S4 - high spatial tuning with low baseline activity

%%
clear 
clc

%% plotting options (panel B)
plot_style = 'errorbars';
% plot_style = 'boxplot';
% plot_style = 'violin';
% plot_style = 'lines';
% plot_style = 'hist';


%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S4_new';
fig_caption_str = 'In-field spikes';
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
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');
pause(0.2); % workaround to solve matlab automatically changing the axes positions...

% create panels
panel_A = axes('position', [ 3 20 3 3]);
panel_B = axes('position', [ 7.5 20 2.5 3]);


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
cells = cellfun(@(c)(cell_load_data(c,'details','meanFR')), {cells.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];
cells_details = [cells.details];
cells_ID = {cells_details.cell_ID};
meanFR = [cells.meanFR];
cells_ID([meanFR.all]>prm.inclusion.interneuron_FR_thr)=[];
clear cells stats cells_details cells_t
cells = cellfun(@(c)(cell_load_data(c,'details','stats','meanFR','stats','inclusion','signif','fields','FR_map','FE')), cells_ID, 'UniformOutput',0);
cells = [cells{:}];

%% arrange data
in_field_spikes_prc = nan(length(cells),2);
InOutFieldFR = repelem(struct(),length(cells),2);
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    for ii_dir = 1:2
        if ~cell.signif(ii_dir).TF
            continue;
        end
        %% in-field spikes percentage
        fields = cell.fields{ii_dir};
        fields([fields.in_low_speed_area])=[];
        fields_spikes_ts = [fields.spikes_ts];
        FE_spikes_ts = [cell.FE{ii_dir}.spikes_ts];
        FE_spikes_pos = [cell.FE{ii_dir}.spikes_pos];
        invalid_IX =( FE_spikes_pos < prm.fields.valid_speed_pos(1) | ...
                      FE_spikes_pos > prm.fields.valid_speed_pos(2) );
        FE_spikes_ts(invalid_IX)=[];
        FE_spikes_pos(invalid_IX)=[];
        in_field_spikes = ismember(FE_spikes_ts, fields_spikes_ts);
        in_field_spikes_prc(ii_cell,ii_dir) = 100 * sum(in_field_spikes) / length(in_field_spikes);
        %% in/out of field FR
        FR_map = cell.FR_map(ii_dir).all;
        fields_edges = cat(1,fields.edges_prc);
%         fields_edges = cat(1,fields.edges_href);
        x = FR_map.bin_centers;
        in_field_pos_IX  =  any( x > fields_edges(:,1) & x<fields_edges(:,2),1) & ...
                               ( x > prm.fields.valid_speed_pos(1)) & (x < prm.fields.valid_speed_pos(2));
        out_field_pos_IX = ~any( x > fields_edges(:,1) & x<fields_edges(:,2),1) & ...
                               ( x > prm.fields.valid_speed_pos(1)) & (x < prm.fields.valid_speed_pos(2));
        in_field_FR = FR_map.PSTH(in_field_pos_IX);
        out_field_FR = FR_map.PSTH(out_field_pos_IX);
        in_field_FR_mean = mean(in_field_FR);
        in_field_FR_median = median(in_field_FR);
        out_field_FR_mean = mean(out_field_FR);
        out_field_FR_median = median(out_field_FR);
        InOutFieldFR(ii_cell,ii_dir).in_field_FR = in_field_FR;
        InOutFieldFR(ii_cell,ii_dir).out_field_FR = out_field_FR;
        InOutFieldFR(ii_cell,ii_dir).in_field_FR_mean = in_field_FR_mean;
        InOutFieldFR(ii_cell,ii_dir).in_field_FR_median = in_field_FR_median;
        InOutFieldFR(ii_cell,ii_dir).out_field_FR_mean = out_field_FR_mean;
        InOutFieldFR(ii_cell,ii_dir).out_field_FR_median = out_field_FR_median;
        
    end
end

%% Panel A - percentage of out-of-field spikes
axes(panel_A);
cla
hold on
text(-0.35,1.15, 'A', 'Units','normalized','FontWeight','bold');

h = histogram(in_field_spikes_prc(:));
h.NumBins = 12;
h.FaceColor = 0.5*[1 1 1];
ha=gca;
ha.XLim = [0 100];
ha.XTick = [0 50 100];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.1;
ha.YScale = 'linear';
xlabel('In-field spikes (%)', 'Units','normalized','Position',[0.5 -0.13]);
ylabel('No. of cells')

%% Panel B - percentage of out-of-field spikes
axes(panel_B);
% cla RESET
cla(gca,'reset')
hold on
text(-0.2,1.15, 'B', 'Units','normalized','FontWeight','bold');
box off

x1 = [InOutFieldFR.in_field_FR_median];
x2 = [InOutFieldFR.out_field_FR_median];
values = [x1 x2]';
grps = [ones(size(x1)) 2*ones(size(x2))]';

% plot
switch plot_style 
    case 'errorbars'
        x = unique(grps); % in/out
        y = accumarray(grps, values, [], @nanmean, nan);
        err = accumarray(grps, values, [], @nansem, nan);
        fprintf('in-field mean FR = %.2f +- %.2f\n', y(1), err(1));
        fprintf('out-of-field mean FR = %.2f +- %.2f\n', y(2), err(2));
        h=bar(x,y);
        h.FaceColor = 0.65*[1 1 1];
        h=errorbar(x,y,err);
        h.LineStyle='none';
        h.Color = 'k';
        h.LineWidth = 1.1;
    case 'boxplot'
        boxplot(values, grps);
    case 'violin'
        rng(0);
        h=violinplot(values, grps);
%         [h.ViolinColor] = disperse( {'k';'k'} );
%         hs = [h.ScatterPlot];
%         [hs.Marker] = disperse( {'.';'.'} );
    case 'lines'
        plot([x1;x2],'.-')
    case 'hist'
        h=histogram(x1-x2);
        h.BinLimits = [-20 20];
        h.NumBins = 20;
end

ha=gca;
switch plot_style
    case {'errorbars', 'boxplot', 'violin', 'lines'}
        ha.XLim = [0.5 2.5];
        ha.XTick = [1 2];
        ha.XTickLabel = {'In-field';'Out-of-field'};
        ha.TickDir='out';
        ha.TickLength = [0.025 0.025];
        % ha.XRuler.TickLabelGapMultiplier = -0.35;
        % ha.YRuler.TickLabelGapMultiplier = 0.1;
        ylabel('Firing rate (Hz)')
    case 'hist'
        xlabel({'{\Delta}Firing rate (Hz)';'In vs. Out of field'})
        ylabel('No. of cells')
end

%% print/save the figure
fig_name_out = fullfile(res_dir, [fig_name_str '_style_' plot_style]);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');



%%





