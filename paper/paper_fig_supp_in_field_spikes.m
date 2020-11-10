%% Large Scale - Fig. S5 - high spatial tuning with low baseline activity

%%
clear 
clc

%% plotting options (panel B)
plot_style = 'errorbars';
% plot_style = 'boxplot';
% plot_style = 'violin';
% plot_style = 'lines';
% plot_style = 'hist';
burst_ISI_thr = 6; % 6ms  as in Mizuseki 2011
% burst_ISI_thr = 7; % 7ms as Csicsvari, J., Hirase, H., Czurko, A., Buzsáki, G., 1998. Reliability and state dependence of pyramidal cell–interneuron synapses in the hippocampus: an ensemble approach in the behaving rat. Neuron 21, 179–189.


%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S5';
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
panel_B = axes('position', [ 8 20 2.9 3]);
panel_C = axes('position', [ 3 15 3 3]);
panel_D(1) = axes('position', [ 8 15 2.9 3]);
% panel_D(2) = axes('position', [ 13 15 2.9 3]);
panel_E = axes('position', [ 3 10 3 3]);
panel_F = axes('position', [ 8 10 2.9 3]);

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
% take only signif cells
signif = cat(1,cells.signif);
signif = arrayfun(@(x)(x.TF),signif);
signif = any(signif,2);
cells(~signif)=[];
clear signif

%% arrange data
in_field_spikes_prc = nan(length(cells),2);
InOutFieldFR = repelem(struct(),length(cells),2);
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    for ii_dir = 1:2
        if ~cell.signif(ii_dir).TF
            continue;
        end
        %% add cell num to each field
        for ii_field = 1:length(cell.fields{ii_dir})
            cells(ii_cell).fields{ii_dir}(ii_field).cell_num = cell.details.cell_num;
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

%% arrange pop data
signif = cat(1,cells.signif);
signif = arrayfun(@(x)(x.TF),signif);
stats = [cells.stats];
stats_all = [stats.all];
stats_dir = cat(1,stats.dir);
fields = cat(1,cells.fields);
fields = fields(signif);
fields = [fields{:}];
fields([fields.in_low_speed_area])=[];

%% calc burstyness stats (per field)
ISI_all = [];
spikes_in_burst_all = [];
spikes_in_burst_prc_all_fields = [];
for ii_field = 1:length(fields)
    field = fields(ii_field);
    ISI = diff(field.spikes_ts);
    ISI = ISI .* 1e-3; % convert to ms
    ISI_all = [ISI_all ISI];
    spikes_in_burst = zeros(1,length(ISI)+1);
    TF = ISI < burst_ISI_thr;
    spikes_in_burst([TF false]) = 1;
    spikes_in_burst([false TF]) = 1;
    spikes_in_burst_prc = sum(spikes_in_burst)/length(spikes_in_burst);
    spikes_in_burst_all = [spikes_in_burst_all spikes_in_burst];
    spikes_in_burst_prc_all_fields = [spikes_in_burst_prc_all_fields spikes_in_burst_prc];
end

%% calc burstyness stats (per cell)
subs = [fields.cell_num];
vals = spikes_in_burst_prc_all_fields;
spikes_in_burst_prc_all_cells = accumarray(subs',vals',[],@mean,nan);
spikes_in_burst_prc_all_cells(isnan(spikes_in_burst_prc_all_cells)) = [];

%% Panel A - percentage of out-of-field spikes
axes(panel_A);
cla
hold on
text(-0.35,1.15, 'A', 'Units','normalized','FontWeight','bold');

h = histogram(in_field_spikes_prc(:));
h.NumBins = 12;
h.FaceColor = 0.5*[1 1 1];
ha=gca;
ha.XLim = [0 101];
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
text(-0.35,1.15, 'B', 'Units','normalized','FontWeight','bold');
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
        h.BarWidth = 0.7;
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

%% within-field ISI hist
axes(panel_C);
cla('reset')
hold on
text(-0.35,1.15, 'C', 'Units','normalized','FontWeight','bold');

hh=histogram(ISI_all,[0:1:1000]);
hh.EdgeColor='k';
hh.FaceColor='k';
hh.Normalization='pdf';
% hl=xline(burst_ISI_thr);
% hl.Color='r';
% hl.LineWidth=1;

ha=gca;
ha.XLim = [0 300];
ha.YLim = [0 0.025];
ha.XTick = [0:100:1000];
ha.YTick = ha.YLim;
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.5;
ha.YRuler.TickLabelGapOffset = -1;
ha.YScale = 'linear';
xlabel({'In-field inter-spike interval (ms)'}, 'Units','normalized','Position',[0.5 -0.13]);
ylabel({'Probability';'density function'},'Units','normalized','Position',[-0.13 0.5]);

%% percentage of bursty spikes (per field)
axes(panel_D(1));
cla('reset')
hold on
text(-0.35,1.15, 'D', 'Units','normalized','FontWeight','bold');

hh=histogram(spikes_in_burst_prc_all_fields.*100);
hh.BinLimits = [0 100];
hh.NumBins = 100;
hh.Normalization='count';
hh.FaceColor = 0.5*[1 1 1];
ha=gca;
ha.XLim = [0 100];
ha.XTick = [0:25:100];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.1;
ha.YScale = 'log';
xlabel('In-field bursty spikes (%)', 'Units','normalized','Position',[0.5 -0.13]);
ylabel('No. of fields');
% text(0.5,0.85,"n="+length(spikes_in_burst_prc_all_fields),'Units','normalized');

%% percentage of bursty spikes (per cell)
% % % % % % % % % % axes(panel_D(2));
% % % % % % % % % % cla('reset')
% % % % % % % % % % hold on
% % % % % % % % % % % text(-0.35,1.15, 'D', 'Units','normalized','FontWeight','bold');
% % % % % % % % % % 
% % % % % % % % % % hh=histogram(spikes_in_burst_prc_all_cells.*100);
% % % % % % % % % % hh.BinLimits = [0 100];
% % % % % % % % % % hh.NumBins = 100;
% % % % % % % % % % hh.Normalization='count';
% % % % % % % % % % hh.FaceColor = 0.5*[1 1 1];
% % % % % % % % % % ha=gca;
% % % % % % % % % % ha.YScale = 'log';
% % % % % % % % % % ha.XLim = [0 100];
% % % % % % % % % % ha.XTick = [0:25:100];
% % % % % % % % % % ha.TickDir='out';
% % % % % % % % % % ha.TickLength = [0.03 0.03];
% % % % % % % % % % ha.XRuler.TickLabelGapMultiplier = -0.35;
% % % % % % % % % % ha.YRuler.TickLabelGapMultiplier = 0.1;
% % % % % % % % % % ha.YScale = 'linear';
% % % % % % % % % % xlabel('In-field bursty spikes (%)', 'Units','normalized','Position',[0.5 -0.13]);
% % % % % % % % % % ylabel('No. of cells');
% % % % % % % % % % text(0.85,0.85,"n="+length(spikes_in_burst_prc_all_cells),'Units','normalized');

%% mean firing rate hist
axes(panel_E);
cla('reset')
hold on
text(-0.35,1.15, 'E', 'Units','normalized','FontWeight','bold');

hh = histogram([stats_all.meanFR_flight]);
% hh.NumBins = 12;
hh.FaceColor = 0.5*[1 1 1];
ha=gca;
ha.XLim = [0 10];
ha.XTick = [0 5 10];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.1;
ha.YScale = 'linear';
xlabel('Mean firing rate (Hz)', 'Units','normalized','Position',[0.5 -0.13]);
ylabel('No. of cells');

%% field peak rate vs field size
axes(panel_F);
cla('reset')
hold on
text(-0.35,1.15, 'F', 'Units','normalized','FontWeight','bold');
x = [fields.width_prc];
y = [fields.peak];
plot(x,y,'.k');
% lm = fitlm(x,y);
% [rho, pval] = corr(x',y','type','pearson');
[rho, pval] = corr(x',y','type','spearman');
text(0.75,1.05, ['{\rho}' sprintf(' = %.2f',rho)] ,'units','normalized','FontSize',8);
if pval == 0
    text(0.75,0.9, 'P < 10^{ -300}' ,'units','normalized','FontSize',8);
else
    text(0.75,0.9, sprintf('P = %.2f',pval) ,'units','normalized','FontSize',8);
end

ha=gca;
ha.XLim = [0 35];
ha.XTick = [0 10 20 30];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.1;
ha.YScale = 'linear';
xlabel('Field size (m)', 'Units','normalized','Position',[0.5 -0.13]);
ylabel('Field peak rate (Hz)');



%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str+"_style_"+plot_style+"_burst_thr_"+burst_ISI_thr);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');



%%





