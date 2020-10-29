%% Large Scale - fig. supp - spike sorting control

%%
clear 
clc

%% options
% corr_type = 'pearson';
corr_type = 'spearman';

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S9';
fig_caption_str = 'Spike sorting control';
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
panel_A = axes('position', [3  20 4 4]);
panel_B = axes('position', [10 20 4 4]);
panel_C = [];
panle_C_size = [0.8 1.2];
% ch_offsets_x = [0:4].*0.95;
% ch_offsets_x(end) = ch_offsets_x(end-1);
ch_offsets_x = [0 0.95 0 0.95 2];
ch_offsets_y = [1.3 1.3 0 0 0];
for ii=1:3
    for jj=1:3
        for ch=1:4+1
            pos_x = 3 + (ii-1)*5 + ch_offsets_x(ch);
            pos_y = 6 + (jj-1)*4 + ch_offsets_y(ch);
            panel_C(ii,jj,ch) = axes('position', [pos_x pos_y panle_C_size]);
%             text(0.5,0.5,"ch"+ch)
%             plot( squeeze(wvfs(:,1,:)) )
%             set(gca,'Visible','off');
        end
    end
end
panel_C = flipdim(panel_C,2);
panel_C = reshape(panel_C,[9 5]);


%% load population data
% get list of significant cells (at least in one direction)
% ---------------------------------------------------------
prm = PARAMS_GetAll();
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
cells_details = [cells.details];

% load chosen cells data
cells = cellfun(@(c)(cell_load_data(c,'details','stats','cluster_quality','signif')), {cells_details.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];
CQ = [cells.cluster_quality];
stats = [cells.stats];
stats_all = [stats.all];
stats_dir = cat(1,stats.dir);


%% panel A - No. of fields vs. Isolation index
axes(panel_A);
cla
hold on
text(-0.27,1.1, 'A', 'Units','normalized','FontWeight','bold');

x = [CQ.Isolation_dis_mean];
x = [x;x]';
% y = [stats_all.field_num];
signif_per_dir = arrayfun(@(x)(x.TF), cat(1,cells.signif));
y = arrayfun(@(x)(x.field_num), stats_dir);
y(~signif_per_dir) = nan;

x = x(:)';
y = y(:)';

plot(x, y, '.k');
xlim([0 100])
[r,pval_r]     = corr(x',y','rows','pairwise','type','Pearson');
[rho,pval_rho] = corr(x',y','rows','pairwise','type','Spearman');
fprintf('num_field vs. IsoDist corr: r=%.2f, P=%.2f\n',r,pval_r);
fprintf('num_field vs. IsoDist corr: rho=%.2f, P=%.2f\n',rho,pval_rho);
switch corr_type
    case 'pearson'
        text(1,1,   {sprintf('r = %.3f',r);sprintf('P = %.2f',pval_r)}, ...
            'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    case 'spearman'
        text(0.94,1, {['{\rho}' sprintf(' = %.3f',rho)];sprintf('P = %.3f',pval_rho)}, ...
            'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
end

xlabel('Isolation distance','Units','normalized','Position',[0.5 -0.13]);
ylabel({'No. of fields per direction'},'Units','normalized','Position',[-0.15 0.5]);

fprintf('panel A number of cells=%d\n',sum(~isnan(y)))

%% panel B - field size ratio vs. Isolation index
axes(panel_B);
cla
hold on
text(-0.35,1.1, 'B', 'Units','normalized','FontWeight','bold');

x = [CQ.Isolation_dis_mean];
y = [stats_all.field_ratio_LS];
plot(x, y, '.k');
xlim([0 100])
[r,pval_r]     = corr(x',y','rows','pairwise','type','Pearson');
[rho,pval_rho] = corr(x',y','rows','pairwise','type','Spearman');
fprintf('field size ratio vs. IsoDist corr: r=%.2f, P=%.2f\n',r,pval_r);
fprintf('field size ratio vs. IsoDist corr: rho=%.2f, P=%.2f\n',rho,pval_rho);
switch corr_type
    case 'pearson'
        text(1,1,   {sprintf('r = %.3f',r);sprintf('P = %.2f',pval_r)}, ...
            'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    case 'spearman'
        text(0.94,1, {['{\rho}' sprintf(' = %.3f',rho)];sprintf('P = %.3f',pval_rho)}, ...
            'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
end
xlabel('Isolation distance','Units','normalized','Position',[0.5 -0.13]);
ylabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[-0.15 0.5]);

fprintf('panel B number of cells=%d\n',sum(~isnan(y)))

%% panel C - examples: waveforms of spikes from different fields
cell_examples = {
433;  56;  51;
609; 419; 477;
 57; 628; 337;
};

loggerType_2_fs = containers.Map(...
    {'SpikeLog-16';'MiniBatLog-16'},...
    [29296.875      31250]            );

for ii_cell = 1:length(cell_examples)
    %% load cell data
    cell_ID = cell_examples{ii_cell};
    cell = cell_load_data(cell_ID,'details','fields','stats','spikes','cluster_quality');
    exp = exp_load_data(cell.details.exp_ID,'details');
    c = prm.graphics.colors.flight_directions;
    
    %% arrange spikes data by fields
    fields = cell.fields;
    fields_all = [];
    for ii_dir = 1:2
        if isempty(fields{ii_dir})
            continue
        end
        fields_to_add = fields{ii_dir};
        fields_all = [fields_all fields_to_add];
    end
    fields_all([fields_all.in_low_speed_area]) = [];
    fields = fields_all;
    
    [~,maxch] = max(max(mean(cell.spikes.waveforms,3),[],1));
    % add spikes waveforms for field struct
    for ii_field = 1:length(fields)
        field = fields(ii_field);
        IX = ismember(cell.spikes.ts,field.spikes_ts);
        fields(ii_field).spikes_wvfs = cell.spikes.waveforms(:,:,IX);
        fields(ii_field).spikes_wvfs_avg = mean(cell.spikes.waveforms(:,:,IX),3);
        fields(ii_field).spikes_wvfs_avg_max = mean(cell.spikes.waveforms(:,maxch,IX),3);
    end

    %% plot waveforms per field (4 channles)
    wvfs = cat(3,fields.spikes_wvfs_avg);
    fs = loggerType_2_fs(exp.details.NeuralLoggerType);
    dt = 1/fs;
    t = (0:31).*dt;
    t = t*1e3; % convert to ms
    xlimits = [0 1];
    m1 = min(wvfs,[],'all');
    m2 = max(wvfs,[],'all');
    ylimits = [m1 m2];
    for ch=1:4
        axes(panel_C(ii_cell,ch));
        cla
        plot(t, squeeze(wvfs(:,ch,:)) )
        set(gca,'Visible','off');
        xlim(xlimits)
        ylim(ylimits)
    end
    
    % scale bars
    axes(panel_C(ii_cell,5));
    cla
    hold on
    
    scale_ms = 1;
    scale_uV = 100;
    if scale_uV > ylimits(2)
        error('largest channel amplitude is smaller than the chosen scale bar!')
    end
    plot([0 1]*scale_ms,[0 0], 'k-', 'LineWidth', 1.5);
    plot([1 1]*scale_ms,[0 scale_uV], 'k-', 'LineWidth', 1.5);
    text(0.5*scale_ms,  -0.18*range(ylimits),  sprintf('%dms',scale_ms),      'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize', 8);
    text(1.1*scale_ms,  0.5*scale_uV,         sprintf('%d  {\\mu}V',scale_uV), 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize', 8);
    set(gca,'Visible','off');
    xlim(xlimits)
    ylim(ylimits)
    
    axes(panel_C(ii_cell,1));
    text(0.1,1.8, sprintf('Cell %d',ii_cell),...
            'Units','normalized','HorizontalAlignment','left','FontSize',8);
    text(0.1,1.5, sprintf('Isolation distance = %.3g', cell.cluster_quality.Isolation_dis_mean),...
        'Units','normalized','HorizontalAlignment','left','FontSize',8);
	text(0.1,1.2, sprintf('N fields = %d', size(wvfs,3) ),...
            'Units','normalized','HorizontalAlignment','left','FontSize',8);
end
axes(panel_C(1,2));
text(-2.5,2, 'C', 'Units','normalized','FontWeight','bold');


%% print/save the figure
fig_name_out = fullfile(res_dir, sprintf('%s__corr_type_%s',fig_name_str,corr_type));
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');


