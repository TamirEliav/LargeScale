%% Large Scale - fig. supp - spike sorting control

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_supp_spike_sorting_control';
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
panle_B_size = [1.5 1.5];
ch_offsets_x = [0 0 2 2];
ch_offsets_y = [0 2 0 2];
panel_A = axes('position', [3 21 4 4]);
for ii=1:3
    for jj=1:3
        for ch=1:4
            pos_x = 3 + (ii-1)*5 + ch_offsets_x(ch);
            pos_y = 5 + (jj-1)*5 + ch_offsets_y(ch);
            panel_B(ii,jj,ch) = axes('position', [pos_x pos_y panle_B_size]);
        end
    end
end
panel_B = flipdim(panel_B,2);
panel_B = reshape(panel_B,[9 4]);

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
cells = cellfun(@(c)(cell_load_data(c,'details','stats')), {cells.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];
stats = [cells.stats];
stats = [stats.all];
cells([stats.meanFR_all]>prm.inclusion.interneuron_FR_thr)=[];
% choose only signif cells (in any direction)
cells_details = [cells.details];
cells = cellfun(@(c)(cell_load_data(c,'details','signif')), {cells_details.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];
signif = cat(1,cells.signif);
signif = arrayfun(@(x)(x.TF), signif);
signif = any(signif,2);
cells(~signif)=[];

cells = cellfun(@(c)(cell_load_data(c,'details','stats','cluster_quality')), {cells_details.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];
CQ = [cells.cluster_quality];
stats = [cells.stats];
stats_all = [stats.all];


%% panel A - field size ratio vs. Isolation index
axes(panel_A);
cla
hold on
text(-0.27,1.1, 'A', 'Units','normalized','FontWeight','bold');

x = [CQ.Isolation_dis_mean];
y = [stats_all.field_ratio_LS];
plot(x, y, '.');
xlim([0 100])
[r,pval_r]     = corr(x',y','rows','pairwise','type','Pearson');
[rho,pval_rho] = corr(x',y','rows','pairwise','type','Spearman');
text(1,1,   {sprintf('r=%.2f',r);sprintf('P=%.2f',pval_r)}, ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
% text(1,0.8, {sprintf('rho=%.2f',rho);sprintf('P=%.2f',pval_rho)}, ...
%         'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    
xlabel('Isolation index','Units','normalized','Position',[0.5 -0.15]);
ylabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[-0.15 0.5]);


%%
cell_examples = {
433; 56; 51;
609; 67; 477;
419; 628; 337;
};
for ii_cell = 1:length(cell_examples)
    %% load cell data
    cell_ID = cell_examples{ii_cell};
    cell = cell_load_data(cell_ID,'details','fields','stats','spikes','cluster_quality');
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

    %% plot
    for ch=1:4
        axes(panel_B(ii_cell,ch));
        cla
        wvfs = cat(3,fields.spikes_wvfs_avg);
        plot( squeeze(wvfs(:,ch,:)) )
    end
    axes(panel_B(ii_cell,2));
    text(1.2,1.2, sprintf('cell %d (IsoDist: %.3g)',cell_ID, cell.cluster_quality.Isolation_dis_mean),...
            'Units','normalized','HorizontalAlignment','center','FontSize',8);
end
axes(panel_B(1,2));
text(-0.7,1.35, 'B', 'Units','normalized','FontWeight','bold');



%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');


