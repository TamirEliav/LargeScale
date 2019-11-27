function cell_plot_figure_choose_examples(cell_ID)

%%
% clear
% clc
% cell_ID = 'b0034_d180314_TT4_SS03';

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','FR_map','fields','stats','FE','spikes','cluster_quality');
exp = exp_load_data(cell.details.exp_ID,'details','LM');
prm = PARAMS_GetAll();
dir_colors = prm.graphics.colors.flight_directions;
fig_name_str = sprintf('%s(%d)',cell.details.cell_ID, cell.details.cell_num);
dir_out = 'L:\Analysis\Results\cells\choose_examples';

%% create figure
figure_size_cm = [20 14];
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
annotation('textbox', [0.5 1 0 0], 'String', fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');
pause(0.2); % workaround to solve matlab automatically changing the axes positions...

%% create panels
panel_A_size_raster = [5 1];
panel_A_size_FR_map = [5 1];
panel_A = [];
panel_A_pos = [1.5 9];
panel_A(1) = axes('position', [panel_A_pos+[0 2.25]  panel_A_size_FR_map]);
panel_A(2) = axes('position', [panel_A_pos+[0 1   ]  panel_A_size_raster]);
panel_A(3) = axes('position', [panel_A_pos+[0 0   ]  panel_A_size_raster]);

panel_B = [];
panle_B_size = [0.8 1.2];
panle_B_pos = [3 4];
ch_offsets_x = [0 0.95 0 0.95 2];
ch_offsets_y = [1.3 1.3 0 0 0];
for ch=1:4+1
    pos_x = ch_offsets_x(ch);
    pos_y = ch_offsets_y(ch);
    panel_B(ch) = axes('position', [panle_B_pos+[pos_x pos_y] panle_B_size]);
%     text(0.5,0.5,"ch"+ch)
%     plot( squeeze(wvfs(:,1,:)) )
%     set(gca,'Visible','off');
end

panel_C = axes('position', [8 0 10 13]);


%% map+fields
axes(panel_A(1));
cla
maps=[cell.FR_map.all];
x = maps(1).bin_centers;
y = cat(1,maps.PSTH);
h=plot(x,y);
[h.Color] = disperse(dir_colors);
box off
h=gca;
h.TickDir = 'out';
h.XTick = [];
m = round(max(y(:)));
h.YTick = [0 m];
h.YLim = [0 m+1];
h.XLim = [0 200];

%% fields
for ii_dir=1:2
    fields = cell.fields{ii_dir};
    if isfield(fields,'in_low_speed_area')
        fields([fields.in_low_speed_area])=[];
    end
    for ii_field = 1:length(fields)
        dir_offsets = [-0.1 -0.17]+0.015;
        [xaf,yaf] = ds2nfu(fields(ii_field).edges_prc, repelem(dir_offsets(ii_dir)*range(h.YLim),2));
        annotation('line',xaf,yaf,'Linewidth', 2, 'Color', dir_colors{ii_dir});
    end
end
% details
cell_stats_str = {  sprintf('max=%.1fm', cell.stats.all.field_largest);...
                    sprintf('min=%.1fm', cell.stats.all.field_smallest);...
                    sprintf('ratio=%.1f', cell.stats.all.field_ratio_LS);...
                };
text(1,1.1, cell_stats_str,...
    'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',6);

%% rasters
FEs = [cell.FE];
for ii_dir=1:2
    axes(panel_A(ii_dir+1));
    cla
    FE = FEs{ii_dir};
    x = [FE.spikes_pos];
    [FE.number2] = disperse(1:length(FE));
    y = arrayfun(@(FE)(FE.number2*ones(1,FE.num_spikes)),FE,'UniformOutput',0);
    y = [y{:}];
    plot(x,y,'.','Color',dir_colors{ii_dir},'MarkerSize',0.05);
    box off
    h=gca;
    m = length(FE);
    h.YTick = [1 m];
    h.XLim = [0 200];
    h.YLim = [0 m+1];
    switch ii_dir
        case 1
            h.XTick = [];
            h.YTickLabel = {'',num2str(m)};
        case 2
            h.XTick = 0:50:200;
            h.XRuler.TickLabelGapOffset = -2;
            h.YTickLabel = {'1',num2str(m)};
            h.TickDir = 'out';
    end
end
axes(panel_A(3));
xlabel('Position (m)', 'Units','normalized','Position',[0.5 -0.35]);
ylabel('Flight no.',   'Units','normalized','Position',[-0.1 1]);
axes(panel_A(1));
ylabel({'Firing rate';'(Hz)'},   'Units','normalized','Position',[-0.07 0.42]);



%% waveforms
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
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
xlimits = [1 32];
m1 = min(wvfs,[],'all');
m2 = max(wvfs,[],'all');
ylimits = [m1 m2];
for ch=1:4
    axes(panel_B(ch));
    cla
    plot( squeeze(wvfs(:,ch,:)) )
    set(gca,'Visible','off');
    xlim(xlimits)
    ylim(ylimits)
end

% scale bars
axes(panel_B(5));
cla
hold on
plot([1 32],[0 0], 'k-', 'LineWidth', 2);
plot([1 1],[0 200], 'k-', 'LineWidth', 2);
set(gca,'Visible','off');
xlim(xlimits)
ylim(ylimits)

axes(panel_B(1));
text(0.1,1.3, sprintf('IsoDist: %.3g', cell.cluster_quality.Isolation_dis_mean),...
        'Units','normalized','HorizontalAlignment','left','FontSize',8);


    
%% table of stats
axes(panel_C)
axis off
% FR_map = cell.FR_map;
% fields = cell.fields;
stats_dir = cell.stats.dir;
stats_all = cell.stats.all;
stats_table_data =  {...
    'SI_bits_spike',    stats_dir.SI_bits_spike,   [], 'bit/spike';...
    'SI_bits_sec',      stats_dir.SI_bits_sec,     [], 'bit/sec';...
    'sparsity',         stats_dir.sparsity,        [], '';...
    'AC.width',         stats_dir.AC_width,        [], 'm';...
    'corr even/odd',    stats_dir.corr_odd_even,   [], '';...
    'corr begin/end',   stats_dir.corr_begin_end,  [], '';...
    'corr all/full',    stats_dir.corr_all_full,   [], '';...
    'map signif',       stats_dir.map_signif,      [], '';...
    'map shuff signif', stats_dir.map_signif_shuffle ,      [], '';...
%     '',                 [], [], [], '';...
    
    'num fields',       stats_dir.field_num,        stats_all.field_num,        '';...
    'Largest',          stats_dir.field_largest,    stats_all.field_largest,    'm';...
    'Smallest',         stats_dir.field_smallest,   stats_all.field_smallest,   'm';...
    'ratio L/S',        stats_dir.field_ratio_LS,   stats_all.field_ratio_LS,   '';...
    'CV field size',    stats_dir.field_CV,         stats_all.field_CV,         '';...
    '#spikes air',      stats_dir.spikes_num_air,   stats_all.spikes_num_air,   '';...
    '#spikes field',    stats_dir.spikes_num_field, stats_all.spikes_num_field, '';...
    '%spikes field',    stats_dir.spikes_prc_field, stats_all.spikes_prc_field, '%';...
	'',                 [], [], [],'';...
	'Iso Dist',         [], [], stats_all.IsoDist,       '';...
    'L-Ratio',          [], [], stats_all.L_Ratio,       '';...
    'meanFR all',       [], [], stats_all.meanFR_all,    'Hz';...
    'meanFR air',       [], [], stats_all.meanFR_flight, 'Hz';...
    'num flights',      stats_dir.num_full_flights, stats_all.num_full_flights, '';...
    'ClusterQuality',   [], [], cell.details.ClusterQuality,                    '';...
};
columnname =   {'Parameter', 'dir1', 'dir2', 'all', 'units'};
columnformat = {'char', 'numeric', 'numeric', 'numeric', 'char'}; 
t = uitable('Units','centimeters',...
            'Position', panel_C.Position,...
            'Data', stats_table_data,... 
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'RowName',[],...
            'ColumnWidth', {110 50 50 50 70},...
            'FontSize', 8);

% change numberic precision (manually...)
for ii_col = 1:length(t.ColumnFormat)
    if ~contains(t.ColumnFormat{ii_col},'numeric')
        continue
    end
    t.ColumnFormat{ii_col} = 'char';
    for ii_row = 1:size(t.Data,1)
        if t.Data{ii_row,ii_col} == round(t.Data{ii_row,ii_col})
            t.Data{ii_row,ii_col} = sprintf('%d', t.Data{ii_row,ii_col});
        else
            t.Data{ii_row,ii_col} = sprintf('%.2f', t.Data{ii_row,ii_col});
        end
    end
end
% set table background color
% first let's set the alternating colors
t.BackgroundColor(1:2:size(t.Data,1),:) = 1;
t.BackgroundColor(2:2:size(t.Data,1),:) = 0.94;
t.BackgroundColor(24,:) = interp1([0 1 2],[1 0 0 ; 1 1 0; 0 1 0], cell.details.ClusterQuality);




%% save figure
fig_filename = fullfile(dir_out, fig_name_str );
saveas(gcf,fig_filename,'tif')









%%









end