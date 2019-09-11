function cell_plot_map_fields(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID);
exp = exp_load_data(cell.details.exp_ID,'details','LM');
prm = PARAMS_GetAll();
dir_colors = prm.graphics.colors.flight_directions;

%% create figure
figure('Units','normalized','Position',[0 0 1 1]);
pnl = panel();
pnl.pack('h',[35 45 20])
pnl(1).pack('v',2)
pnl(1,1).pack('v',2)
pnl(1,2).pack('v',2)
pnl(2).pack('v',2)
pnl(3).pack('v',[80 20])
pnl.margin = [15 25 15 15];
pnl(1).margin = 10;
pnl(1).de.margin = 20;
pnl(1,1).margin = 10;
pnl(1,1).de.margin = 10;
pnl(1,2).margin = 10;
pnl(1,2).de.margin = 5;


%% FR map + fields + LM
for ii_dir = 1:2
    pnl(1,1,ii_dir).select();hold on
    
    % FR map
    plot(cell.FR_map(ii_dir).all.bin_centers, cell.FR_map(ii_dir).all.PSTH,...
        'Color', dir_colors{ii_dir}, 'LineWidth',1.5);
    
    % fields
    fields = cell.fields{ii_dir};
    for ii_field = 1:length(fields)
        field = fields(ii_field);
        if field.in_low_speed_area
            c = 'r';
        else
            c = 'k';
        end
        plot(field.loc, field.peak, '*', 'MarkerSize', 4, 'Color', c)
        plot(field.edges_href, repelem(prm.fields.width_href * field.peak,2), 'Color', c)
        text(field.loc, field.peak,...
            sprintf('{\\color{magenta}%d}\n{\\color{red}%2.0f%%}\n%2.2f\n%2.2f',....
            length(field.spikes_ts),...
            field.num_flights_with_spikes_prc*100,...
            field.width_href,...
            field.width_prc),...
            'FontSize',7,'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
    
    % plot Landmarks
    plot_LM(exp.LM);
    
    % labels
%     xlabel('Position (m)')
    ylabel('Firing rate (Hz)')
    set(gca,'TickDir','out');
    if ii_dir == 1
        set(gca,'XTickLabel',[]);
    end
end

%% trajectory + spikes
ti = exp_get_sessions_ti(cell.details.exp_ID, 'Behave');
for ii_dir = 1:2
    pnl(1,2,ii_dir).select();hold on
    FE = cell.FE{ii_dir};
    plot([FE.pos],       [FE.ts],        '.', 'Color', 0.9.*[1 1 1])
    plot([FE.spikes_pos],[FE.spikes_ts], '.', 'Color', dir_colors{ii_dir})
    xlimits = get(gca, 'xlim');
    if ~isempty(cell.details.stable_ts)
        plot(repmat(xlimits,[2 1])', repmat(cell.details.stable_ts, [ 2 1] ),'--', 'Color','k');
    end
    ylim(ti)
    rescale_plot_data('y',[1e-6/60 ti(1)])
    
    % plot Landmarks
    plot_LM(exp.LM);
    
    % labels
%     xlabel('Position (m)')
    ylabel('Time (min)')
    set(gca,'TickDir','out');
    if ii_dir == 1
        set(gca,'XTickLabel',[]);
    end
    if ii_dir == 2
        xlabel('Position (m)')
    end
end

%% link position axes
linkaxes(pnl(1).de.axis, 'x');

%% add arrows for directions
arrow_pos_X = [.01 .05;
               .05 .01]+0.04;
arrow_pos_Y = [repelem(0.825, 2);
               repelem(0.815, 2)]+0.15;
for ii_dir = 1:2
    c = dir_colors{ii_dir};
    ah=annotation('arrow',arrow_pos_X(ii_dir,:),arrow_pos_Y(ii_dir,:),'Color',c);
    ah.LineWidth = 1.5;
    set(ah,'HeadStyle','cback1','HeadWidth',5);
end


%% spikes waveforms (per field)
cell = fields_add_spikes_waveforms(cell);
for ii_dir = 1:2
    %%
    fields = cell.fields{ii_dir};
    if length(fields)==0
        continue;
    end
    
    %% create axes for fields spikes waveforms
    pnl(2,ii_dir).pack('h',[25 75]);
    max_col = 5;
    max_row = 4;
    max_fields_to_plot = max_col * max_row;
    ncol = min(length(fields),max_col);
    nrow = min(ceil( length(fields) / ncol ), max_row);
    pnl(2,ii_dir,2).pack(nrow,ncol);
    pnl(2,ii_dir,2).select('all');
    pnl(2,ii_dir,2).de.margin = 5;
    fields_axes = pnl(2,ii_dir,2).de;
    fields_axes(cellfun(@(x)(isempty(x.axis)),fields_axes)) = [];
    if length(fields) > max_fields_to_plot
        warning(sprintf('Too many fields!!! plotting only %d out of %d fields',max_fields_to_plot,length(fields)));
    end
    
    %% compare mean waveforms of the strongest channel
    pnl(2,ii_dir,1).select(); hold on;
    all_fields_spikes_wvfrms = cat(3,fields.spikes_wvfrm);
    [~,ch_max] = max(max(mean(all_fields_spikes_wvfrms,3)),[],2);
    plot( mean(all_fields_spikes_wvfrms(:,ch_max,:),3) , 'k', 'LineWidth',1.5);
    for ii_field = 1:length(fields)
        field = fields(ii_field);
        field_spikes_wvfrms = [field.spikes_wvfrm];
        plot( mean(field_spikes_wvfrms(:,ch_max,:),3));
    end
    xlim([1 32])
    legend_str = cellfun(@(x)(sprintf('#%d',x)),num2cell(1:length(fields)),'UniformOutput',false);
    legend_str = {'all' legend_str{:}};
    leg_max_row = 8;
    leg_ncol = ceil(length(legend_str)/leg_max_row);
    hleg = columnlegend(leg_ncol, legend_str, 'location', 'NorthEast','boxoff');
    hax = gca;
    hleg.Position([1:2]) = hax.Position([1:2])+[0.06 0.14];
    hleg.Position([3:4]) = [0.03 0.2];
    title('max ch. all fields')
    xlabel('Samples')
    ylabel('\muV')
    
    %% plot each field spikes seperaetly
    for ii_field = 1:min(length(fields),max_fields_to_plot)
        field = fields(ii_field);
%         ax = pnl(2,ii_dir,2).de(ii_field);
        ax = fields_axes{ii_field};
        ax.select(); hold on;
        field_spikes_wvfrms = [field.spikes_wvfrm];
        plot(mean(field_spikes_wvfrms,3))
        xlim([1 32])
        title(sprintf('#%d',ii_field))
        axis off
    end
    % remove extra axes
    for ii_ax = (length(fields)+1) : length(fields_axes)
        ax = fields_axes{ii_ax};
        ax.select(); hold on;
        axis off
    end
    
end
if length(cell.fields{2})>0 & length(cell.fields{2})>0
    linkaxes(pnl(2).de.axis, 'xy');
end
        
%% recording stability
pnl(3,2).select();hold on

% plot avg FR in large chuncks
bar(cell.RecStability.bin_centers, cell.RecStability.FR);
ylimits = get(gca,'ylim');
for ii_session = 1:length(exp.details.session_names)
    sname = exp.details.session_names{ii_session};
    ti = exp.details.session_ts(ii_session,:);
    plot(repelem(ti(1),2), ylimits,'-m');
    plot(repelem(ti(2),2), ylimits,'-m')
    text(mean(ti), ylimits(2)+0.05*diff(ylimits), sname,'HorizontalAlignment','center','FontWeight','bold');
end
xlimits = exp.details.session_ts([1 end]);
xlimits = xlimits + [-1 1].*range(xlimits)*0.025;
xlim(xlimits);
% add stability ts
ylimits = get(gca,'ylim');
if ~isempty(cell.details.stable_ts)
    plot( repmat(cell.details.stable_ts, [2 1]), repmat(ylimits,[2 1])', '--', 'Color','k');
end
% align time to behave session (time zero)
behave_ti = exp_get_sessions_ti(exp.details.exp_ID,'Behave');
rescale_plot_data('x',[1e-6/60 behave_ti(1)])
xlimits = get(gca,'xlim'); % workaround
xlabel('Time (min)')
ylabel('Firing rate (Hz)')

% plot spikes peak vs. time
wvfms=[cell.spikes.waveforms];
[~,max_ch_IX] = max(max(squeeze(mean(wvfms,3))));
spikes_peaks = squeeze(range(wvfms(:,max_ch_IX,:),1));
yyaxis right
hold on
plot(cell.spikes.ts, spikes_peaks,'.','MarkerSize',4, 'Color',[1 0 0]);
ylim([0 max(spikes_peaks)*1.1])
rescale_plot_data('x',[1e-6/60 behave_ti(1)])
ylabel('Spikes size(\muV)')

% set xlim again (workaround)
xlim(xlimits);


%% results table
if isfield(cell,'stats')
pnl(3,1).select();hold on
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
t = uitable('Units','normalized',...
            'Position', pnl(3,1).axis.Position,...
            'Data', stats_table_data,... 
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'RowName',[],...
            'ColumnWidth', {110 50 50 50 70},...
            'FontSize', 12);

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
end
%% set table background color
% first lets set the alternating colors
t.BackgroundColor(1:2:size(t.Data,1),:) = 1;
t.BackgroundColor(2:2:size(t.Data,1),:) = 0.94;
t.BackgroundColor(24,:) = interp1([0 1 2],[1 0 0 ; 1 1 0; 0 1 0], cell.details.ClusterQuality);

%% save figure
h=pnl.title(sprintf('%s (%d)',cell.details.cell_ID,cell.details.cell_num)); h.FontSize=16; h.Interpreter='none';h.Position=[0.5 1.03];
fig_filename = fullfile('L:\Analysis\Results\cells\figures', sprintf('%s(%d)_map_fields',cell.details.cell_ID,cell.details.cell_num));
saveas(gcf,fig_filename,'tif')
% saveas(gcf,fig_filename,'fig')



end



%% local function to plot lines for landmarks
function plot_LM(LM)

ylimits = get(gca,'ylim');
for ii_LM=1:length(LM)
    x = LM(ii_LM).pos_proj;
    plot(repelem(x,2) , ylimits, '-', 'color', 0.9.*[1 1 1], 'LineWidth',0.5)
%     text(x, ylimits(2)+0.02*diff(ylimits), LM(ii_LM).name, 'Rotation', 45, 'FontSize',4)
end

end






%%
