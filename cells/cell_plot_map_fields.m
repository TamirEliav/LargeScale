function cell_plot_map_fields(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','FE','FR_map','Ipos','fields','RecStability','spikes','signif');
exp = exp_load_data(cell.details.exp_ID,'details');
prm = PARAMS_GetAll();
dir_colors = prm.graphics.colors.flight_directions;

%% create figure
figure('Units','normalized','Position',[0 0 1 1]);
pnl = panel();
pnl.pack('h',[30 50 20])
pnl(1).pack('v',2)
pnl(1,1).pack('v',2)
pnl(1,2).pack('v',2)
pnl(2).pack('v',2)
pnl(3).pack('v',[60 40])
pnl.margin = [15 25 10 15];
pnl(1).margin = 10;
pnl(1).de.margin = 10;
% pnl(1,1).margin = 10;
pnl(1,1).de.margin = 7;
pnl(1,2).margin = 20;
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
        plot(field.loc, field.peak, 'k*', 'MarkerSize', 4)
        plot(field.edges_href, repelem(prm.fields.width_href * field.peak,2), 'k')
        text(field.loc, field.peak, sprintf('%2.2f\n%2.2f',field.width_href,field.width_prc),...
            'FontSize',7,'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
    
    % Landmarks (TODO: add!)
    
    % labels
    xlabel('Position (m)')
    ylabel('Firing rate (Hz)')
end

%% trajectory + spikes
ti = exp_get_sessions_ti(cell.details.exp_ID, 'Behave');
for ii_dir = 1:2
    pnl(1,2,ii_dir).select();hold on
    FE = cell.FE{ii_dir};
    plot([FE.pos],       [FE.ts],        '.', 'Color', 0.9.*[1 1 1])
    plot([FE.spikes_pos],[FE.spikes_ts], '.', 'Color', dir_colors{ii_dir})
    ylim(ti)
    rescale_plot_data('y',[1e-6/60 ti(1)])
    
    % labels
    xlabel('Position (m)')
    ylabel('Time (min)')
end

%% link position axes
linkaxes(pnl(1).de.axis, 'x');

%% spikes waveforms (per field)
cell = fields_add_spikes_waveforms(cell);
for ii_dir = 1:2
    %%
    fields = cell.fields{ii_dir};
    
    %% create axes for fields spikes waveforms
    pnl(2,ii_dir).pack('h',[25 75]);
    ncol = min(length(fields),8);
    nrow = ceil( length(fields) / ncol );
    pnl(2,ii_dir,2).pack(nrow,ncol);
    pnl(2,ii_dir,2).select('all');
    fields_axes = pnl(2,ii_dir,2).de;
    fields_axes(cellfun(@(x)(isempty(x.axis)),fields_axes)) = [];
    
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
    legend_str = cellfun(@(x)(sprintf('#%d',x)),num2cell(1:length(fields)),'UniformOutput',false);
    legend_str = {'all fields' legend_str{:}};
    legend(legend_str )
    title('max ch. all fields')
    xlabel('Samples')
    ylabel('\muV')
    
    %% plot each field spikes seperaetly
    for ii_field = 1:length(fields)
        field = fields(ii_field);
%         ax = pnl(2,ii_dir,2).de(ii_field);
        ax = fields_axes{ii_field};
        ax.select(); hold on;
        field_spikes_wvfrms = [field.spikes_wvfrm];
        plot(mean(field_spikes_wvfrms,3))
        title(sprintf('#%d',ii_field))
    end
end
linkaxes(pnl(2).de.axis, 'xy');

%% results table
pnl(3,1).select();hold on
axis off
stats_table_data =  {...
    'SI_bits_spike',    cell.FR_map(1).all.SI_bits_spike, cell.FR_map(2).all.SI_bits_spike, 'bit/spike';...
    'SI_bits_sec',      cell.FR_map(1).all.SI_bits_sec, cell.FR_map(2).all.SI_bits_sec, 'bit/sec';...
    'sparsity',         cell.FR_map(1).all.sparsity, cell.FR_map(2).all.sparsity, '';...
    'AC.width',         cell.FR_map(1).all.AC.width, cell.FR_map(2).all.AC.width, 'm';...
    'corr even/odd',    cell.FR_map(1).corr_odd_even.rho, cell.FR_map(2).corr_odd_even.rho, '';...
    'corr begin/end',   cell.FR_map(1).corr_begin_end.rho, cell.FR_map(2).corr_begin_end.rho, '';...
    'corr all/full',    cell.FR_map(1).corr_all_full.rho, cell.FR_map(2).corr_all_full.rho, '';...
    'map signif',       cell.signif(1).TF, cell.signif(2).TF, '';...
    ' ',                [], [], ' ';...
    
    'num fields',       length(cell.fields{1}), length(cell.fields{2}), '';...
    'Largest',          max([cell.fields{1}.width_prc]), max([cell.fields{2}.width_prc]), '';...
    'Smallest',         min([cell.fields{1}.width_prc]), min([cell.fields{2}.width_prc]), '';...
    'ratio L/S',        max([cell.fields{1}.width_prc])/min([cell.fields{1}.width_prc]), ...    
                        max([cell.fields{2}.width_prc])/min([cell.fields{2}.width_prc]), '';...
    'CV field size',    std([cell.fields{1}.width_prc])/mean([cell.fields{1}.width_prc]), ...    
                        std([cell.fields{2}.width_prc])/mean([cell.fields{2}.width_prc]), '';...
	' ',                [], [], ' ';...
	'Iso Dist.',        cell.spikes.Isolation_dis, [], '';...
    'L-Ratio.',         cell.spikes.L_Ratio, [], '';...
                        
                        
};
columnname =   {'Parameter', 'dir1', 'dir2', 'units'};
columnformat = {'char', 'numeric', 'numeric'}; 
t = uitable('Units','normalized',...
            'Position', pnl(3,1).axis.Position,...
            'Data', stats_table_data,... 
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'RowName',[]);

%% recording stability
pnl(3,2).select();hold on
behave_ti = exp_get_sessions_ti(exp.details.exp_ID,'Behave');
bar(cell.RecStability.bin_centers, cell.RecStability.FR)
ylimits = get(gca,'ylim');
for ii_session = 1:length(exp.details.session_names)
    sname = exp.details.session_names{ii_session};
    ti = exp.details.session_ts(ii_session,:);
    plot(repelem(ti(1),2), ylimits,'--m');
    plot(repelem(ti(2),2), ylimits,'--m')
    text(mean(ti), ylimits(2), sname,'HorizontalAlignment','center','FontWeight','bold');
end
rescale_plot_data('x',[1e-6/60 behave_ti(1)])
xlabel('Time (min)')
ylabel('Firing rate')

%% save figure
h=pnl.title(cell_ID); h.FontSize=16; h.Interpreter='none';h.Position=[0.5 1.03];
fig_filename = fullfile('L:\Analysis\Results\cells\figures', [cell_ID '_map_fields']);
saveas(gcf,fig_filename,'tif')
% saveas(gcf,fig_filename,'fig')



end










%%
