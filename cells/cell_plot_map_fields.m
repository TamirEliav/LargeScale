function cell_plot_map_fields(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','FE','FR_map','Ipos','fields','RecStability');
exp = exp_load_data(cell.details.exp_ID,'details');
prm = PARAMS_GetAll();
dir_colors = prm.graphics.colors.flight_directions;

%% create figure
figure('Units','normalized','Position',[0 0 1 1]);
pnl = panel();
pnl.pack('h',[30 70])
pnl(1).pack('v',2)
pnl(1,1).pack('v',2)
pnl(1,2).pack('v',2)
pnl(2).pack('v',2)
pnl.margin = [15 25 10 15];
pnl(1).margin = 100;
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
for ii_dir = 1:2
    fields = cell.fields{ii_dir};
    pnl(2,ii_dir).pack('h',length(fields)+1);
    
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
    
    %% plot each field spikes seperaetly
    for ii_field = 1:length(fields)
        field = fields(ii_field);
        pnl(2,ii_dir,ii_field+1).select(); hold on;
        field_spikes_wvfrms = [field.spikes_wvfrm];
        plot(mean(field_spikes_wvfrms,3))
    end
end
linkaxes(pnl(2).de.axis, 'xy');

%% save figure
h=pnl.title(cell_ID); h.FontSize=16; h.Interpreter='none';h.Position=[0.5 1.03];
fig_filename = fullfile('L:\Analysis\Results\cells\figures', [cell_ID '_map_fields']);
saveas(gcf,fig_filename,'tif')
% saveas(gcf,fig_filename,'fig')



end










%%
