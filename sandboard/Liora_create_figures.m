%%
clear
clc

%%
cell_ID = 'b0034_d180312_TT4_SS01';
% cell_ID = 'b0079_d160928_TT2_SS01';

%% load cell/exp data
cell = cell_load_data(cell_ID);
exp = exp_load_data(cell.details.exp_ID,'details','LM');
prm = PARAMS_GetAll();
dir_colors = prm.graphics.colors.flight_directions;

%% create figure
figure('Units','centimeters','Position',[5 5 15 15]);
pnl = panel();
pnl.pack('v',[30 60])
pnl(2).pack('v',2)
pnl.margin = [15 5 5 5];
pnl(1).margin = 10;
pnl(1).de.margin = 10;

%% FR map + fields + LM
pnl(1).select();hold on
for ii_dir = 1:2
    % FR map
    plot(cell.FR_map(ii_dir).all.bin_centers, cell.FR_map(ii_dir).all.PSTH,...
        'Color', dir_colors{ii_dir}, 'LineWidth',1.5);
    % labels
%     xlabel('Position (m)')
    ylabel('Firing rate (Hz)')
end

%% trajectory + spikes
ti = exp_get_sessions_ti(cell.details.exp_ID, 'Behave');
for ii_dir = 1:2
    pnl(2,ii_dir).select();hold on
    FE = cell.FE{ii_dir};
%     plot([FE.pos],       [FE.ts],        '.', 'Color', 0.9.*[1 1 1])
    plot([FE.spikes_pos],[FE.spikes_ts], '.', 'Color', dir_colors{ii_dir})
    xlimits = get(gca, 'xlim');
    if ~isempty(cell.details.stable_ts)
        plot(repmat(xlimits,[2 1])', repmat(cell.details.stable_ts, [ 2 1] ),'--', 'Color','k');
    end
    ylim(ti)
    rescale_plot_data('y',[1e-6/60 ti(1)])
    % labels
    ylabel('Time (min)')
    if ii_dir==1
        set(gca,'xtick',[])
    else
        xlabel('Position (m)')
    end
    set(gca,'tickdir','out')
end

%% add arrows for directions
arrow_pos_X = [.01 .10;
               .10 .01]+0.84;
arrow_pos_Y = [repelem(0.825, 2);
               repelem(0.805, 2)]+0.15;
for ii_dir = 1:2
    c = dir_colors{ii_dir};
    ah=annotation('arrow',arrow_pos_X(ii_dir,:),arrow_pos_Y(ii_dir,:),'Color',c);
    ah.LineWidth = 1.5;
    set(ah,'HeadStyle','cback1','HeadWidth',5);
end

%% save figure
fig_filename = fullfile('L:\Analysis\Results\posters_presentations\Liora_20190523', [cell_ID '_map_fields']);
saveas(gcf,fig_filename,'tif')
saveas(gcf,fig_filename,'pdf')





