%%
clear
clc

%%
dir_out = 'L:\Analysis\Results\cells\repeating_cells';
mkdir(dir_out);

%% load cells
cells_t = DS_get_cells_summary();
cells_t(~ismember(cells_t.bat, [79,148,34,9861,2289] ),:) = [];
prm = PARAMS_GetAll();

%% filter cells - brain region / sorting quality
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.details];
cells(~contains({cells.brain_area}, 'CA1')) = [];
cells(~ismember([cells.ClusterQuality], [2])) = [];
% cells(cellfun(@isempty, {cells.stable_ts})) = [];
cells_t = cells_t({cells.cell_ID},:);

%% filter cells - mean FR
cells = cellfun(@(c)(cell_load_data(c,'stats')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.stats];
cells = [cells.all];
cells_t([cells.meanFR_all]>prm.inclusion.interneuron_FR_thr,:) = [];

%% disp final cells table
% cells_t
whos cells_t


%% load cells FR maps
cells = cellfun(@(c)(cell_load_data(c,'FR_map')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = arrayfun(@(c)([c.FR_map(1).all.PSTH c.FR_map(2).all.PSTH]), cells, 'UniformOutput',0);
FR_maps_all = cat(1,cells{:});
% calc map correlations
FR_maps_all(isnan(FR_maps_all)) = 0;
maps_corr = corr(FR_maps_all');

%% load cells details (mask only same bat+TT)
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.details];
[X,Y] = meshgrid([cells.bat]);
same_bat_mask = (X==Y);
[X,Y] = meshgrid([cells.TT]);
same_TT_mask = (X==Y);
same_bat_TT_mask = same_TT_mask & same_bat_mask;

%% create null distribution for setting a threshold
mask = tril(~same_bat_TT_mask,-1);
maps_corr_diff_TT = maps_corr(mask);
mask = tril(same_bat_TT_mask,-1);
maps_corr_same_TT = maps_corr(mask);
figure('units','normalized','outerposition',[0 0 1 1])
hold on
clc
histogram(maps_corr_diff_TT,'Normalization','pdf')
histogram(maps_corr_same_TT,'Normalization','pdf')
p = [0.05 0.01 0.001];
p_quantiles = quantile(maps_corr_diff_TT, 1-p);
ylimits = get(gca,'ylim');
plot(repelem(p_quantiles',1,2), ylimits, '--m')
text(p_quantiles, repelem(ylimits(1)+0.95*range(ylimits),length(p_quantiles)), "" + p_quantiles)
legend({'different TT';'same TT'},'Location','eastoutside');
xlabel('Map Correlation')
ylabel('pdf')
title('cell pairs map correlations distribution')
% h=gca;
% h.YScale = 'log';
% max(maps_corr_diff_TT)
% sort(sdf,'ascend')
% ylim([0 2])
figname = fullfile(dir_out, 'map_corr_hist');
saveas(gcf,figname, 'tif')
saveas(gcf,figname, 'fig')

%% plot corr matrix
figure
corr_thr = 0.7;
imagesc( same_bat_TT_mask.*tril(maps_corr,-1) )
colorbar
set(gca,'CLim',[corr_thr 1])

%% create figures of putative same cell pair
corr_thr_plot = 0.5;
mask = tril(same_bat_TT_mask,-1);
maps_corr_same_TT = maps_corr .* mask;
[I J] = find( maps_corr_same_TT > corr_thr_plot  );
whos I J
cell_IDs = {cells.cell_ID};
pairs_IX = [I J];
pairs = cell_IDs(pairs_IX);
for ii_pair = 1:size(pairs,1)
    %% re-order pair by date
    cell1 = cell_load_data(pairs{ii_pair,1},'details');
    cell2 = cell_load_data(pairs{ii_pair,2},'details');
    [~,sort_IX] = sort([cell1.details.date cell2.details.date]','ascend');
    pairs(ii_pair,:) = pairs(ii_pair,sort_IX);
    cell_ID1 = pairs{ii_pair,1};
    cell_ID2 = pairs{ii_pair,2};
    pair_map_corr = maps_corr_same_TT(I(ii_pair), J(ii_pair));
    fprintf('%3d:\tcorr=%.2f\t%s vs. %s\n', ii_pair , pair_map_corr , cell_ID1 , cell_ID2);
%     continue
    %% create figure with panel
    figure('units','normalized','outerposition',[.1 .1 .8 .7])
    pnl = panel();
    pnl.pack(2,2,2);
%     pnl.margin = 5;
    pnl(1).de.margin = 5;
    pnl(2).de.margin = 5;
%     pnl.select('all')
%     pnl.identify()
    
    %% 
    depths = [];
    for ii_cell = 1:2
        cell_ID = pairs{ii_pair,ii_cell};
        cell = cell_load_data(cell_ID,'details','FR_map','FE');
        depths(ii_cell) = cell.details.depth;
        ti = exp_get_sessions_ti(cell.details.exp_ID, 'Behave');
        for ii_dir = 1:2
            c = prm.graphics.colors.flight_directions{ii_dir};
            % FR map
            pnl(ii_cell,ii_dir,1).select(); hold on
            plot(cell.FR_map(ii_dir).all.bin_centers, cell.FR_map(ii_dir).all.PSTH,...
                'Color', c);
            xlim([0 200])
            % trajectory + spikes
            pnl(ii_cell,ii_dir,2).select(); hold on
            plot([cell.FE{ii_dir}.pos], [cell.FE{ii_dir}.ts], '.', 'color', 0.9*[1 1 1], 'MarkerSize',0.5);
            plot([cell.FE{ii_dir}.spikes_pos], [cell.FE{ii_dir}.spikes_ts], '.', 'Color', c);
            rescale_plot_data('y',[1e-6/60 ti(1)]);
            n_flights = sum([cell.FE{ii_dir}.distance] > prm.flight.full_min_distance);
            text(0.99,0.9,num2str(n_flights), 'Units','normalized','FontWeight','bold','HorizontalAlignment','right');
            xlim([0 200])
        end
        h=pnl(ii_cell).title(cell_ID);
        h.Position = [0.5 1];
        h.Interpreter='none';
        h.FontSize = 14;
    end
    h=pnl.title(sprintf('corr=%.2f, depth_diff=%d', pair_map_corr, diff(depths)));
    h.Position = [0.95 1];
    h.Interpreter='none';
    h.FontSize = 12;
    figname = sprintf('corr_%.2f_diff_depth_%dum_%s(%d)_vs_%s(%d)', ...
        pair_map_corr, diff(depths),...
        cell_ID1, cell1.details.cell_num,...
        cell_ID2, cell2.details.cell_num);
    figname = strrep(figname, '.', '_');
    figname = fullfile(dir_out, figname);
    saveas(gcf, figname, 'tif')
    close(gcf)
end




%%









%%

