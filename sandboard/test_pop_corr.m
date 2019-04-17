%% cell analysis pipline
clear
clc

%% load cells summary and choose cells
bats = [9861];
clstr_Q_filt = [1];
% map_type = 'FR'; % FR / peak-norm FR / mean-norm FR / Ipos
map_type = 'Ipos'; % FR / peak-norm FR / mean-norm FR / Ipos
map_normalization = 'no'; % no / range
corr_type = 'spearman'; % pearson / spearman / ?
% corr_type = 'pearson'; % pearson / spearman / ?

cells_t = DS_get_cells_summary();
cells_t(~ismember(cells_t.bat, bats),:) = [];
cells_t(~ismember(cells_t.ClusterQuality, clstr_Q_filt),:) = [];
Ncells = height(cells_t);
fprintf('num cells: %d\n', Ncells );

%% load data for all cells
cell = cell_load_data(cells_t.cell_ID{1},'FR_map');
% maps = nan(Ncells , 2*length(cell.FR_map(1).all.PSTH));
maps = [];
for ii_cell = 1:Ncells
    cell = cell_load_data(cells_t.cell_ID{ii_cell}, 'FR_map','Ipos');
    switch map_type 
        case 'FR'
            maps(ii_cell, :) = [cell.FR_map(1).all.PSTH cell.FR_map(2).all.PSTH];
        case 'Ipos'
            maps(ii_cell, :) = [cell.Ipos.data(1).Ipos cell.Ipos.data(2).Ipos];
    end
end
switch map_normalization 
    case 'no'
    case 'range'
        maps = normalize(maps,2,'range');
end

%%
switch corr_type
    case 'spearman'
        maps_corr = corr(maps, 'Type', 'Spearman');
    case 'pearson'
        maps_corr = corr(maps, 'Type', 'Pearson');
end

%%
figure
mesh(maps)
colorbar
% set(gca, 'CLim', [0 20])

%%
figure
imagesc(maps_corr)
axis equal
axis xy
colorbar
colormap jet



%%
