%%
clear
clc

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
clear cells
% cells_t
whos cells_t

%% load cells details
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells_details = [cells.details];
cells_exp_ID = {cells_details.exp_ID}';
cells_ID = {cells_details.cell_ID}';

%% calc num cells per exp
[C,IA,IC] = unique(cells_exp_ID);
cells_per_exp = accumarray(IC,ones(size(IC)));
[cells_per_exp,sort_IX] = sort(cells_per_exp,'descend');
exp_list = C(sort_IX);
figure
histogram(cells_per_exp)
% histogram(categorical(cells_per_exp))

%% 
prm = PARAMS_GetAll();
for ii_exp = 1:3
    figure('Units','normalized','Position',[0 0 1 1])
    pnl = panel();
    pnl.pack(3,2);
    exp_ID = exp_list{ii_exp};
    sim_cells = cells_ID(ismember(cells_exp_ID, exp_ID));
    for ii_cell = 1:length(sim_cells)
        cell_ID = sim_cells{ii_cell};
        cell = cell_load_data(cell_ID, 'details', 'FR_map','Ipos');
        for ii_dir=1:2
            pnl(1,ii_dir).select(); hold on
            x = cell.FR_map(ii_dir).all.bin_centers;
            y = cell.FR_map(ii_dir).all.PSTH;
            plot(x,y);
            
            pnl(2,ii_dir).select(); hold on
            x = cell.FR_map(ii_dir).all.bin_centers;
            y = cell.FR_map(ii_dir).all.PSTH;
            y = y./max(y);
            plot(x,y);
            
            pnl(3,ii_dir).select(); hold on
            x = cell.Ipos.data(ii_dir).pos_bins_centers;
            y = cell.Ipos.data(ii_dir).Ipos;
            plot(x,y)
        end
    end
end







%%












%%
