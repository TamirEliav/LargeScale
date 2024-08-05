%%
clear
clc
close all

%%
load('L:\paper_replay\figures\single_unit_replay_tuning_corr.mat')

%%
ccc_all = [];
cell_IDs_all = [];
cell_num_all = [];
map_dir_all = [];
for ii_dir = 1:2
    ccc_dir = ccc_data{ii_dir};
    details_dir = details(cells_inc_TF_per_dir{ii_dir});
    cell_IDs_dir = {details_dir.cell_ID}';
    cell_num_dir = [details_dir.cell_num]';
    map_dir_dir = ones(size(ccc_dir)).*ii_dir;
    ccc_all = [ccc_all; ccc_dir];
    cell_IDs_all = [cell_IDs_all; cell_IDs_dir];
    cell_num_all = [cell_num_all; cell_num_dir];
    map_dir_all = [map_dir_all; map_dir_dir];
end

%% sort by corr 
[~,sort_IX] = sort(ccc_all,'descend');
ccc_all = ccc_all(sort_IX);
cell_IDs_all = cell_IDs_all(sort_IX);
cell_num_all = cell_num_all(sort_IX);
map_dir_all = map_dir_all(sort_IX);

%% add more info
num_fields_all = zeros(size(ccc_all));
SI_all = zeros(size(ccc_all));
for ii_cell = 1:length(cell_IDs_all)
    cell_ID = cell_IDs_all{ii_cell};
    ii_dir = map_dir_all(ii_cell);
    cell = cell_load_data(cell_ID,'stats');
    num_fields_all(ii_cell) = cell.stats.dir(ii_dir).field_num_include_near_balls;
    SI_all(ii_cell) = cell.stats.dir(ii_dir).field_num_include_near_balls;
end

%% choose cells
TF = true(size(ccc_all));
TF = TF & ccc_all>0.5;
TF = TF & num_fields_all>1;
sum(TF)
cell_IX_list = find(TF);

%%
% cell_IX_list = find(ismember(cell_num_all,...
%     [766 1444 1379 1020 1329 1059 1017 830 1490 1540 1016 72 746 1022 178 1383 853 1395 1560 1309 1552 1331 1064 1540 211 216 1554 499 ]));
whos cell_IX_list


%% plot
ncol = 7;
nrow = 12;
nExamples2plot = nrow*ncol;
ii_cat = 3;
fig = figure();
fig.WindowState = 'maximized';
h=tiledlayout(nrow,ncol,'TileSpacing','tight');
h.Padding = 'compact';
for ii_ex = 1:length(cell_IX_list)
    ii_cell = cell_IX_list(ii_ex);
    nexttile
    hold on
    cell_ID = cell_IDs_all{ii_cell};
    cell_num = cell_num_all(ii_cell);
    cell = cell_load_data(cell_ID,'FR_map','replay_FR_map');
    ii_dir = map_dir_all(ii_cell);
    x = cell.FR_map(ii_dir).all.bin_centers;
    y = cell.FR_map(ii_dir).all.PSTH;
    y = y./max(y);
    plot(x,y,'k-','LineWidth',3);
    x = cell.replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).bin_centers;
    y = cell.replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).PSTH;
    y = y./max(y);
    plot(x,y,'r-','LineWidth',1.5);
    title_str = {cell_ID;sprintf('cell %d, dir %d, r=%.2f', cell_num, map_dir_all(ii_cell), ccc_all(ii_cell))};
    text(.5,1.15, title_str ,'HorizontalAlignment','center','units','normalized', 'FontSize',6,'Interpreter','None');
    axis off
end
res_dir =  'L:\paper_replay\figures';
fig_name = 'fig_2_replay_tuning_corr_examples.pdf';
file_out = fullfile(res_dir, fig_name);
saveas(fig, file_out, 'pdf')
exportgraphics(fig,file_out)

%%








%%