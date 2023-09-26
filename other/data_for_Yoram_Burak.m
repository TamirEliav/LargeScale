%%
clear
clc

%%
load('L:\processed_data_structs\cells_bat_200m.mat');

%%
signif = arrayfun(@(x)(x.TF), cat(1,cells.signif));
FR_maps_all = cat(1,cells.FR_map);
FR_maps_all = reshape([FR_maps_all.all],size(FR_maps_all,1),size(FR_maps_all,2),[]);
M = cat(1,FR_maps_all.PSTH);
M = reshape(M,size(FR_maps_all,1),size(FR_maps_all,2),[]);
signif = repmat(signif,1,1,size(M,3));
M(~signif) = nan;

%% arrange data
signif = cat(1,cells.signif);
TF = arrayfun(@(x)(x.TF),signif);
maps = cat(1,cells.FR_map);
PSTHs = [];
for ii_cell=1:size(maps,1)
    for ii_dir = 1:2
        PSTHs(ii_cell, ii_dir,:) = maps(ii_cell, ii_dir).all.PSTH;
        if ~TF(ii_cell, ii_dir)
            PSTHs(ii_cell, ii_dir,:) = nan;
        end
    end
end
whos PSTHs

%%
