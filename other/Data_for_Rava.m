
%% data for Rava
clear
clc
load('L:\processed_data_structs\cells_bat_200m.mat')

%%
data = struct();

%%
signif = cat(1,cells.signif);
signif = arrayfun(@(x)(x.TF),signif);
maps = cat(1,cells.FR_map);
maps = maps(signif);
maps = [maps.all];
fields_per_dir = cat(1,cells(:).fields);
fields_per_dir = fields_per_dir(signif);

%%
data.bin_centers = maps(1).bin_centers;
data.maps = cat(1,maps.PSTH);
data.maps01 = zeros(size(data.maps));
for ii_map = 1:length(maps)
    fields = fields_per_dir{ii_map};
    fields([fields.in_low_speed_area]) = [];
    fields_per_dir{ii_map} = fields;
    fields_edges = cat(1,fields.edges_prc);
    IX = get_data_in_ti(data.bin_centers,fields_edges);
    data.maps01(ii_map,IX) = 1;
end

%% order by number of fields
nfields = cellfun(@length, fields_per_dir);
[~,sorted_IX] = sort(nfields);
data.maps = data.maps(sorted_IX,:);
data.maps01 = data.maps01(sorted_IX,:);
% set nans also in the binary maps
data.maps01(isnan(data.maps)) = nan;

%% save the data
save('L:\DATA_for_people\Simon_Rava\data','data')












%% load LM data
exp_ID = 'b2289_d180615';
exp = exp_load_data(exp_ID,'LM');
LM = exp.LM;
% LM( contains({LM.name},{'ball','enter'}) ) = [];

%% plot maps
figure
subplot(311)
imagesc(data.bin_centers,[],data.maps)
subplot(312)
imagesc(data.bin_centers,[],data.maps01)
subplot(313)
yyaxis left
plot(data.bin_centers,sum(data.maps))
ylabel('PSTH sum')
yyaxis right
plot(data.bin_centers,sum(data.maps01))
ylabel('fields sum')

%%
ccc=corr(data.maps);
ccc01=corr(data.maps01);

%%
figure
subplot(221)
hold on
imagesc(data.bin_centers,data.bin_centers,ccc);
% plot_LM(LM)
axis square
subplot(222)
imagesc(data.bin_centers,data.bin_centers,ccc01)
% plot_LM(LM)
axis square
subplot(223)
yyaxis left
plot(data.bin_centers,mean(data.maps));
ylabel('mean FR')
yyaxis right
plot(data.bin_centers,nanmean(ccc));
ylabel('mean corr')
subplot(224)
yyaxis left
plot(data.bin_centers,mean(data.maps01));
ylabel('mean FR')
yyaxis right
plot(data.bin_centers,nanmean(ccc01));
ylabel('mean corr')

%%
figure
hs=scatter(nanmean(ccc01), mean(data.maps01), '.');
hs.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('position',data.bin_centers);
xlabel('mean corr')
ylabel('mean FR')


%%
figure
hs=scatter3(data.bin_centers, nanmean(ccc01), mean(data.maps01), 'o');
% hs.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('position',data.bin_centers);
xlabel('position')
ylabel('mean corr')
zlabel('mean FR')


%%






%%
function plot_LM(LM)
h=arrayfun(@xline, [LM.pos_proj]);
[h.Color] = disperse(repelem('r',length(LM)));
[h.LineWidth]= disperse(repelem(0.1,length(LM)))
h=arrayfun(@yline, [LM.pos_proj]);
[h.Color] = disperse(repelem('r',length(LM)));
[h.LineWidth]= disperse(repelem(0.1,length(LM)))
% text()
end


