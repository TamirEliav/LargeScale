%%
clear 
clc

%% run over cells from Shir analysis 
dir_IN = 'L:\Analysis\from_Shir_20181011\Wild_cells_rearranged';
files = dir([dir_IN '\cell_data*b*d*TT*SS*'])

%%
FR_maps_all = [];
tic
for ii_file = 1:length(files)
    
    %% load cell data
    disp( sprintf('loading cell %d/%d',ii_file,length(files)) )
    filename = fullfile(files(ii_file).folder, files(ii_file).name);
    cell_data = load(filename, 'data');
    
    %% re arrange data struct
%     cell_data.data = [cell_data.data.direction1 cell_data.data.direction2];
%     cell_data.cell_ID = sprintf('b%s_d%s_T%s_SS%s',...
%                                 cell_data.bat,...
%                                 cell_data.day,...
%                                 cell_data.TT,...
%                                 cell_data.SS );
                            
    %% do some calculations...
%     FR_calc_local_SI_map(cell_data)
%     close all

    %% 
    FR_maps_all(ii_file,:,:) = cat(1,cell_data.data.PSTH);
    
end
toc

%% plot Firing rate maps
figure
pnl=panel();
pnl.pack(2,2)
% p.select('all')
% p.identify()
pnl(1,1).select()
plot(squeeze(FR_maps_all(:,1,:))')
pnl(2,1).select()
plot(squeeze(nanmean(FR_maps_all(:,1,:)) ))
pnl(1,2).select()
plot(squeeze(FR_maps_all(:,2,:))')
pnl(2,2).select()
plot(squeeze(nanmean(FR_maps_all(:,2,:)) ))


%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% run over cells from Shir analysis 
% bat = '0148';
% bat = '0034';
% bat = '9861';
bat = '0079';
dir_IN = 'L:\Analysis\Results\Ipos';
files = dir([dir_IN '\*b' bat '*Ipos.mat'])

%%
Ipos_maps_all = [];
% Ipos_maps_all = nan( length(Ipos_maps_all), 2, 200);
for ii_file = 1:length(files)
    
    %% load cell data
    ii_file
    filename = fullfile(files(ii_file).folder, files(ii_file).name);
    clear Ipos
    load(filename);
    
    %% 
    for ii_dir = 1:2
        Ipos_maps_all(ii_file,ii_dir, :) = Ipos.data(ii_dir).Ipos;
    end
    
    
end

%% plot Ipos maps
figure
pnl=panel();
set(pnl.figure, 'units', 'centimeters', 'Position', [2 2 35 22])
pnl.pack(6,2)
pnl.margintop = 10;
pnl.de.margin = [5 10 5 5];
h=pnl.title(['bat ' bat]); set(h,'Fontsize', 16);
% pnl.select('all')
% pnl.identify()
for ii_dir = 1:2
    maps = squeeze(Ipos_maps_all(:,ii_dir,:));
    sdf_diff = diff(maps,1,2);
    sdf_diff = [sdf_diff nan(size(maps,1),1)];
%     sdf_diff = [nan(size(sdf,1),1) sdf_diff ];
    sdf_diff_smooth = arrayfun(@smooth, sdf_diff);
    x = Ipos.data.pos_bins_centers;
    ncells = size(maps,1);
    
    pnl(1,ii_dir).select()
    plot(x, maps')
    if ii_dir == 1
        ylabel('Ipos')
    end
    
    pnl(2,ii_dir).select(); hold on
    plot(x, nanmean(maps))
    if ii_dir == 1
        ylabel('Avg. Ipos')
    end
    
    pnl(3,ii_dir).select(); hold on
    imagesc(x,1:ncells, log2(maps))
    axis tight
%     colorbar
%     colormap jet
%     set(gca,'CLim',[0.2 max(maps(:))])
    if ii_dir == 1
        ylabel('Ipos')
    end
    
    pnl(4,ii_dir).select(); hold on
    plot(x, mean(abs(sdf_diff).*maps))
%     plot(x, abs(sdf_diff).*maps)
    draw_LM(LM)
    if ii_dir == 1
        ylabel('Avg abs diff Ipos')
    end
    
    pnl(5,ii_dir).select(); hold on
    plot(x, sdf_diff_smooth')
    
    pnl(6,ii_dir).select(); hold on
    plot(x, mean(abs(sdf_diff_smooth)))
    xlabel('Position (m)')
end
linkaxes(pnl.de.axis,'x')
xlim([0 200])
pnl(1,1).select()
text(0.5,1.1,'>>>>>', 'Units', 'normalized')
pnl(1,2).select()
text(0.5,1.1,'<<<<<', 'Units', 'normalized')

filename = ['Ipos_pop_b' bat];
dir_out = dir_IN;
file_out = fullfile(dir_out, filename);
saveas(gcf, file_out, 'tif')
saveas(gcf, file_out, 'fig')

%%
pop_vec = [squeeze(Ipos_maps_all(:,1,:)) squeeze(Ipos_maps_all(:,2,:))];
figure
imagesc(corr(pop_vec))

%%
function draw_LM(LM)
%     IX = strcmp( LM.name,'polygal');
%     IX = 1:size(LM,1);

    y = get(gca,'ylim');

    % special mark important LM
    IX = [4 8 11 13 14];
    x = LM.pos_linearized(IX,:);
    
    strs = LM.name(IX,:);
    plot(repmat(x,1,2), y, 'm', 'LineWidth', 2);
    text( x, repmat(y(2),1,length(x)), strs , 'HorizontalAlignment', 'center', 'Rotation', 30, 'FontSize', 8, 'verticalalignment', 'bottom');
    
    % less prominent LM
    IX = setdiff(1:size(LM,1),IX);
    x = LM.pos_linearized(IX,:);
    y = get(gca,'ylim');
    strs = LM.name(IX,:);
    plot(repmat(x,1,2), y, 'c');
    text( x, repmat(y(2)+0.05*diff(y),1,length(x)), strs , 'HorizontalAlignment', 'center', 'Rotation', 30, 'FontSize', 6);
end















%%
