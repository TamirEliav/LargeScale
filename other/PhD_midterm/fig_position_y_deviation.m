%% create all exp position data
clear
clc

%% load cells summary and choose cells
cells_t = DS_get_cells_summary();
cells_t(~ismember(cells_t.bat, [79,148,34,9861] ),:) = [];
cells_t(~strcmp(cells_t.brain_area, 'dCA1'),:)=[];

%% only signif cells
cells = cellfun(@(x)(cell_load_data(x,'signif')),cells_t.cell_ID);
signif = cat(1,cells.signif);
signif = arrayfun(@(x)(x.TF), signif);
signif = any(signif,2);
cells_t = cells_t(signif,:);

%% get the relevant exp list
cells = cellfun(@(x)(cell_load_data(x,'details')),cells_t.cell_ID);
exp_list = unique(cellfun(@(x)(x.exp_ID),{cells.details},'uniformoutput',false));

%%
exps = cellfun(@(x)(exp_load_data(x,'pos_y_std')),exp_list);
ystd_median_all = arrayfun(@(x)(x.ystd_median), cat(1,exps.pos_y_std));

%% create figure!
prm = PARAMS_GetAll();
dir_colors = prm.graphics.colors.flight_directions;
example_exp_ID = 'b9861_d180709';

figure('Units','centimeters','Position',[5 5 16 4])
pnl=panel();
pnl.pack('h',[80 20])

% example
pnl(1).select();hold on
exp = exp_load_data(example_exp_ID,'pos_y_std');
for ii_dir = 1:2
    % plot the actual data points
    plot(exp.pos_y_std(ii_dir).xy(:,1),...
         exp.pos_y_std(ii_dir).xy(:,2),...
         '.', 'color', dir_colors{ii_dir},'MarkerSize',1e-20);
    
    % plot the mean+std
    h=shadedErrorBar(exp.pos_y_std(ii_dir).bin_centers,...
                     exp.pos_y_std(ii_dir).ymean,...
                     exp.pos_y_std(ii_dir).ystd,...
                     'lineprops',{'color', dir_colors{ii_dir}});
    h.patch.FaceAlpha=0.05;
	ylim([-2 2])
end
ylabel('Y Position (m)')


% population
pnl(2).select();hold on
h1=histogram(ystd_median_all(:,1),'FaceColor',dir_colors{1},'Normalization','probability','FaceAlpha',0.5);
h2=histogram(ystd_median_all(:,2),'FaceColor',dir_colors{2},'Normalization','probability','FaceAlpha',0.5);
xlabel({'Y position';...
        'deviation (m)'})
ylabel('Prob.');


%% save figure
dir_out = 'L:\Analysis\Results\midterm';
mkdir(dir_out);
filename = fullfile(dir_out, 'fig_pos_y_deviation')
set(gcf, 'Renderer', 'painters');
% set(gcf, 'Renderer', 'zbuffer');
% set(gcf, 'Renderer', 'opengl');
saveas(gcf,filename,'fig')
saveas(gcf,filename,'tif')
saveas(gcf,filename,'pdf')

%%









%%

