function exp_plot_pos_std_y(exp_ID)

%% load exp data
exp = exp_load_data(exp_ID,'details','pos_y_std');
prm = PARAMS_GetAll();
dir_colors = prm.graphics.colors.flight_directions;

%%
figure('Units','normalized','Position',[0.1 0.1 0.8 0.8])
pnl = panel();
pnl.pack('h',[80 20],'v',2);
pnl.margin = [15 25 10 15];
pnl.de.margin = 10;
h=pnl.title(exp_ID); h.FontSize=16; h.Interpreter='none';h.Position=[0.5 1.02];

%% TZAMOT
pnl(1,1).select(); hold on
for ii_dir = 1:2
    % plot the actual data points
    plot(exp.pos_y_std(ii_dir).xy(:,1),...
         exp.pos_y_std(ii_dir).xy(:,2),...
         '.', 'color', dir_colors{ii_dir},'MarkerSize',1);
    
    % plot the mean+std
    h=shadedErrorBar(exp.pos_y_std(ii_dir).bin_centers,...
                     exp.pos_y_std(ii_dir).ymean,...
                     exp.pos_y_std(ii_dir).ystd,...
                     'lineprops',{'color', dir_colors{ii_dir}});
	ylim([-2 2])
end
ylabel('Y Position (m)')

%% y std + y range
pnl(1,2).select(); hold on
for ii_dir = 1:2
    plot(exp.pos_y_std(ii_dir).bin_centers,...
         exp.pos_y_std(ii_dir).ystd,...
         '-','color', dir_colors{ii_dir}, 'LineWidth',2);
     
     plot(exp.pos_y_std(ii_dir).bin_centers,...
          diff(exp.pos_y_std(ii_dir).yrange),...
         '--','color', dir_colors{ii_dir}, 'LineWidth',2);
end
columnlegend(2,{...
    sprintf('std %d%% trim',exp.pos_y_std(ii_dir).trim_prc),...
    sprintf('range %d-%d%%',[exp.pos_y_std(ii_dir).range_prc(:)]),...
    sprintf('std %d%% trim',exp.pos_y_std(ii_dir).trim_prc),...
    sprintf('range %d-%d%%',[exp.pos_y_std(ii_dir).range_prc(:)]),...
     } );
xlabel('X Position (m)')
ylabel('Y-deviation (m)')

%% y std + y range
pnl(2,2).select(); hold on
for ii_dir = 1:2
    dev_std = exp.pos_y_std(ii_dir).ystd;
    dev_range = diff(exp.pos_y_std(ii_dir).yrange);
    [N1,edges1] = histcounts(dev_std,  'Normalization','probability');
    [N2,edges2] = histcounts(dev_range, 'Normalization','probability');
    plot(N1, edges2centers(edges1), '-',  'Color', dir_colors{ii_dir}, 'LineWidth',2);
    plot(N2, edges2centers(edges2), '--', 'Color', dir_colors{ii_dir}, 'LineWidth',2);
    h=refline(0, nanmedian(dev_std));   h.LineStyle = '-';  h.Color = dir_colors{ii_dir};
    h=refline(0, nanmedian(dev_range)); h.LineStyle = '--'; h.Color = dir_colors{ii_dir};
end
xlabel('Prob.')

%% linkaxes
linkaxes([pnl(2,2).axis pnl(1,2).axis],'y')

%% save figure
fig_filename = fullfile('L:\Analysis\Results\exp\pos_y_std', [exp_ID '_pos_y_std']);
saveas(gcf,fig_filename,'tif')
% saveas(gcf,fig_filename,'fig')



end



