function wingbeat_plot_map(exp_ID)

%%
exp=exp_load_data(exp_ID);
prm = PARAMS_GetAll();

%%
figure('Units','normalized','Position',[0 0 1 1]);
directions = [-1 1];
for ii_dir = 1:2
    direction = directions(ii_dir);
    IX = ismember([exp.flight.FE.direction], direction );
    flights = exp.flight.FE(IX);
    IX = ismember(sign(exp.wingbeat.vel), direction );
    wingbeat.ts = exp.wingbeat.ts(IX);
    wingbeat.pos = exp.wingbeat.pos(IX);
    
    c = prm.graphics.colors.flight_directions{ii_dir};
    
    subplot(3,1,ii_dir+1); hold on;
    plot([flights.pos],  [flights.ts],  '.', 'Color',0.8*[1 1 1]);
    plot([wingbeat.pos], [wingbeat.ts], '.', 'Color', c)
    xlabel('Position (m)')
    ylabel('Time (min)')
    
    subplot(3,1,1); hold on;
    [f,xi] = ksdensity(wingbeat.pos, 0:0.2:200, 'Bandwidth', 0.2);
    plot(xi,f)
    xlabel('Position (m)')
    ylabel('density')
end
h=suptitle(exp.details.exp_ID)
h.Interpreter = 'none';
h.FontSize = 16;

%%
filename = fullfile('L:\Analysis\Results\exp\wingbeat',[exp_ID '_exp_wingbeat_map']);
saveas(gcf, filename,'tif')
saveas(gcf, filename,'fig')
close(gcf)
