function POS_plot_flight(exp_ID)

%% load data
exp = exp_load_data(exp_ID, 'details', 'pos', 'flight');
prm = PARAMS_GetAll();

%% arrange data
FE = exp.flight.FE;
directions = [1 -1];
session_ti = exp_get_sessions_ti(exp_ID, 'Behave');
vel_limits = [-1 1].*12;
pos_limits = [0 200];
dir_colors = prm.graphics.colors.flight_directions;

%%
figure('Units','normalized','Position',[0 0 1 1])
pnl = panel();
pnl.pack('h',[60 40])
pnl(1).pack('h',2)
pnl(2).pack('v',3)
pnl(2,3).pack('h',3)
pnl.margin = [20 30 20 15];
pnl.de.margin = 15;
% h=pnl.title( sprintf('%s             bsp tag ID: %d, fs=%0.2f', exp_ID, exp.details.bsp_tag_ID, exp.pos.raw.fs) );
h=pnl.title(exp_ID);
h.Interpreter = 'none'; h.FontSize=16; h.Position = [0.5 1.03];

%% plot - position by flight number
pnl(1,1).select();
hold on
for ii_dir = 1:2
    IX = find([FE.direction] == directions(ii_dir));
    x = [FE(IX).pos];
    y = repelem([FE(IX).number],cellfun(@length,{FE(IX).ts}));
    plot(x, y,'.', 'Color', dir_colors{ii_dir});
end
xlim(pos_limits)
xlabel('X (m)')
ylabel('flight#')
% h=title('raw data');
% h.Units = 'normalized'; h.Position = [0.05 1]; h.HorizontalAlignment = 'Left'; h.FontSize = 12;

%% plot - position by time
pnl(1,2).select();
hold on
for ii_dir = 1:2
    IX = find([FE.direction] == directions(ii_dir));
    x = [FE(IX).pos];
    y = [FE(IX).ts];
    plot(x, y,'.', 'Color', dir_colors{ii_dir});
end
rescale_plot_data('y',[1e-6/60 session_ti(1)])
xlim(pos_limits)
xlabel('X (m)')
ylabel('Time (min)')

% h=title('raw data');
% h.Units = 'normalized'; h.Position = [0.05 1]; h.HorizontalAlignment = 'Left'; h.FontSize = 12;

%% plot - velocity by time
pnl(2,1).select();
hold on
for ii_dir = 1:2
    IX = find([FE.direction] == directions(ii_dir));
    x = [FE(IX).ts];
    y = [FE(IX).vel];
    plot(x, y,'.', 'Color', dir_colors{ii_dir});
end
rescale_plot_data('x',[1e-6/60 session_ti(1)])
ylim(vel_limits)
xlabel('Time (min)')
ylabel('Velocity (m/s)')

%% plot - velocity vs. position
pnl(2,2).select();
hold on
for ii_dir = 1:2
    IX = find([FE.direction] == directions(ii_dir));
    x = [FE(IX).pos];
    y = [FE(IX).vel];
    plot(x, y,'.', 'Color', dir_colors{ii_dir});
end
ylim(vel_limits)
xlim(pos_limits)
xlabel('Position (m)')
ylabel('Velocity (m/s)')

%% Distance histograms
pnl(2,3,1).select();
cla
hold on
for ii_dir = 1:2
    IX = find([FE.direction] == directions(ii_dir));
    h=histogram([FE(IX).distance]);
    h.BinEdges = linspace(0,200,14);
    h.FaceColor = dir_colors{ii_dir};
end
ha = gca;
ha.XTick = 0:50:200;
xlabel('Flight distance (m)')
ylabel('Counts')

%% Duration histograms
pnl(2,3,2).select();
cla
hold on
for ii_dir = 1:2
    IX = find([FE.direction] == directions(ii_dir));
    durations = [FE(IX).duration];
    h=histogram(durations);
    m = 35;
%     m = max(m,max(durations));
    h.BinEdges = linspace(0,m,14);
    h.FaceColor = dir_colors{ii_dir};
end
ha = gca;
ha.XTick = 0:10:30;
xlabel('Flight duration (s)')
ylabel('Counts')


%% add text details
pnl(2,3,3).select();
h=text(0.1,0.9,{...
    sprintf('speed low thr = %d m/s', exp.flight.speed_low_thr),...
    sprintf('speed high thr = %d m/s', exp.flight.speed_high_thr),...
    sprintf('high speed min duration = %d s', exp.flight.high_speed_min_duration),...
    });
h.Units = 'normalized';
axis off

%% save figure
fig_filename = fullfile('L:\Analysis\Results\exp\flight', [exp_ID '_exp_flight']);
saveas(gcf,fig_filename,'tif')
% saveas(gcf,fig_filename,'fig')




end




