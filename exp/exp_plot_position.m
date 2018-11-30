function exp_plot_position(exp_ID)

%% load exp data
exp = exp_load_data(exp_ID,'details','position');
prm = PARAMS_GetAll();

%%
figure('Units','normalized','Position',[0 0 1 1])
pnl = panel();
pnl.pack('v',[20 20 30 30])
pnl.margin = [20 30 20 15];
pnl.de.margin = 15;
h=pnl.title( sprintf('%s             bsp tag ID: %d, fs=%0.2f', exp_ID, exp.details.bsp_tag_ID, exp.pos.raw.fs) );
% h=pnl.title(exp_ID);
h.Interpreter = 'none'; h.FontSize=16; h.Position = [0.5 1.03];

%% raw
pnl(1).select();
plot(exp.pos.raw.pos(:,1), exp.pos.raw.pos(:,2), '.k')
xlim([1320 1500])
ylim([2460 2505]);
h=gca;
h.XDir = 'reverse';
xlabel('X (m)')
ylabel('Y (m)')
h=title('raw data');
h.Units = 'normalized'; h.Position = [0.05 1]; h.HorizontalAlignment = 'Left'; h.FontSize = 12;

%% projected
pnl(2).select();
hold on
directions = sign(diff(exp.pos.proj.pos(:,1)));
dir_colors = prm.graphics.colors.flight_directions;
dir_colors{3} = [0 0 0];
dir_sign = [1 -1 0];
for ii_dir = 1:3
    IX = directions == dir_sign(ii_dir);
    plot(exp.pos.proj.pos(IX,1), exp.pos.proj.pos(IX,2), '.', 'Color', dir_colors{ii_dir})
end
xlim([0 200])
ylim([-1 1].*2)
xlabel('X (m)')
ylabel('Y (m)')
h=legend({'dir 1: -->';'dir 2: <--';'static'});
h.Location = 'northeastoutside';
h.Position(3:4) = [0.05 0.05];
h=title('project raw to tunnel midline');
h.Units = 'normalized'; h.Position = [0.05 1]; h.HorizontalAlignment = 'Left'; h.FontSize = 12;

%% processed - velocity vs. position
pnl(3).select();
hold on
plot(exp.pos.proc_1D.pos, exp.pos.proc_1D.vel,'-k')
plot(exp.pos.proc_1D.pos_csaps, exp.pos.proc_1D.vel_csaps,'-r')
xlim([0 200])
ylim([-1 1].*12)
xlabel('Position (m)')
ylabel('velocity (m/s)')
h=legend({'raw';'csaps';});
h.Location = 'northeastoutside';
h.Position(3:4) = [0.05 0.05];
h=title('processed (interp/extrap/upsample) - velocity vs. position');
h.Units = 'normalized'; h.Position = [0.05 1]; h.HorizontalAlignment = 'Left'; h.FontSize = 12;

%% processed - position vs. time
pnl(4).select();
hold on
ib = find_nearest_point(exp.pos.proc_1D.ts, exp.pos.raw.ts_nlg_usec);
dist_from_orig_ts = abs(exp.pos.proc_1D.ts - exp.pos.raw.ts_nlg_usec(ib)');
IX = dist_from_orig_ts > (median(diff(exp.pos.raw.ts_nlg_usec)));
plot(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos,'.k')
plot(exp.pos.proc_1D.ts(IX), exp.pos.proc_1D.pos(IX),'.r')
rescale_plot_data('x',[1e-6/60 exp.pos.proc_1D.ts(1)])
xlabel('Time (min)')
ylabel('Position (m)')
h=legend({'original';'interp/extrap';});
h.Location = 'northeastoutside';
h.Position(3:4) = [0.05 0.05];
h=title('processed (interp/extrap/upsample) - position vs. time');
h.Units = 'normalized'; h.Position = [0.05 1]; h.HorizontalAlignment = 'Left'; h.FontSize = 12;

%% save figure
fig_filename = fullfile('L:\Analysis\Results\exp\position', [exp_ID '_exp_position'])
saveas(gcf,fig_filename,'tif')
saveas(gcf,fig_filename,'fig')



end




