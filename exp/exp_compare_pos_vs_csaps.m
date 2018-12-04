function exp_compare_pos_vs_csaps(exp_ID)

%% load exp data
exp=exp_load_data(exp_ID,'details','position');

%% plot
figure
pnl=panel();
pnl.pack('v',2);
pnl.margin = [15 15 10 10];
pnl.de.margin = 12;
pnl(1).select()
err = exp.pos.proc_1D.pos - exp.pos.proc_1D.pos_csaps;
plot(exp.pos.proc_1D.pos_csaps, err,'.')
xlim([0 200])
ylim([-1 1])
xlabel('Position (m)')
ylabel('Pos_{raw} - Pos_{csaps} (m)')

pnl(2).select()
err = exp.pos.proc_1D.pos - exp.pos.proc_1D.pos_csaps;
plot(exp.pos.proc_1D.vel_csaps, err,'.')
xlim([-10 10])
ylim([-1 1])
xlabel('Veolocty (m/s)')
ylabel('Pos_{raw} - Pos_{csaps} (m)')

h=pnl.title([exp_ID ' Compare raw vs. csaps position']);
h.Interpreter='none';h.Position=[0.5 1.01];h.FontSize=14;

%% save figure
fig_filename = fullfile('L:\Analysis\Results\exp\position', [exp_ID '_exp_position_raw_vs_csaps']);
saveas(gcf,fig_filename,'tif')
% saveas(gcf,fig_filename,'fig')

end