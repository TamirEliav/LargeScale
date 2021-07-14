function MUA_plot(exp_ID)

%% get exp info
exp = exp_load_data(exp_ID,'details','MUA');

%% get data
trig = exp.MUA.trig;

%% plot MUA triggered FR
figure
subplot(211)
shadedErrorBar(trig.t, trig.FR_mean,trig.FR_sem)
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
subplot(212)
shadedErrorBar(trig.t, trig.zFR_mean,trig.zFR_sem)
xlabel('Time (ms)')
ylabel('Firing rate (z)')

sgtitle({exp.details.exp_ID;'MUA event triggered Firing Rate'},'Interpreter','None');

fig_filename = fullfile('L:\Analysis\Results\exp\MUA', [exp_ID '_MUA_events']);
saveas(gcf,fig_filename,'jpeg')



end