%%
clear
clc
close all
exp_ID = 'b0184_d191127';
exp_create_details(exp_ID);
% exp_detect_rest(exp_ID);
check_data(exp_ID);
% decoding_prepare_exp_data(exp_ID);


%%
function check_data(exp_ID)

%%
[~, LFP_ts] = LFP_load(exp_ID,1,"band",'delta');
exp=exp_load_data(exp_ID,'details','path','pos','flight','flight_6m','rest');
FE = [exp.flight.FE];
FE_6m = [exp.flight_6m.FE];
fig=figure;
fig.WindowState = 'maximized';
hold on
plot(exp.pos.proc_1D.ts,exp.pos.proc_1D.pos,'.k')
plot([FE.ts],[FE.pos],'.r')
plot([FE_6m.ts],[FE_6m.pos],'.c')
ylimits = ylim;
for ii_session = 1:length(exp.details.session_names)
    area(exp.details.session_ts(ii_session,:),ylimits([2 2]),'FaceAlpha',0.15);
end
xline(LFP_ts([1 end]),'r','rec');
% rescale_plot_data('x',[1e-6/60 exp.details.session_ts(1)]);
sgtitle(exp_ID,'Interpreter','none');
zoom on
end