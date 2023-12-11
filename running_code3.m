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
% data = load('F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & replay distance_gt_10.mat');
data = load('F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115.mat');
IX = find(data.TF);
win_s = 0.5;
params_opt = 11;
epoch_type = 'rest';
replays = [];
for ii_event = IX
    exp_ID = data.data_info.exp_ID(ii_event)
    event_num = data.data_info.evnet_num(ii_event)
    replay = struct('exp_ID',exp_ID,'epoch_type',epoch_type,'params_opt',params_opt,'event_num',event_num);
    replays = [replays replay]
end
struct2table(replays)
for ii_replay = 1:length(replays)
    ii_replay
    replay = replays(ii_replay);
    addFieldsToWorkspace(replay);
    paper_replay_fig_single_replay_example(char(exp_ID),epoch_type,params_opt,event_num,win_s);
    close all
end


%%

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