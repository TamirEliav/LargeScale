%%
clear
clc

%%
exp_list_2bats = {
    'b2299_d191202',
    'b2299_d191203',
    'b2299_d191204',
    'b2299_d191205',
    'b2299_d191208',
    'b2299_d191209',
    'b2299_d191210',
    'b2299_d191213',
    };

%%
clear res_all
balls_locs = [];
for ii_exp = 1:length(exp_list_2bats)
    exp_ID = exp_list_2bats{ii_exp};
    exp = exp_load_data(exp_ID,'rest');
    balls_locs(ii_exp,:) = exp.rest.balls_loc;
    filename = fullfile('F:\sequences\decoded_figs\flight\conf_mat\', sprintf('%s_flight_decoding_opt_4.mat',exp_ID));
    res = load(filename);
    res_all(ii_exp) = res.res;
%     res.res.mean_err_prob
%     res.res.pos_err_median_prc
end
100-[res_all.mean_err_prob]'.*100
[res_all.pos_err_median_prc]'.*100
% what distance from the balls do we 
fig_3_replay_inc_range = [25 115];
balls_locs - fig_3_replay_inc_range
mean(abs(balls_locs - fig_3_replay_inc_range),'all')

%%
x = [seqs_pooled.env_size];
thr = 150;
figure
histogram(x)
mean(x(x<thr))
mean(x(x>thr))
median(x(x<thr))
median(x(x>thr))
std(x(x<thr))
std(x(x>thr))

%%
L = [];
for ii_exp = 1:height(T)
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'rest');
    L(ii_exp) = diff(exp.rest.balls_loc);
end


%%
for ii_ex = 1:length(replay_examples)
    ex = replay_examples(ii_ex);
    addFieldsToWorkspace(ex);
    exp = exp_load_data(exp_ID,'ripples','details');
    exp.details.TT_loc{exp.ripples.stats.best_TT}
end







%%

