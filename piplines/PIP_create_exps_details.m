%% create all exp details
clear
clc

%%
exp_t = DS_get_exp_summary();
for ii_exp = 1:height(exp_t)
    exp_ID = exp_t.exp_ID{ii_exp};
    fprintf('%d/%d %s\n', ii_exp, height(exp_t), exp_ID);
    exp_create_details(exp_ID);
end