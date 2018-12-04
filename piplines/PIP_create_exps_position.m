%% create all exp position data
clear
clc

%% load exp summary and choose exps
exp_t = DS_get_exp_summary();
% exp_t(~(exp_t.batNum == 148),:) = [];

%%
for ii_exp = 1:height(exp_t)
    %%
    exp_ID = exp_t.exp_ID{ii_exp};
    fprintf('%d/%d %s\n', ii_exp, height(exp_t), exp_ID);
    
    %%
    try
        
%     exp_create_position(exp_ID);
%     POS_detect_flight(exp_ID);
%     exp_plot_position(exp_ID);
%     POS_plot_flight(exp_ID);
    exp_compare_pos_vs_csaps(exp_ID)
    
    catch err
        disp(err);
    end
	
    close all
end
