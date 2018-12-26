%% create all exp position data
clear
clc

%% load exp summary and choose exps
exp_t = DS_get_exp_summary();
exp_t(~contains(exp_t.recordingArena, '200m'),:) = [];
exp_t(exp_t.position_data_exist==0,:) = [];
exp_t(~ismember(exp_t.batNum, [79,148,34,9861,2289] ),:) = [];
% exp_t(~contains(exp_t.exp_ID, 'b0034_d180312'),:) = [];

%%
for ii_exp = 1:height(exp_t)
    %%
    exp_ID = exp_t.exp_ID{ii_exp};
    datetime
    fprintf('%d/%d %s\n', ii_exp, height(exp_t), exp_ID);
    
    %%
try
%     exp=exp_load_data(exp_ID);
%     bsp_extract_data(exp.path.bsp);
    exp_create_details(exp_ID);    
    exp_sync_bsp2nlg(exp_ID);

    exp_create_position(exp_ID);
    POS_detect_flight(exp_ID);
    exp_plot_position(exp_ID);
    POS_plot_flight(exp_ID);

    exp_calc_pos_std_y(exp_ID);
    exp_plot_pos_std_y(exp_ID);
    
catch err
    disp(err);
end
	
    close all
end
