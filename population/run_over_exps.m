%%
clear all
clc

%% get exp list
exp_list_source = 2; % 1 - list; 2- xls
switch exp_list_source
    case 1
        exp_IDs = {...
            'b0079_d160921',...
            'b0148_d170611',...
            'b0148_d170618',...
            };
    case 2
        exp_data_file = 'L:\Analysis\Code\inclusion_lists\recordings_summary.xlsx';
        exp_t = readtable(exp_data_file, 'Sheet', 'Experiments', 'ReadRowNames',1);
        exp_t = exp_t(exp_t.batNum==9861,:);
        exp_IDs = exp_t.Properties.RowNames;
end

% exp_IDs = exp_IDs(1);

%% run over days
for ii_exp = 1:length(exp_IDs)
try    
    %% get exp details
    exp_ID = exp_IDs{ii_exp};
    [exp_path exp_info] = DS_get_path(exp_ID);
    disp(exp_ID)
    
    %%
%     PRE_run_exp(exp_ID);
% % % % 
% % % % if ~exist(exp_path.position, 'dir')
% % % %     bsp_extract_data(exp_path.bsp)
% % % %     Nlg2Nlx(exp_path.raw) % TODO: insert params from here, rather than change them inside Nlg2Nlx function, also consider to remove some params (like ref channel...). this code should ONLY be a reader/parser code!
% % % %     PRE_sync_bsp_to_nlg(exp_path.bsp, exp_path.nlx, exp_path.sync);
% % % %     POS_pre_process(exp_ID);
% % % % end

    %%

    %%
%     PRE_calc_flight_rhythm(exp_ID);
    LFP_filter_theta(exp_ID);

    %%
%     POS_detect_flight(exp_ID)
%     behave_detect_rest(exp_ID,1);
%     TT = 2;
%     plot_main_low_freq_LFP(exp_ID, TT)


    %%
%     wingbeat_detect( exp_ID );
%     wingbeat_calc_phase_map(exp_ID)
%     wingbeat_calc_map( exp_ID );
%     wingbeat_plot_map( exp_ID );

catch err
    disp(err)
end

end















%%