%% exp analysis pipline
clear
clc
%% open log file
log_name_str = ['exp_analysis ' datestr(clock, 'yyyy-mm-dd HH-MM-SS') '.txt'];
log_name_out = fullfile('L:\Analysis\Results\pipelines', log_name_str );
diary off; diary(log_name_out); diary on
% TODO: save script code to log

%%
exp_list = {
'b0148_d170712',...
'b0148_d170713',...
'b0148_d170718',...
'b0148_d170803',...
'b0148_d170807',...
'b0034_d180311',...
'b0034_d180313',...
'b9861_d180604',...
'b9861_d180605',...
'b9861_d180606',...
'b9861_d180607',...
'b9861_d180613',...
'b9861_d180614',...
'b9861_d180618',...
'b9861_d180629',...
'b9861_d180706',...
};

%% load exp summary and choose exps
exp_t = DS_get_exp_summary();
exp_t(~contains(exp_t.recordingArena, '200m'),:) = [];
exp_t(exp_t.position_data_exist==0,:) = [];
exp_t(exp_t.neural_data_exist==0,:) = [];
% exp_t(~ismember(exp_t.batNum, [79,148,34,9861,2289] ),:) = [];
exp_t(~ismember(exp_t.batNum, [34] ),:) = [];
% exp_t(~contains(exp_t.exp_ID, exp_list),:) = [];
exp_t

%%
forcecalc = 1;
err_list = {};
for ii_exp = 1:height(exp_t)
    %%
    exp_ID = exp_t.exp_ID{ii_exp};
    disp(datetime)
    fprintf('%d/%d %s\n', ii_exp, height(exp_t), exp_ID);
    
    %%
try
    exp=exp_load_data(exp_ID,'details','path');
    Nlg2Nlx(exp.path.raw,forcecalc);
%     PRE_filter_CSCs(exp_ID, forcecalc);
%     exp_calc_CSC_RAW_stats(exp_ID)
    
%     bsp_extract_data(exp.path.bsp);
%     exp_create_details(exp_ID);    
%     exp_sync_bsp2nlg(exp_ID);

%     exp_create_position(exp_ID);
%     POS_detect_flight(exp_ID);
%     exp_plot_position(exp_ID);
%     POS_plot_flight(exp_ID);
% 
%     exp_calc_pos_std_y(exp_ID);
%     exp_plot_pos_std_y(exp_ID);
    
catch err
    getReport(err)
    err_list{ii_exp}.exp_ID = exp_ID;
    err_list{ii_exp}.err = err;
end
	
    close all
end
err_list = [err_list{:}];

%% disp all exp with errors
disp('------------------------------------------------')
disp('------------------------------------------------')
disp('------------------------------------------------')
disp('List of errors')
for ii_exp = 1:length(err_list)
    disp(err_list(ii_exp).exp_ID)
    getReport(err_list(ii_exp).err)
end
disp('------------------------------------------------')
fprintf('%d/%d exps had error!\nSee details above:\n', length(err_list), height(exp_t));
{err_list.exp_ID}'



%% close diary!
diary off



