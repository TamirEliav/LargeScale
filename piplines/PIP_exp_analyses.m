%% exp analysis pipline
clear
clc
%% open log file
log_name_str = ['exp_analysis ' datestr(clock, 'yyyy-mm-dd HH-MM-SS') '.txt'];
log_name_out = fullfile('L:\Analysis\Results\pipelines', log_name_str );
diary off; diary(log_name_out); diary on
% log script code
disp('-------------------------------------------------------------------')
p = mfilename('fullpath')
code_txt = fileread([p '.m'])
disp('-------------------------------------------------------------------')

%%
exp_list = {
'b0148_d170711',...
};

%% load exp summary and choose exps
exp_t = DS_get_exp_summary();
exp_t(~contains(exp_t.recordingArena, '200m'),:) = [];
exp_t(exp_t.position_data_exist==0,:) = [];
exp_t(exp_t.neural_data_exist==0,:) = [];
exp_t(~ismember(exp_t.batNum, [79,148,34,9861,2289] ),:) = [];
% exp_t(~ismember(exp_t.batNum, [34] ),:) = [];
% exp_t(exp_t.date < datetime('08/06/2018','InputFormat','dd/MM/yyyy'),:) = [];
% exp_t(exp_t.date > datetime('17/06/2018','InputFormat','dd/MM/yyyy'),:) = [];
% exp_t(~contains(exp_t.exp_ID, exp_list),:) = [];
% exp_t = flip(exp_t);
exp_t 
whos exp_t 

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
%     exp_create_details(exp_ID);
%     exp=exp_load_data(exp_ID,'details','path');
%     Nlg2Nlx(exp.path.raw,forcecalc);
%     PRE_filter_CSCs(exp_ID, forcecalc);
%     exp_calc_CSC_RAW_stats(exp_ID);
%     PRE_detect_spikes(exp_ID,forcecalc);
%     PRE_calc_write_artifacts(exp_ID);
    
%     bsp_extract_data(exp.path.bsp);
%     exp_create_details(exp_ID);    
%     exp_sync_bsp2nlg(exp_ID);

%     exp_create_position(exp_ID);
%     POS_detect_flight(exp_ID);
%     exp_calc_pos_std_y(exp_ID)
    exp_calc_speed_traj(exp_ID)
%     exp_plot_position(exp_ID);
%     POS_plot_flight(exp_ID);
% 
%     exp_calc_pos_std_y(exp_ID);
%     exp_plot_pos_std_y(exp_ID);

%     wingbeat_detect(exp_ID)
%     wingbeat_plot_map(exp_ID)

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
fprintf('%d/%d exps had error!\nSee details above\n', length(err_list), height(exp_t));
if ~isempty(err_list)
    {err_list.exp_ID}'
end


%% close diary!
diary off



