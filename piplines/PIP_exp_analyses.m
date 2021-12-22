
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
% I still need to run the decoding spikes detection for these days (they
% had some error) - see error details in:
% load('L:\Analysis\Results\pipelines\exp_analysis 2021-06-28 23-58-19.mat')
% exp_list = ...
%     {'b9861_d180622'
%     'b9861_d180628'
%     'b9861_d180704'
%     'b9861_d180705'
%     'b9861_d180706'
%     'b9861_d180710'};

% exp_list = {...
%     'b9861_d180528'
%     'b9861_d180529'
%     'b9861_d180530'
%     'b9861_d180531'
%     'b9861_d180601'
%     'b9861_d180603'
%     'b9861_d180604'};

% exp_list = {...
%     'b0184_d191130'
%     'b0184_d191201'
%     'b9861_d180526'
%     'b9861_d180527'
%     'b0034_d180313'
%     };

% list of days to run decoding analysis
exp_list = {...
    'b0184_d191126'
    'b0184_d191127'
    'b0184_d191128'
    'b0184_d191129'
    'b0184_d191130'
    'b0184_d191201'
    'b0184_d191202'
    'b0184_d191203'
    'b0184_d191204'
    'b0184_d191205'
    'b0184_d191208'
    'b0184_d191209'
    
%     'b9861_d180519'
%     'b9861_d180521'
%     'b9861_d180522'
%     'b9861_d180523'
%     'b9861_d180524'
%     'b9861_d180525'
%     'b9861_d180526'
%     'b9861_d180527'
%     'b9861_d180529'
%     'b9861_d180530'
%     
%     'b0034_d180305'
%     'b0034_d180306'
%     'b0034_d180308'
%     'b0034_d180310'
%     'b0034_d180311'
%     'b0034_d180312'
%     'b0034_d180313'
%     'b0034_d180314'
%     'b0034_d180315'

%     'b2289_d180518'
%     'b2289_d180520'
    
%     'b0079_d161003'
% TODO: run filtering for bat 79    
% TODO: look for additional good days from bat 79    
    
%     'b0148_d170606'
%     'b0148_d170607'
%     'b0148_d170608'
%     'b0148_d170625' % fig 1 sorting example from the paper
%     'b0148_d170626'
%     'b0148_d170627'
% %     'b0148_d170628' % had error 
%     'b0148_d170703'
%     'b0148_d170718'
%     'b0148_d170720'
%     'b0148_d170723'
%     'b0148_d170801'
%     'b0148_d170802'
% %     'b0148_d170803' % had error 
%     'b0148_d170806'
%     'b0148_d170807'
    
    };

% list of days that had problems running the ripples/MUA/PE code
% exp_list = {
%     'b0034_d180305'
%     'b0034_d180306'
%     'b9861_d180521'
%     'b9861_d180524'
%     'b9861_d180526'
% };

exp_list = {...
    'b9861_d180526'
    'b9861_d180527'
    'b0184_d191130'
    'b0184_d191201'
    };

%% load exp summary and choose exps
exp_t = DS_get_exp_summary();
% exp_t(~contains(exp_t.recordingArena, '200m'),:) = [];
exp_t(exp_t.position_data_exist==0,:) = [];
exp_t(exp_t.neural_data_exist==0,:) = [];
% exp_t(~ismember(exp_t.batNum, [79,148,34,9861,2289] ),:) = [];
% exp_t(~ismember(exp_t.batNum, [9861] ),:) = [];
exp_t(~ismember(exp_t.batNum, [184] ),:) = [];
% exp_t(~contains(exp_t.TT_loc,{'CA1','CA3'}),:) = [];
% exp_t(exp_t.date < datetime('08/06/2018','InputFormat','dd/MM/yyyy'),:) = [];
% exp_t(exp_t.date > datetime('17/06/2018','InputFormat','dd/MM/yyyy'),:) = [];
% exp_t(contains(exp_t.exp_ID, exp_list),:) = [];
% exp_t(~contains(exp_t.exp_ID, exp_list),:) = [];
% exp_t = flip(exp_t);
exp_t 
whos exp_t 

%%
forcecalc = 0;
err_list = {};
for ii_exp = 1:height(exp_t)
    %%
    exp_ID = exp_t.exp_ID{ii_exp};
    disp(datetime)
    fprintf('%d/%d %s\n', ii_exp, height(exp_t), exp_ID);
    
    %%
try
%     exp_create_details(exp_ID);
%     exp=exp_load_data(exp_ID,'details','path','rest','MUA');
%     util_fix_ncs(exp.path.nlx);
%     Nlg2Nlx(exp.path.raw,forcecalc);
%     PRE_filter_CSCs(exp_ID, forcecalc);
%     PRE_filter_LFP_bands(exp_ID, forcecalc);
%     exp_calc_CSC_RAW_stats(exp_ID);
%     PRE_detect_spikes(exp_ID,forcecalc);
%     PRE_calc_write_artifacts(exp_ID);

%     bsp_extract_data(exp.path.bsp);
%     exp_create_details(exp_ID);
%     exp_sync_bsp2nlg(exp_ID);

%     exp_create_position(exp_ID);
%     POS_detect_flight(exp_ID);
%     exp_calc_pos_std_y(exp_ID)
%     exp_calc_speed_traj(exp_ID)
%     exp_plot_position(exp_ID);
%     POS_plot_flight(exp_ID);

%     POS_test_join_short_flights_with_gap(exp_ID)

%     exp_calc_speed_traj(exp_ID);
%     exp_calc_pos_std_y(exp_ID);
%     exp_plot_pos_std_y(exp_ID);
%     exp_plot_xyz_pos(exp_ID)

%     wingbeat_detect(exp_ID)
%     wingbeat_plot_map(exp_ID)

%     exp_detect_rest(exp_ID);

%     decoding_detect_spikes(exp_ID,forcecalc);
%     decoding_prepare_exp_data(exp_ID);

%     ripples_detect(exp_ID);
%     MUA_detect(exp_ID);
%     PE_plot_ripples_vs_MUA(exp_ID); % I put the detection here...
% %     ripples_MUA_PE_save_to_nlx(exp_ID,forcecalc);

%     ripples_trigger_LFP(exp_ID);
%     ripples_trigger_MUA(exp_ID);
%     ripples_xcorr(exp_ID)
    

%     for params_opt = 4%1:6
%         decoding_plot_flight_conf_mat(exp_ID, params_opt);
%         decoding_plot_flight_posterior(exp_ID, params_opt);
%     end
    
% % %     epoch_type = 'sleep';
% % % %     epoch_type = 'rest';
% % % %     epoch_type = 'flight';
% % %     params_opts = [8:14];
% % % %     event_type = 'PE';
% % %     event_type = 'posterior';
% % %     for params_opt = params_opts
% % %         fprintf('params_opt: %d\n', params_opt);
% % % %         decoding_plot_MAP(exp_ID, epoch_type, params_opt);
% % % %         decoding_detect_posterior_events(exp_ID, epoch_type, params_opt);
% % % %         decoding_seq_quantify(exp_ID, epoch_type, params_opt, event_type);
% % % %         decoding_seq_quantify_plot(exp_ID, epoch_type, params_opt, event_type);
% % % %         close all
% % % %         decoding_plot_PE_posterior(exp_ID, epoch_type, params_opt, event_type);
% % %         decoding_xcorr_ripples_MUA_PE_vs_posterior_events(exp_ID, epoch_type, params_opt);
% % %     end
% % % %     decoding_compare_replay_speeds(exp_ID, epoch_type, params_opts, event_type);
% % % %     close all

catch err
    getReport(err)
    err_list{ii_exp}.exp_ID = exp_ID;
    err_list{ii_exp}.err = err;
end

    save(strrep(log_name_out,'.txt','.mat'), 'err_list');

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





%%



