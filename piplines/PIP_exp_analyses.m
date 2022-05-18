%% exp analysis pipline
close all
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
% list of days without pos/flight structs for bats 9845/194
% exp_list = {
%     'b0194_d180429'
% %     'b9845_d170212'
%     'b9845_d170213'
% %     'b9845_d170215'
% %     'b9845_d170606'
%     'b9845_d170608'
%     'b9845_d170615'
%     }

%% load exp summary and choose exps
% exp_t = DS_get_exp_summary();
% exp_list = decoding_get_inclusion_list();
% % exp_t(~contains(exp_t.recordingArena, ['200m']),:) = [];
% exp_t(~contains(exp_t.recordingArena, {'200m','120m'}),:) = [];
% exp_t(exp_t.position_data_exist==0,:) = [];
% exp_t(exp_t.neural_data_exist==0,:) = [];
% % exp_t(~ismember(exp_t.batNum, [79,148,34,9861,2289] ),:) = []; % 5 wild bats in the paper
% % exp_t(~ismember(exp_t.batNum, [9861 34 2289 79 148 184 194 2382 9845] ),:) = [];
% exp_t(~ismember(exp_t.batNum, [9861 34 2289 148 184 194 2382] ),:) = []; % bats good for decoding
% % exp_t(~ismember(exp_t.batNum, [194] ),:) = [];
% % exp_t(~ismember(exp_t.batNum, [184] ),:) = [];
% % exp_t(~ismember(exp_t.batNum, [9845] ),:) = [];
% exp_t(~contains(exp_t.TT_loc,{'CA1','CA3'}),:) = [];
% % exp_t(exp_t.date < datetime('08/06/2018','InputFormat','dd/MM/yyyy'),:) = [];
% exp_t(ismember(exp_t.exp_ID,'b0194_d180429'),:)=[]; % bad position data
% % exp_t(~ismember(exp_t.exp_ID,'b2382_d190814'),:)=[]; % rest decoding need to finish (session with long rest epoch)
% % exp_t(ismember(exp_t.exp_ID,exp_list),:)=[];
% exp_t(~ismember(exp_t.exp_ID,exp_list),:)=[];
% 
% exp_t 
% whos exp_t 

%% load exp summary and choose exps - preproc 2bats data
exp_t = DS_get_exp_summary();
% exp_t(~ismember(exp_t.batNum, [2336] ),:) = []; % Ayelet
exp_t(~ismember(exp_t.batNum, [2299] ),:) = []; % shaked
% exp_t(~ismember(exp_t.batNum, [2382] ),:) = []; % shir
% exp_t(exp_t.date >= datetime('23/06/2019','InputFormat','dd/MM/yyyy'),:) = [];
% exp_t = flip(exp_t);
exp_t 
whos exp_t 

%%
% clc
% exp_ID = exp_t.exp_ID{end-0}
% exp_create_details(exp_ID);
% close all
% check_data(exp_ID);
% exp_detect_rest(exp_ID);

%% run some pop analysis
% exp_list = exp_t.exp_ID;
% decoding_flight_pop_analysis(exp_list);

%%
forcecalc = 0;
err_list = {};
for ii_exp = 1:height(exp_t)
    %%
    exp_ID = exp_t.exp_ID{ii_exp};
    fprintf('%d/%d %s \t\t (start run: %s)\n', ii_exp, height(exp_t), exp_ID, datetime);
    
    %%
try
    exp_create_details(exp_ID);
%     exp=exp_load_data(exp_ID,'details','path');
%     exp=exp_load_data(exp_ID,'details','path','pos','flight');
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

    decoding_detect_spikes(exp_ID,forcecalc);
%     decoding_prepare_exp_data(exp_ID);

%     ripples_detect(exp_ID);
%     MUA_detect(exp_ID);
%     PE_plot_ripples_vs_MUA(exp_ID); % I put the detection here...
%     ripples_MUA_PE_save_to_nlx(exp_ID,forcecalc);
%     exp_calc_MUA_FR_map(exp_ID);

%     ripples_trigger_LFP(exp_ID);
%     ripples_trigger_MUA(exp_ID);
%     ripples_xcorr(exp_ID);
    

%     epoch_type = 'flight';
%     for params_opt = 4%1:6
%         decode = decoding_load_data(exp_ID, epoch_type, params_opt);
%         decoding_plot_flight_conf_mat(exp_ID, params_opt);
%         decoding_plot_MAP(decode);
% %         decoding_plot_flight_posterior(exp_ID, params_opt);
%     end
    

% %     epoch_type = 'sleep';
%     epoch_type = 'rest';
% %     params_opts = [8:14];
% %     params_opts = [11];
%     params_opts = [11];
% %     event_type = 'PE';
%     event_type = 'posterior';
%     for params_opt = params_opts
%         fprintf('params_opt: %d\n', params_opt);
%         decode = decoding_load_data(exp_ID, epoch_type, params_opt);
% %         decoding_plot_MAP(decode);
% %         decoding_detect_posterior_events(decode);
% %         decoding_seq_quantify(decode, event_type);
% %         decoding_seq_quantify_add_info(exp_ID, epoch_type, params_opt , event_type);
% %         decoding_seq_quantify_plot(exp_ID, epoch_type, params_opt, event_type); 
%         decoding_plot_PE_posterior(decode, event_type);
% %         decoding_plot_session_seqs(exp_ID, epoch_type, params_opt, event_type);
% %         decoding_xcorr_ripples_MUA_PE_vs_posterior_events(decode);
%         close all
%     end
% %     decoding_compare_replay_speeds(exp_ID, epoch_type, params_opts, event_type);

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
function check_data(exp_ID)

%%
% [~, LFP_ts] = LFP_load(exp_ID,1,"band",'delta');
exp=exp_load_data(exp_ID,'details','path','pos','flight','flight_6m','rest');
FE = [exp.flight.FE];
% FE = [exp.flight_6m.FE];
fig=figure;
fig.WindowState = 'maximized';
hold on
plot(exp.pos.proc_1D.ts,exp.pos.proc_1D.pos,'.k')
plot([FE.ts],[FE.pos],'.r')
ylimits = ylim;
for ii_session = 1:length(exp.details.session_names)
    area(exp.details.session_ts(ii_session,:),ylimits([2 2]),'FaceAlpha',0.15);
end
% xline(LFP_ts([1 end]),'r','rec');
% rescale_plot_data('x',[1e-6/60 exp.details.session_ts(1)]);
sgtitle(exp_ID,'Interpreter','none');
zoom on
end


