%%
clear
clc

%% load exp summary
exp_t = DS_get_exp_summary();
exp_t(~contains(exp_t.recordingArena, '200m'),:) = [];
exp_t(exp_t.position_data_exist==0,:) = [];
exp_t(exp_t.neural_data_exist==0,:) = [];

%%
% bat_num = 34;
% bat_num = 9861;
bat_num = 2289;
TTs = [1 2 3 4];
% TT_turning_file = 'L:\DATA\0034_Ace\TT_Turning\TT_turning_0034_backup_20180504_final_time_fixed.xlsx';
% TT_turning_file = 'L:\DATA\9861_Somo\TT_turning\TT_turning_9861__20180712_FINAL_time_fixed.xlsx';
TT_turning_file = 'L:\DATA\2289_Sami\TT_turning\TT_turning_2289__20180712_FINAL_time_fixed.xlsx';
TT_data = {};
for TT = TTs
    T = readtable(TT_turning_file,'Sheet',sprintf('TT%d',TT));
    T.time = datetime(T.time,'ConvertFrom','excel','Format','HH:mm:ss');
    T.date.Format = 'dd/MM/yyyy HH:mm:ss';
    T.date = T.date + timeofday(T.time);
    TT_data{TT} = T;
end

%%
exp_t(~ismember(exp_t.batNum, bat_num ),:) = [];
for ii_exp = 1:height(exp_t)
    exp = exp_t(ii_exp,:);
    session_ts = eval(exp.session_ts{1});
    session_ts = session_ts(1);
    exp_ts = exp.date + duration(0,0,session_ts*1e-6);
    exp_ts.Format = 'dd/MM/yyyy HH:mm:ss';
    exp_t.session_time{ii_exp} = exp_ts;
    
%     depth = zeros(size(TTs));
%     for TT=TTs
%         T = TT_data{TT};
%         IX = find(exp_ts > T.date, 1,'last');
%         depth(TT) = T.depth(IX);
%     end
%     exp_t.depth{ii_exp} = char("["+join(""+depth,";")+"]");
end