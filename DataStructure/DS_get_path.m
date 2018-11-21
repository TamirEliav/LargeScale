function [exp_path exp_info] = DS_get_exp_path(exp_ID)
%% 
% TODO: change top only one output exp with exp_path inside as a field

%%
exp_t = DS_get_exp_summary();
exp = exp_t(exp_ID,:);
exp = table2struct(exp);
exp.exp_ID = exp_ID;

%% generaet data path strings
% date_str = datestr(datevec(char(exp.date),'dd/mm/yyyy'), 'yyyymmdd' );
% date_str = datestr(datevec(exp.date,'dd/mm/yyyy'), 'yyyymmdd' );
date_str = datestr(exp.date, 'yyyymmdd');
batNum_str = sprintf('%04d',exp.batNum);
% batNum_str = exp.batNum;

% RAW
exp_path.raw = fullfile('L:\DATA', [batNum_str '_' exp.batName], date_str);
exp_path.nlg = fullfile('L:\DATA', [batNum_str '_' exp.batName], date_str, 'nlg');
exp_path.nlx = fullfile('L:\DATA', [batNum_str '_' exp.batName], date_str, 'nlx');
exp_path.bsp = fullfile('L:\DATA', [batNum_str '_' exp.batName], date_str, 'bsp', 'client');
exp_path.audio = fullfile('L:\DATA', [batNum_str '_' exp.batName], date_str, 'audio');
% exp_path.landmarks_file = fullfile('L:\Analysis\Code\calib', exp.LandmarksFile);
exp_path.calib_landmarks_file = 'L:\DATA\9892_Rocky\calib\landmarks.xlsx';
% exp_path.tunnel_calib_file = fullfile('L:\Analysis\Code\calib', exp.tunnelCalib);
exp_path.calib_tunnel_file = 'L:\DATA\9892_Rocky\calib\calib_tunnel.mat';

% PRE-PROCESS
exp_path.sync = fullfile('L:\DATA', [batNum_str '_' exp.batName], date_str, 'sync');
exp_path.spikes_raw = fullfile('L:\Analysis\pre_proc', batNum_str, date_str, 'spikes_raw');
exp_path.spikes_detection = fullfile('L:\Analysis\pre_proc', batNum_str, date_str, 'spikes_detection');
exp_path.spikes_clustring = fullfile('L:\Analysis\pre_proc', batNum_str, date_str, 'spikes_clustering');
exp_path.LFP = fullfile('L:\Analysis\pre_proc', batNum_str, date_str, 'LFP');
exp_path.LFP_theta = fullfile('L:\Analysis\pre_proc', batNum_str, date_str, 'LFP_theta');
exp_path.LFP_reref = fullfile('L:\Analysis\pre_proc', batNum_str, date_str, 'LFP_reref');
exp_path.ripples = fullfile('L:\Analysis\pre_proc', batNum_str, date_str, 'ripples');
exp_path.ripples_reref = fullfile('L:\Analysis\pre_proc', batNum_str, date_str, 'ripples_reref');
exp_path.position = fullfile('L:\Analysis\pre_proc', batNum_str, date_str, 'position');
exp_path.cells_data = fullfile('L:\Analysis\pre_proc', batNum_str, date_str, 'cells_data');

% RESULTS


%% just rename...
exp_info = exp;

end