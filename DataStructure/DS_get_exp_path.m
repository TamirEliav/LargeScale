function exp_path = DS_get_exp_path(exp_ID)
%% 
% TODO: change top only one output exp with exp_path inside as a field

%% load exp details
exp = exp_load_data(exp_ID,'details');

%% generaet data path strings
date_str = datestr(exp.details.date, 'yyyymmdd');
batNum_str = sprintf('%04d',exp.details.batNum);
batName_str = exp.details.batName;
batNumName_str = [batNum_str '_' batName_str];

% base folders
data_dir = 'L:\DATA';
preproc_dir = 'L:\Analysis\pre_proc';
if ismember(exp.details.batNum,[2299 2336])
    data_dir = 'G:\Tamir\work\PROJECTS\LargeScale\DATA';
    preproc_dir = 'G:\Tamir\work\PROJECTS\LargeScale\Analysis\pre_proc';
end

% RAW
exp_path.raw          = fullfile(data_dir, batNumName_str, date_str);
exp_path.nlg          = fullfile(data_dir, batNumName_str, date_str, 'nlg');
exp_path.nlx          = fullfile(data_dir, batNumName_str, date_str, 'nlx');
exp_path.bsp          = fullfile(data_dir, batNumName_str, date_str, 'bsp', 'client');
exp_path.bsp_tag_file = fullfile(data_dir, batNumName_str, date_str, 'bsp', 'client',['bsp_pos_tag_' num2str(exp.details.bsp_tag_ID) '.mat']);
exp_path.audio        = fullfile(data_dir, batNumName_str, date_str, 'audio');
exp_path.bin          = fullfile(data_dir, batNumName_str, date_str, 'bin');

% CALIB
exp_path.calib_landmarks_file = 'L:\DATA\9892_Rocky\calib\landmarks.xlsx';
exp_path.calib_tunnel_file = 'L:\DATA\9892_Rocky\calib\calib_tunnel.mat';

% PRE-PROCESS
exp_path.sync = fullfile(data_dir, batNumName_str, date_str, 'sync');
exp_path.spikes_raw = fullfile(preproc_dir, batNum_str, date_str, 'spikes_raw');
exp_path.spikes_detection = fullfile(preproc_dir, batNum_str, date_str, 'spikes_detection');
exp_path.decoding_spikes_detection = fullfile(preproc_dir, batNum_str, date_str, 'decoding_spikes_detection');
exp_path.spikes_sorting = fullfile(preproc_dir, batNum_str, date_str, 'spikes_sorting');
exp_path.LFP = fullfile(preproc_dir, batNum_str, date_str, 'LFP');
exp_path.LFP_bands = fullfile(preproc_dir, batNum_str, date_str, 'LFP_bands');
exp_path.LFP_theta = fullfile(preproc_dir, batNum_str, date_str, 'LFP_theta');
exp_path.LFP_reref = fullfile(preproc_dir, batNum_str, date_str, 'LFP_reref');
exp_path.ripples = fullfile(preproc_dir, batNum_str, date_str, 'ripples');
exp_path.ripples_reref = fullfile(preproc_dir, batNum_str, date_str, 'ripples_reref');

% RESULTS

end