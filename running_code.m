%%
clear
clc

%% git commmands
% git status
% git add -A
% git commit -a -m "comment"
% git push
% git pull
% git fetch -v --dry-run
% git remote update

%%
clear
clc
exp_t = DS_get_exp_summary();
sessions_names = {'Sleep1';'Behave';'Sleep2'};
% exp_t.sessions_names = {};
T = table();
for ii_exp = 1:height(exp_t)
    exp = exp_t(ii_exp,:);
    
    sessions_ts = [
    exp.session_sleep1_start_ts exp.session_sleep1_end_ts;    
    exp.session_behave_start_ts exp.session_behave_end_ts;
    exp.session_sleep2_start_ts exp.session_sleep2_end_ts;
    ];
    T.sessions_ts{ii_exp} = mat2str(sessions_ts);
    T.sessions_names{ii_exp} = sessions_names;
end
writetable(T,'testTable','FileType','spreadsheet')

%%
clear
clc
exp_t = DS_get_exp_summary();

%% 19/11/2108 disp events (to fill the session ts in exp_summary)
clc
exp_ID = 'b0079_d160902';

[exp_path exp_info] = DS_get_path(exp_ID);

if ~exist(exp_path.nlx,'dir')
    disp('create EVENT files')
    Nlg2Nlx(exp_path.raw);
end
PRE_event_disp(exp_ID);

filename_template = sprintf('L:\\DATA\\%04d_%s\\%s\\nlg\\LOG*.txt',...
    exp_t.batNum(exp_ID),...
    exp_t.batName{exp_ID},...
    datestr(exp_t.date(exp_ID),'yyyymmdd'));
file = dir(filename_template);
fid = fopen(fullfile(file.folder,file.name),'rt');
while ~feof(fid)
    tline = fgetl(fid);
    if contains(tline, {'***','Sent string for remote log'})
        disp(tline)
    end
end


%% 20/11/2018 rename sorted NTT files from Shir
clear
clc
bat = '9861';
dir_main = ['L:\Analysis\pre_proc\SpikeSorting\bat' bat '\2*'];
folders = dir(dir_main);
for ii_folder = 1:length(folders)
    dir_files = fullfile(folders(ii_folder).folder, folders(ii_folder).name,'spikes_NTT');
    files = dir(dir_files);
    files = files(~[files.isdir]);
    file_temp = ['spikes_b' bat '_d' folders(ii_folder).name(3:end) '_'];
    for ii_file = 1:length(files)
        filename_orig = fullfile(files(ii_file).folder, files(ii_file).name);
        filename_new = strrep(filename_orig, 'spikes__', file_temp);
        [status,msg,msgID] = movefile(filename_orig,filename_new)
    end
end















%%




