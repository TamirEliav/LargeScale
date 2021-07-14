function ripples_MUA_PE_save_to_nlx(exp_ID,forcecalc)

%% defaults
if nargin==1; forcecalc = 0; end

%% get exp info
exp = exp_load_data(exp_ID,'details','path','ripples','MUA','PE');

%% check if file already exist
dir_OUT = 'L:\Analysis\Results\exp\ripples_MUA_PE_nlx';
dir_OUT = fullfile(dir_OUT, exp_ID);
if exist(dir_OUT,'dir')
    if forcecalc
        % delete existing file
        warning('ripples/MUA/PE folder already exist and you chose to override it, deleting old folder!')
        rmdir(dir_OUT,'s')
    else
        warning('zpripple file already exist, use forcecalc to override it!');
        return;
    end
end
mkdir(dir_OUT)

%% write zpripple (.ncs)
file_name = fullfile(dir_OUT, [exp_ID '_zpripple' '.ncs']);
nlx_csc_write(file_name, exp.ripples.zpripple_all, exp.ripples.t);

%% write MUA FR (.ncs)
file_name = fullfile(dir_OUT, [exp_ID '_MUA_zFR' '.ncs']);
nlx_csc_write(file_name, exp.MUA.zFR, exp.MUA.t, exp.MUA.fs);

%% write ripples events (.nev)
for TT = 1:length(exp.ripples.by_TT)
    rpl = exp.ripples.by_TT{TT};
    if isempty(rpl)
        continue;
    end
    ripples_TT_ts = [rpl.peak_ts];
    Filename = fullfile(dir_OUT, [exp_ID '_ripples_byTT' num2str(TT) '.nev']);
    Mat2NlxEV( Filename , 0, 1, [], [1 0 0 0 0 0], ripples_TT_ts);
end
rpl = exp.ripples.all;
ripples_all_ts = [rpl.peak_ts];
Filename = fullfile(dir_OUT, [exp_ID '_ripples' '.nev']);
Mat2NlxEV( Filename, 0, 1, [], [1 0 0 0 0 0], ripples_all_ts);

%% write MUA events (.nev)
Filename = fullfile(dir_OUT, [exp_ID '_MUA' '.nev']);
Mat2NlxEV( Filename, 0, 1, [], [1 0 0 0 0 0], [exp.MUA.events.peak_ts]);

%% write PE events (.nev)
Filename = fullfile(dir_OUT, [exp_ID '_PE_all' '.nev']);
Mat2NlxEV( Filename, 0, 1, [], [1 0 0 0 0 0], [exp.PE.all.peak_ts]);
Filename = fullfile(dir_OUT, [exp_ID '_PE_strong' '.nev']);
Mat2NlxEV( Filename, 0, 1, [], [1 0 0 0 0 0], [exp.PE.strong.peak_ts]);
Filename = fullfile(dir_OUT, [exp_ID '_PE_thr' '.nev']);
Mat2NlxEV( Filename, 0, 1, [], [1 0 0 0 0 0], [exp.PE.thr.peak_ts]);

end