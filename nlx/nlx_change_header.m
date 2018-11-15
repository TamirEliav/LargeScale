function nlx_change_header(file_IN, file_OUT)


%% sanity checks
if ~exist(file_IN,'file')
    warning('Input file does not exist!')
    return
end
if isempty(file_OUT)
    [dir_OUT,NAME,EXT] = fileparts(file_IN);
    file_OUT = fullfile(dir_OUT, [NAME '_new_header' EXT] );
end
if exist(file_OUT,'file')
    warning('File OUT already exist, make sure to delete it before running this function');
    return 
end
[dir_OUT,NAME,EXT] = fileparts(file_OUT);
if ~exist(dir_OUT,'dir')
    mkdir(dir_OUT);
end

%%
switch EXT
    case {'.NTT','.ntt'}
        [Timestamps, ScNumbers, CellNumbers, Features, Samples, header_old] =...
            Nlx2MatSpike(file_IN, [1 1 1 1 1], 1, 1, [] );
        disp('Old header:');
        disp(header_old);
        header_file = 'Nlx_header_NTT.txt';
        header_new = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');
        disp('New header:');
        disp(header_new);
        Mat2NlxSpike(file_OUT, 0, 1, [], [1 1 1 1 1 1], Timestamps,...
              ScNumbers, CellNumbers, Features, Samples, header_new);
              
    case {'.NCS','.ncs'}
        header_file = 'Nlx_header_NCS.txt';
        
    otherwise
        disp(['File type ' EXT ' not supported'])
end
    




%% save
