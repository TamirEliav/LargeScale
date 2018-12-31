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
        
        %% read old file
        [Timestamps, ScNumbers, CellNumbers, Features, Samples, header_old] =...
            Nlx2MatSpike(file_IN, [1 1 1 1 1], 1, 1, [] );
        disp('Old header:');
        disp(header_old);
        
        %% update header
        header_file = 'Nlx_header_NTT.txt';
        header_new = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');

        ADMaxValue = 32767;
%         InputRange = 3000;
        InputRange = max(Samples(:));
%         InputRange = prctile(Samples(:),99.9) * 1.2;
        ADC = InputRange / ADMaxValue / 1e6;
        Samples = Samples ./ ADC ./ 1e6;
%         Samples = round(Samples); % TODO: check if rouding is neccesary
        ADC_str = sprintf('%.24f',ADC);
        InputRange_str = sprintf('%g',InputRange);
        ADC_str_IX = contains(header_new, 'ADBitVolts');
        InputRange_str_IX = contains(header_new, 'InputRange');
        header_new{ADC_str_IX} = sprintf('-ADBitVolts %s %s %s %s', ADC_str, ADC_str, ADC_str, ADC_str);
        header_new{InputRange_str_IX} = sprintf('-InputRange %s %s %s %s', InputRange_str, InputRange_str, InputRange_str, InputRange_str);

        disp('New header:');
        disp(header_new);
        
        %% write new file
        Mat2NlxSpike(file_OUT, 0, 1, [], [1 0 1 1 1 1], ...
            Timestamps, CellNumbers, Features, Samples, header_new);
              
    case {'.NCS','.ncs'}
        header_file = 'Nlx_header_NCS.txt';
        
    otherwise
        disp(['File type ' EXT ' not supported'])
end
    
end






%%
