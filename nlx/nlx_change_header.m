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
    error('File OUT already exist, make sure to delete it before running this function');
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
        
        %% update header + data
        header_file = 'Nlx_header_NTT.txt';
        header_new = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');

        ADMaxValue = double(intmax('int16'));
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
        %% read old file
        [Timestamps, Samples, header_old] = ...
            Nlx2MatCSC( file_IN, [1 0 0 0 1], 1, 1, []);
        disp('Old header:');
        disp(header_old);
        
        %% parse and change data
        fs               = sscanf(header_old{contains(header_old,'SamplingFrequency')}, '-SamplingFrequency %g');
        InputInvertedStr = sscanf(header_old{contains(header_old,'InputInverted')},     '-InputInverted %s');

        %% change data and build new header
        ADMaxValue = double(intmax('int16'));
        InputRange = 7000;
        ADBitVolts = InputRange / ADMaxValue / 1e6;
        % convert uVolt to nlx bit
        Samples = Samples ./ ADBitVolts ./ 1e6;
        [~,CSC_name]=fileparts(file_IN);
        SampleFrequencies = repelem(fs,length(Timestamps));
        
        %% update header + data
        header_file = 'Nlg2Nlx_header.txt';
        % read header template
        header_new = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');
        header_new{contains(header_new,'SamplingFrequency')}  = sprintf('-SamplingFrequency %g', fs);
        header_new{contains(header_new,'ADMaxValue')}         = sprintf('-ADMaxValue %d', ADMaxValue);
        header_new{contains(header_new,'InputRange')}         = sprintf('-InputRange %g', InputRange);
        header_new{contains(header_new,'ADBitVolts')}         = sprintf('-ADBitVolts %.24f', ADBitVolts);
        header_new{contains(header_new,'InputInverted')}      = sprintf('-InputInverted %s', InputInvertedStr);
        header_new{contains(header_new,'AcqEntName')}         = sprintf('-AcqEntName %s', CSC_name);
        
        disp('New header:');
        disp(header_new);
        
        %% write new file
        Mat2NlxCSC(file_OUT, 0, 1, 1, [1 0 1 0 1 1], Timestamps, SampleFrequencies, Samples, header_new);
        
    otherwise
        disp(['File type ' EXT ' not supported'])
end
    
end






%%
