function nlx_csc_write(file_name, signal, ts_usec, fs)

%% dummy inputs
% % file_name = 'D:\Tamir\test\spikes_TT1_ch1.ncs';
% % fs = 1000.0;
% % ts = (0:1/fs:18) + 20;
% % ts_usec = round(ts.*1e6);
% % freq = linspace(10,100,length(ts));
% % amp = linspace(1000,2000,length(ts));
% % signal = amp.*sin(2*pi*freq.*ts);

%%
% TODO: add check of file_name end with .ncs
% TODO: check if file exist (because there are problems with nlx files if they exist - because of appending bug they have...)

%% arrange data
signal = signal(:)';
ts_usec = ts_usec(:)';
padding_n = 512 - mod(length(signal),512);
padding = zeros(1,padding_n);
signal_blocks = reshape( [signal padding], 512, []);
num_blocks = size(signal_blocks,2);
blocks_ts = round(ts_usec(1:512:end));
SampleFrequencies = repmat(fs,1,num_blocks);
% Header = {...
%     '-FileType: CSC ';...
%     '-RecordSize 1044 ';...
% ...%     '-ADCounts 1';...
%     '-ADBitVolts 1.0';...
%     ['-SamplingFrequency ' sprintf('%.30f',fs)];...
%     };

header_file = which('Nlg2Nlx_header.txt');
Header = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');
Header = strrep(Header,'-SamplingFrequency', ['-SamplingFrequency ' char(vpa(fs))]);
% Header = strrep(Header,'-RecordSize', ['-RecordSize ' num2str(512*num_blocks)]);


%% write file
AppendToFileFlag = 0;
ExportMode = 1;
ExportModeVector = 0;
FieldSelectionFlags = [1 0 0 0 1 1];

% make sure folder exists
[PATHSTR,NAME,EXT] = fileparts(file_name);
if ~exist(PATHSTR, 'dir')
    mkdir(PATHSTR)
end

Mat2NlxCSC( file_name,...
            AppendToFileFlag,...
            ExportMode, ExportModeVector,...
            FieldSelectionFlags,...
            blocks_ts,...
...%             SampleFrequencies,...
            signal_blocks,...
            Header);

end