function nlx_csc_write(file_name, signal, ts_usec, fs, header)

%% dummy inputs
% % file_name = 'D:\Tamir\test\spikes_TT1_ch1.ncs';
% % fs = 1000.0;
% % ts = (0:1/fs:18) + 20;
% % ts_usec = round(ts.*1e6);
% % freq = linspace(10,100,length(ts));
% % amp = linspace(1000,2000,length(ts));
% % signal = amp.*sin(2*pi*freq.*ts);

%% sanity checks
% TODO: add check of file_name end with .ncs
% TODO: check if file exist (because there are problems with nlx files if they exist - because of appending bug they have...)
if exist(file_name,'file')
    error('File exist, cannot override!');
end

%% parse input
if nargin<3
    error('Not enough input variables')
end
if nargin<4
    fs = 1e6/median(diff(ts_usec));
end
% % if nargin<5
% %     is_header_exist = 0;
% % else if isfield(params,'header')
% %         is_header_exist = 1;
% %     else
% %         is_header_exist = 0;
% %     end
% % end

%% arrange data
signal = signal(:)';
ts_usec = ts_usec(:)';
block_size = 512;
padding_n = block_size - mod(length(signal),block_size);
padding = zeros(1,padding_n);
signal_blocks = reshape( [signal padding], block_size, []);
num_blocks = size(signal_blocks,2);
blocks_ts = round(ts_usec(1:block_size:end));
SampleFrequencies = repmat(fs,1,num_blocks);

%% Header params
if nargin<5
    header_filename = 'Nlg2Nlx_header.txt';
    header = textread(header_filename, '%s', 'delimiter', '\n', 'whitespace', '');
end

% change fields that are signal specific
ADMaxValue = 32767;
InputRange = max(signal(:));
ADBitVolts = InputRange / ADMaxValue / 1e6;
signal_blocks = signal_blocks  ./ ADBitVolts ./ 1e6;

header{contains(header,'SamplingFrequency')}  = sprintf('-SamplingFrequency %g', fs);
header{contains(header,'ADMaxValue')}         = sprintf('-ADMaxValue %d', ADMaxValue);
header{contains(header,'InputRange')}         = sprintf('-InputRange %g', InputRange);
header{contains(header,'ADBitVolts')}         = sprintf('-ADBitVolts %.24f', ADBitVolts);
% header{contains(header,'InputInverted')}      = sprintf('-InputInverted %s', InputInvertedStr);

%% write file
AppendToFileFlag = 0;
ExportMode = 1;
ExportModeVector = 0;
FieldSelectionFlags = [1 0 1 0 1 1];

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
            SampleFrequencies,...
            signal_blocks,...
            header);

end
