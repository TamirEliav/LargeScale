function plexon_plx2ntt(PLX_filename, NTT_header_filename, NTT_out_filename)

%% arrange files
if nargin==0
    %%
    [PLX_filename PLX_path] = uigetfile('L:/Analysis/pre_proc/*.plx', 'choose plx file (IN)');
    [NTT_header_filename NTT_header_path] = uigetfile(fullfile(PLX_path,'*.NTT'), 'choose original ntt file (header copy)');
    [NTT_out_filename NTT_out_path] = uiputfile(fullfile(PLX_path,strrep(PLX_filename,'.plx','.NTT')),'choose NTT file to save (OUT)');
    
    PLX_filename = fullfile(PLX_path, PLX_filename);
    NTT_header_filename = fullfile(NTT_header_path, NTT_header_filename);
    NTT_out_filename = fullfile(NTT_out_path, NTT_out_filename);
end
if ~exist(PLX_filename,'file')
    disp('PLX input file does not exist. ABORT!')
    return;
end
if ~exist(NTT_header_filename,'file')
    disp('NTT_header input file does not exist. ABORT!')
    return;
end
if exist(NTT_out_filename,'file')
    disp('NTT output file already exist. ABORT!')
    return;
end

%% file out
% if isempty(NTT_filename_out)
%     NTT_filename_out = strrep(PLX_filename,'.plx','.NTT');
% end
% if exist(NTT_filename_out,'file')
%     fprintf('file out: %s\n',NTT_filename_out);
%     error('file out already exist!!!');
% end

%% read info from plx SPIKES file (we assume the file contains spikes... and not cont...)
[n,names] = plx_chan_names(PLX_filename);
[tscounts,~,~,~] = plx_info(PLX_filename, 1);
tscounts
tscounts(:,[1 3:end]) = []; % remove empty/repeating columns (tetrode? probably nlx->plx convertion problem)
tscounts(~tscounts) = []; % remove empty units
nUnits = length(tscounts);
nWvfrms_per_unit = tscounts;
nWvfrms_total = sum(nWvfrms_per_unit);

fprintf('num units: %d\n', nUnits-1);
fprintf('num waveforms per unit:\n');
for ii_unit = 1:nUnits
    fprintf('%d %d\n', ii_unit-1, nWvfrms_per_unit(ii_unit));
end

%% read waveforms for each unit
% 32 X 4 X nSamples
Samples = zeros(32,4,nWvfrms_total);
Timestamps = zeros(1,nWvfrms_total);
CellNumbers = zeros(1,nWvfrms_total);
counter = 1;
for ii_unit = 1:nUnits
    %%
    unit = ii_unit-1;
    wave = zeros(nWvfrms_per_unit(ii_unit),32,4);
%     clear wave
    for ii_ch = 1:4
        [n, npw, ts, wave(:,:,ii_ch)] = plx_waves_v(PLX_filename, ii_ch, unit);
    end
    wave = permute(wave, [2 3 1]);
    wave = wave * 1e3; % convert to uVolt    
    ts = ts .* 1e6; % change to uSec
    
    IX = counter:(counter+nWvfrms_per_unit(ii_unit)-1);
    Samples(:,:,IX) = wave;
    Timestamps(IX) = ts;
    CellNumbers(IX) = repelem(unit, nWvfrms_per_unit(ii_unit));
    counter = counter + nWvfrms_per_unit(ii_unit);
    
end

%% arrange data to write 
% sort by ts
[~,sort_IX] = sort(Timestamps, 'ascend');
Samples = Samples(:,:,sort_IX);
Timestamps = Timestamps(sort_IX);
CellNumbers = CellNumbers(sort_IX);

%% write ntt file
% load ntt header
if isempty(NTT_header_filename)
    header_file = 'Nlx_header_NTT.txt';
    header = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');
else
    header = Nlx2MatSpike( NTT_header_filename, [0 0 0 0 0], 1, 1, []);
end

ADMaxValue = 32767;
InputRange = max(abs(Samples(:)));
ADC = InputRange / ADMaxValue / 1e6;
Samples = Samples ./ ADC ./ 1e6;
ADC_str = sprintf('%.24f',ADC);
InputRange_str = sprintf('%g',InputRange);
ADC_str_IX = contains(header, 'ADBitVolts');
InputRange_str_IX = contains(header, 'InputRange');
header{ADC_str_IX} = sprintf('-ADBitVolts %s %s %s %s', ADC_str, ADC_str, ADC_str, ADC_str);
header{InputRange_str_IX} = sprintf('-InputRange %s %s %s %s', InputRange_str, InputRange_str, InputRange_str, InputRange_str);

Mat2NlxSpike(NTT_out_filename, 0, 1, [], [1 0 1 0 1 1], ...
    Timestamps, CellNumbers, Samples, header);

disp('NTT file succesfully written!');

end




