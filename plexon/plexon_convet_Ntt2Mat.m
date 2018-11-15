function plexon_convet_Ntt2Mat(file_IN, file_OUT)

%%
% clear 
% clc
% file_IN = 'D:\Tamir\PROJECTS\Plexon\Data\0034\20180313\spikes_b0034_d180313_TT3.NTT';
% file_OUT = [];

%% sanity checks
if ~exist(file_IN,'file')
    warning('Input file does not exist!')
    return
end
if isempty(file_OUT)
    [dir_OUT,NAME,EXT] = fileparts(file_IN);
    file_OUT = fullfile(dir_OUT, [NAME '.mat'] );
end

[dir_OUT,NAME,EXT] = fileparts(file_OUT);
if ~exist(dir_OUT,'dir')
    mkdir(dir_OUT);
end


%% read input file
[Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = ...
    Nlx2MatSpike(file_IN, [1 1 1 1 1], 1, 1, [] );
    
%% convert to plexon matlab format
num_ch = 4; % tetrode!
waveforms = cell(num_ch,1);
timestamps = cell(num_ch,1);
% units = cell(num_ch,1);
Samples = permute(Samples,[2 3 1]);
for ii_ch = 1:num_ch
    waveforms{ii_ch} = squeeze(Samples(ii_ch,:,:));
    timestamps{ii_ch} = (Timestamps-Timestamps(1))'.*1e-6;
end

%% save mat file
save(file_OUT, 'waveforms', 'timestamps')




