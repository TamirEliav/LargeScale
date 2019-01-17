% plexon_create_demo_matlab_files

%%
clear
clc

%% params
rng(123);
% cells_FRs = [10 5 3 8 1];
cells_FRs = [10 10];
num_cells = length(cells_FRs);
% waveforms_peaks_val = 200+300.*rand(num_cells,4);
waveforms_peaks_val = [ 100 100 100 100;...
                        200 200 200 200];
waveforms_peaks_std = 10; 
Ttotal = 15*60;
Fs = 32000;
Ts = 1/Fs;

%% spike canonical waveform
load('L:\Analysis\Code\pre_proc\library_of_acceptable_spike_shapes.mat');
spike_waveform = mean(library_of_acceptable_spike_shapes);

%% create demo data

% create timestamps
t = 1:Ts:Ttotal;
clear cells
ch_signals = zeros(4,length(t));
global_amp_change = 100.*sin(0.1*t);
rng(123);
for ii_cell = 1:num_cells
    num_spikes = cells_FRs(ii_cell)*Ttotal;
    spikes_IX = randsample(1:length(t), num_spikes);
    % remove overlapping spikes
    min_ISI = 0.002/Ts;
    spikes_IX( diff(spikes_IX)< min_ISI) = [];
    cells{ii_cell}.spikes_ts = t(spikes_IX);
    cells{ii_cell}.spikes_IX = spikes_IX;
    for ii_ch = 1:4
        ch_signals(ii_ch, spikes_IX) = ...
            ch_signals(ii_ch, spikes_IX) +...
            waveforms_peaks_val(ii_cell,ii_ch) +...
            waveforms_peaks_std.*randn(1,length(spikes_IX)) +...
            global_amp_change(spikes_IX);
    end
end

%% convolve the spike train with the waveform
for ii_ch = 1:4
    ch_signals(ii_ch, :) = conv(ch_signals(ii_ch, :), spike_waveform, 'same');
end

%% extract spikes waveforms
trigger_IX = [];
for ii_cell = 1:num_cells
    spikes_IX = cells{ii_cell}.spikes_IX;
    trigger_IX = [trigger_IX spikes_IX];
end
trigger_IX = sort(trigger_IX, 'ascend');
spikes_ts = t(trigger_IX );
spikes_waveforms = zeros(4, length(trigger_IX), 32);
for ii_ch = 1:4
    signal = ch_signals(ii_ch, :);
    win = [-16 15];
    triggered_signal = trigger_signal_by_IX(signal, trigger_IX, win);
    spikes_waveforms(ii_ch, :, :) = triggered_signal;
end

%% plot
sdf = max(spikes_waveforms,[],3);
figure
for ii_ch = 1:4
    subplot(2,2,ii_ch)
    plot(sdf(ii_ch,:),'.')
    xlabel('time')
    ylabel('peak')
end
suptitle('peak changed over time');
% saveas(gcf, 'D:\Tamir\PROJECTS\Plexon\Data\MATLAB\peak changed over time','tif');

%% export demo data
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% export to mat files

format_opt = 1;

% prepare plexon format
switch format_opt 
    case 1
        waveforms = cell(4,1);
        timestamps = cell(4,1);
        for ii_ch = 1:4
            waveforms{ii_ch} = squeeze(spikes_waveforms(ii_ch,:,:));
            timestamps{ii_ch} = spikes_ts';
        end
    case 2
        waveforms = cell(1,1);
        timestamps = cell(1,1);
        waveforms{1} = spikes_waveforms;
        timestamps{1} = spikes_ts';
end

% save
save('D:\Tamir\PROJECTS\Plexon\Data\MATLAB\contdata.mat', 'ch_signals')
save('D:\Tamir\PROJECTS\Plexon\Data\MATLAB\spikes.mat', 'waveforms', 'timestamps')

%% export to NLX files (NCS)
if 1
for ii_ch = 1:4
    file_name = ['D:\Tamir\PROJECTS\Plexon\Data\MATLAB\CSC\' sprintf('TT1_ch%d.ncs', ii_ch)];
    signal = ch_signals(ii_ch, :);
    time_offset = 60*60*10;
    time_offset = 0;
    ts_usec = (t+time_offset)*1e6;
    fs = Fs;
    
    block_size = 512;
    rsdl = mod(length(signal), block_size);
    signal(end:end-rsdl+1) = [];
    ts_usec(end:end-rsdl+1) = [];
    
    nlx_csc_write(file_name, signal, ts_usec', fs);
end
end
%% export to NLX files (NTT)
if 1
Timestamps = t(trigger_IX)*1e6;
Samples = permute(spikes_waveforms,[3 1 2]);
file_name = ['D:\Tamir\PROJECTS\Plexon\Data\MATLAB\TT1_spikes.NTT'];
header_file = 'plexon_Nlx_header.txt';
header = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');
Mat2NlxSpike(file_name, 0, 1, [], [1 0 0 0 1 1], Timestamps, Samples, header);
end

%% nlx test
% % % % % % file_name = 'D:\Tamir\PROJECTS\Plexon\Data\MATLAB\test.ncs';
% % % % % % ts = 1:Ts:5*60;
% % % % % % ts = ts*1e6;
% % % % % % rng(123);
% % % % % % signal = randn(size(ts));
% % % % % % nlx_csc_write(file_name, signal, ts, Fs);







%%
