%% look for the source of 1kHz noise in spike waveforms


%%
clear
clc

%%
exp_ID = 'b9861_d180519';
exp=exp_load_data(exp_ID,'path','details');

%% re-filter the raw data
t_start_end = [];
t_start_end = exp_get_sessions_ti(exp_ID,'Sleep2');
clear filter_params
filter_params.type = 'custom';
stopfreq = 550;
passfreq = 650;
filter_params.designfilt_input = {'highpassfir',...
    'StopbandFrequency', stopfreq,...
    'PassbandFrequency', passfreq,...
    'StopbandAttenuation', 60,...
    'PassbandRipple', 1};
file_IN = fullfile(exp.path.nlx, 'CSC1.ncs');
% file_OUT = fullfile(exp.path.nlx, sprintf('CSC1__filtered_type=%s_%d-%dHz_order=%d.ncs',filter_params.type,filter_params.passband,filter_params.order));
file_OUT = fullfile(exp.path.nlx, sprintf('CSC1__filtered_type=%s_stop=%d_pass=%d.ncs',filter_params.type,stopfreq,passfreq));
Nlx_filter_CSC2(file_IN, file_OUT, t_start_end, filter_params);

%% load data from NLX
limits_ts = exp_get_sessions_ti(exp_ID,'Sleep2');
limits_ts = limits_ts(2) + [-2*60*1e6 0];
[signal1, ts, fs, params1] = Nlx_csc_read(fullfile(exp.path.nlx, 'CSC1.ncs'), limits_ts);
[signal2, ts, fs, params2] = Nlx_csc_read(fullfile(exp.path.spikes_raw, 'spikes_b9861_d180519_TT1_ch2.ncs'), limits_ts);
[signal3, ts, fs, params3] = Nlx_csc_read(fullfile(exp.path.nlx, 'CSC1__filtered_600-4000Hz_order=300.ncs'), limits_ts);
[signal4, ts, fs, params4] = Nlx_csc_read(fullfile(exp.path.nlx, 'CSC1__filtered_600-10000Hz_order=50.ncs'), limits_ts);
[signal5, ts, fs, params5] = Nlx_csc_read(fullfile(exp.path.nlx, 'CSC1__filtered_type=custom_stop=750_pass=1000.ncs'), limits_ts);
[signal6, ts, fs, params6] = Nlx_csc_read(fullfile(exp.path.nlx, 'CSC1__filtered_type=custom_stop=500_pass=600.ncs'), limits_ts);
[signal7, ts, fs, params7] = Nlx_csc_read(fullfile(exp.path.nlx, 'CSC1__filtered_type=custom_stop=550_pass=650.ncs'), limits_ts);

%%
clc
median(abs(signal2))
median(abs(signal5))
median(abs(signal6))
median(abs(signal7))
std(signal2)
std(signal5)
std(signal6)
std(signal7)

%%
myfilt = ...
designfilt('bandpassfir',...
    'StopbandFrequency1', 500,...
    'PassbandFrequency1', 600,...
    'PassbandFrequency2', 9000,...
    'StopbandFrequency2', 10000,...
    'StopbandAttenuation1', 60,...
    'StopbandAttenuation2', 10,...
    'PassbandRipple', 1,...
    'SampleRate', fs);
% fvtool(myfilt)
filtord(myfilt)

%%
myfilt = ...
designfilt('highpassfir',...
    'StopbandFrequency', 550,...
    'PassbandFrequency', 650,...
    'StopbandAttenuation', 60,...
    'PassbandRipple', 1,...
    'SampleRate', fs);
% fvtool(myfilt)
filtord(myfilt)

%%
tic
sdf=filtfilt(myfilt,signal1);
toc

%% compare filtered signals
figure
hold on
plot(signal2,'k')
plot(signal5,'r')
% xlim([2672919 2673452])

%% calc power spectrum
win = 2^16;
win = fs*12;
[Pxx1,F1] = pwelch(signal1,win,round(win/2),win,fs);
[Pxx2,F2] = pwelch(signal2,win,round(win/2),win,fs);
[Pxx3,F3] = pwelch(signal3,win,round(win/2),win,fs);
[Pxx4,F4] = pwelch(signal4,win,round(win/2),win,fs);
[Pxx5,F5] = pwelch(signal5,win,round(win/2),win,fs);
[Pxx6,F6] = pwelch(signal6,win,round(win/2),win,fs);
[Pxx7,F7] = pwelch(signal7,win,round(win/2),win,fs);

%%
figure
hold on
plot(F1,log10(Pxx1),'k')
plot(F2,log10(Pxx2),'r')
plot(F5,log10(Pxx5),'g')
plot(F6,log10(Pxx6),'b')
plot(F7,log10(Pxx7),'m')



%% ------------------------------------------------------------------------
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% create new spike shape templates and insert them to the library
clear;clc
cells{1} = load('L:\Analysis\Results\pre_proc\comparing_new_filtering\cells\NEW_filtering_OLD_lib\b9861_d180519_TT1_SS01_cell_spikes.mat');
cells{2} = load('L:\Analysis\Results\pre_proc\comparing_new_filtering\cells\NEW_filtering_OLD_lib\b9861_d180519_TT1_SS02_cell_spikes.mat');
cells{3} = load('L:\Analysis\Results\pre_proc\comparing_new_filtering\cells\NEW_filtering_OLD_lib\b9861_d180519_TT1_SS03_cell_spikes.mat');
cells{1}.best_ch = 2;
cells{2}.best_ch = 4;
cells{3}.best_ch = 3;
load('library_of_acceptable_spike_shapes_new.mat');
new_spike_shapes = [];
for ii_cell = 1:length(cells)
    ch = cells{ii_cell}.best_ch;
    wv = cells{ii_cell}.spikes.waveforms;
    wv = squeeze(mean(wv(:,ch,:),3));
    wv = wv ./ max(abs(wv));
    new_spike_shapes(ii_cell,:) = wv;
end
library_of_acceptable_spike_shapes = [library_of_acceptable_spike_shapes; new_spike_shapes];

figure
plot(library_of_acceptable_spike_shapes')
figure
plot(new_spike_shapes')

save('L:\Analysis\Code\pre_proc\library_of_acceptable_spike_shapes_new_Tamir_20190217', 'library_of_acceptable_spike_shapes');

%%
clear;clc

% load cells from different filtering/spike-shape lib
%         cells{1} = load('L:\Analysis\Results\cells\spikes\b9861_d180519_TT1_SS02_cell_spikes.mat');
%         cells{2} = load('L:\Analysis\Results\cells\spikes\b9861_d180519_TT1_SS01_cell_spikes.mat');
cell_cmpr_opt = 4;
switch cell_cmpr_opt
    case 1
        cells{1} = load('L:\Analysis\Results\pre_proc\comparing_new_filtering\cells\OLD_filtering_OLD_lib\b9861_d180519_TT1_SS01_cell_spikes.mat');
        cells{2} = load('L:\Analysis\Results\pre_proc\comparing_new_filtering\cells\NEW_filtering_OLD_lib\b9861_d180519_TT1_SS02_cell_spikes.mat');
        best_ch = 4;
    case 2
        cells{1} = load('L:\Analysis\Results\pre_proc\comparing_new_filtering\cells\OLD_filtering_OLD_lib\b9861_d180519_TT1_SS02_cell_spikes.mat');
        cells{2} = load('L:\Analysis\Results\pre_proc\comparing_new_filtering\cells\NEW_filtering_OLD_lib\b9861_d180519_TT1_SS01_cell_spikes.mat');
        best_ch = 2;
    case 3
        cells{1} = load('L:\Analysis\Results\pre_proc\comparing_new_filtering\cells\OLD_filtering_OLD_lib\b9861_d180519_TT1_SS03_cell_spikes.mat');
        cells{2} = load('L:\Analysis\Results\pre_proc\comparing_new_filtering\cells\NEW_filtering_OLD_lib\b9861_d180519_TT1_SS03_cell_spikes.mat');
        best_ch = 3;
    case 4
        cells{1} = load('L:\Analysis\Results\pre_proc\comparing_new_filtering\cells\OLD_filtering_OLD_lib\b9861_d180519_TT1_SS04_cell_spikes.mat');
        cells{2} = load('L:\Analysis\Results\pre_proc\comparing_new_filtering\cells\NEW_filtering_OLD_lib\b9861_d180519_TT1_SS04_cell_spikes.mat');
        best_ch = 4;
end
legend_str = {'original filtering (bandpass)';'new filtering (highpass)'};

figure('Units','normalized','Position',[0 0 1 1])
%avg wvfrm
subplot(1,3,1)
hold on
for ii_cell = 1:length(cells)
    wvfrms = cells{ii_cell}.spikes.waveforms;
    avg_wv = squeeze( mean(wvfrms(:,best_ch,:),3) );
    range(avg_wv)
    plot(avg_wv);
%     text(0.8,1.1-0.05*ii_cell,sprintf('n=%d',size(wvfrms,3)),'Units','normalized')
end
legend(legend_str)
title('spikes shape')

% (old) lib corr
load('library_of_acceptable_spike_shapes_new_Ayelet_2018.mat');
subplot(1,3,2)
set(gca,'ColorOrderIndex',1)
hold on
for ii_cell = 1:length(cells)
    wvfrms = cells{ii_cell}.spikes.waveforms;
    rrr = PRE_calc_lib_corr(wvfrms, library_of_acceptable_spike_shapes);
    h=histogram(rrr);
    median(rrr)
end
legend(legend_str)
title('lib correlations (old lib)')

% (new) lib corr
load('library_of_acceptable_spike_shapes_new_Tamir_20190217');
subplot(1,3,3)
set(gca,'ColorOrderIndex',1)
hold on
for ii_cell = 1:length(cells)
    wvfrms = cells{ii_cell}.spikes.waveforms;
    rrr = PRE_calc_lib_corr(wvfrms, library_of_acceptable_spike_shapes);
    h=histogram(rrr);
    median(rrr)
%     text(0.8,0.8,sprintf('n=%d',size(wv,3)),'Units','normalized')
end
legend(legend_str)
title('lib correlations (new lib)')

filename_out = sprintf('cell_cmpr_opt_%d',cell_cmpr_opt);
filename_out = fullfile('L:\Analysis\Results\pre_proc\comparing_new_filtering\cells\', filename_out);
saveas(gcf, filename_out, 'tif');
saveas(gcf, filename_out, 'fig');
close all

%%


