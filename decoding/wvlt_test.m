%%
clear
clc

%%
fs=1024*30;
dt=1/fs;
T=5;
% nCh=16;
t=0:dt:T;
f=1e3;
signal = sin(2*pi*f.*t);
signal_gpu = gpuArray(single(signal));
nreps=12*16;
nreps=1;

%% create filter bank
fb = cwtfilterbank( 'SignalLength',length(signal), ...
                    'Wavelet','amor',...
                    'SamplingFrequency',fs,...
                    'VoicesPerOctave',4);
                    
%%
% tic
% for ii=1:nreps
%     [wt1, f1, period1] = cwt(signal,'FilterBank',fb);
% end
% toc

%%
tic
for ii=1:nreps
    [wt2, f2, period2] = cwt(signal_gpu,'FilterBank',fb);
end
toc

%%
figure
cwt(signal,'amor',fs, 'ExtendSignal', 1, 'FrequencyLimits',[2,fs/2], 'VoicesPerOctave',4);

%%
figure
imagesc(period2,f2,abs(wt2))
h=gca;
h.YScale='log';
axis xy



%%
fname = 'L:\Analysis\Results\decoding\raw_datasets\b9861_d180705_pre_proc_dir_1.h5';
wt = h5read(fname, '/inputs/wavelets');
freqs = h5read(fname, '/inputs/fourier_frequencies');
pos = h5read(fname, '/outputs/position');
[BINS, EDGES] = discretize(pos,5);



%






%%