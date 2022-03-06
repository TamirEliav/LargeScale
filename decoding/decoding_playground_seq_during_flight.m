%%
clear
clc

%%
exp_ID='b0184_d191208';
epoch_type ='flight';
params_opt = 10;
exp = exp_load_data(exp_ID,'details','path','pos','flight');
decode = decoding_load_data(exp_ID, epoch_type, params_opt);

%% load LFP (theta filtered)
TT_to_use = find(contains(exp.details.TT_loc, {'CA1','CA3'}));
[LFP, ts, fs, params, ch_valid] = LFP_load(exp_ID,TT_to_use,'band','theta');
wingbeat = nanmean(LFP,[2 3]);
wingbeat_phase = angle(hilbert(wingbeat));

%%
t = decode.time;
pos_real = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos, t, 'linear');
pos_estimate = decode.MAP_pos;
a = -decode.params.pos_bin_size/2;
b =  decode.params.pos_bin_size/2;
noise = a+(b-a).*rand(size(pos_estimate));
% pos_estimate = pos_estimate + noise;
% pos_estimate = smoothdata(pos_estimate,'gaussian',50);
pos_err = pos_estimate-pos_real;
pos_err_smooth = smoothdata(pos_err,'gaussian',0.1*decode.Fs);
Wn   = [0.8 100] / (decode.Fs/2);
b = fir1(300, Wn,'bandpass');
a = 1;
pos_err_filt = filtfilt(b,a,pos_err);
pos_err_phase = angle(hilbert(pos_err_filt));
g1 = interp1(decode.pos,1:length(decode.pos), pos_real, 'nearest','extrap');
g2 = interp1(decode.pos,1:length(decode.pos), pos_estimate, 'nearest','extrap');
g = zeros(size(t));
g(diff(t)>2e6/decode.Fs) = 1;
g = cumsum(g)+1;

%%
figure
subplot(311)
hold on
plot(t,pos_estimate,'r.','MarkerSize',5);
plot(t,pos_real,'k.','MarkerSize',5);
rescale_plot_data('x',[1e-6 t(1)]);
subplot(312)
hold on
plot(t,pos_err,'.')
plot(t,pos_err_smooth,'r')
plot(t,pos_err_filt,'m')
ylim([-5 5])
rescale_plot_data('x',[1e-6 t(1)]);
subplot(313)
plot(ts,wingbeat)
rescale_plot_data('x',[1e-6 t(1)]);
ylim([-500 500])
linkaxes(findall(gcf,'type','axes'),'x');
zoom on

%%
x = interp1(ts,wingbeat_phase,t);
y = pos_err_phase;
figure
nbins = 200;
[N,C]=hist3([x;y]',[nbins nbins]);
N = N ./ mean(N,1);
%     contourf(C{1},C{2},N);
imagesc('CData',N,'XData',C{1},'YData',C{2});
axis equal
axis tight
title('wingbeat-pos_err phase-phase density plot')
xlabel('pos_diff Phase (rad)')
ylabel('wingbeat Phase (rad)')


%%
figure
hold on
x = decode.pos;
y = splitapply(@median, pos_err, g1);
err = splitapply(@nansem, pos_err, g1);
splitapply(@(x,y)(plot(x,y,'.')), pos_real,pos_err,g);
shadedErrorBar(x,y,err)
ylim([-5 5])
zoom on

%% pos_err autocorrelation
win_sec = 20;
win_samples = round(decode.Fs*win_sec);
[ccc,lags] = xcorr(pos_err,win_samples);
[ccc_smooth,lags] = xcorr(pos_err_smooth,win_samples);
[ccc_filt,lags] = xcorr(pos_err_filt,win_samples);
lags_sec = lags./decode.Fs;
figure
hold on
plot(lags_sec,ccc);
plot(lags_sec,ccc_smooth);
plot(lags_sec,ccc_filt);
legend("raw","smooth","filt")

%% error spectrum
figure
hold on
Nfft = 2^14;
pwelch(pos_err,Nfft,Nfft/2,Nfft,decode.Fs);
pwelch(pos_err_smooth,Nfft,Nfft/2,Nfft,decode.Fs);
pwelch(pos_err_filt,Nfft,Nfft/2,Nfft,decode.Fs);
h = findobj('Type','line');
h(1).Color = 'r';
h(2).Color = 'g';
h(3).Color = 'b';
xlim([0 20])
legend("raw","smooth","filt")

%% error spectrum (fit auto-regressive model)
figure
hold on
order = 4;
sys_b = ar(pos_err,order,'burg');
sys_fb = ar(pos_err,order);
spectrum(sys_b,sys_fb)
legend('Burg','Forward-Backward')

%% error spectrum (Autoregressive power spectral density estimate — Yule-Walker method)
y = pos_err;
y(abs(y)>5)=0;
order = 1000;
nfft = 2^12;
pyulear(y,order,nfft,decode.Fs);
xlim([0 20])

%% pos_
figure
hold on
plot(t,pos_err,'k')
plot(t,pos_err_filt,'b')
plot(t,pos_err_smooth,'g')
ylim(5*[-1 1])

%% pos err gradient hist
EDGES = linspace(-5,5,100);
histogram(diff(pos_err_smooth),EDGES)







%%
