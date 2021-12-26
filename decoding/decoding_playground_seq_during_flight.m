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
pos_diff = pos_estimate-pos_real;
pos_diff_smooth = smoothdata(pos_diff,'gaussian',0.1*decode.Fs);
Wn   = [4 10] / (decode.Fs/2);
b = fir1(300, Wn,'bandpass');
a = 1;
pos_diff_filt = filtfilt(b,a,pos_diff);
pos_diff_phase = angle(hilbert(pos_diff_filt));
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
plot(t,pos_diff,'.')
plot(t,pos_diff_smooth,'r')
plot(t,pos_diff_filt,'m')
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
y = pos_diff_phase;
figure
nbins = 200;
[N,C]=hist3([x;y]',[nbins nbins]);
N = N ./ mean(N,1);
%     contourf(C{1},C{2},N);
imagesc('CData',N,'XData',C{1},'YData',C{2});
axis equal
axis tight
title('wingbeat-pos_diff phase-phase density plot')
xlabel('pos_diff Phase (rad)')
ylabel('wingbeat Phase (rad)')


%%
figure
hold on
x = decode.pos;
y = splitapply(@median, pos_diff, g1);
err = splitapply(@nansem, pos_diff, g1);
splitapply(@(x,y)(plot(x,y,'.')), pos_real,pos_diff,g);
shadedErrorBar(x,y,err)
ylim([-5 5])
zoom on

%%








%%
