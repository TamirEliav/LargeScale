%%
clear
clc

%% load exp data
exp_ID = 'b0184_d191130';
dir_IN = 'F:\sequences\proc';
filename = fullfile(dir_IN,exp_ID);
data = load(filename);
% exp = exp_load_data(exp_ID, 'details','MUA','flight','rest');
exp = exp_load_data(exp_ID, 'details','flight','rest');

%% load LFP
[LFP, LFP_ts, LFP_fs, LFP_params, ch_valid] = LFP_load(exp_ID);

%% calc MUA with higher spike thr
thr = 500;
% exp.MUA.FR = sum(squeeze(max(data.multiunits,[],2))>thr,2)';
% sdf = mean(data.multiunits_spikes,2);
sdf = mean(squeeze(max(data.multiunits,[],2))>thr,2);
exp.MUA.FR = data.fs .*  smoothdata(sdf, 'gaussian', 37.5)';
exp.MUA.t = data.t;
exp.MUA.fs = data.fs;

%% filter LFP / MUA 
passband = [3 10];
forder = 1000;

b = fir1(forder,passband./exp.MUA.fs*2);
exp.MUA.FR_filtered = filtfilt(b,1,exp.MUA.FR);
exp.MUA.phase = angle(hilbert(exp.MUA.FR_filtered));
exp.MUA.FR_MA = smoothdata(exp.MUA.FR,'movmedian',exp.MUA.fs*5);

b = fir1(forder,passband./LFP_fs*2);
noise = mean(LFP,[2 3]);
noise_filtered = filtfilt(b,1,noise);
noise_phase = angle(hilbert(noise_filtered));

%% get flight indexes
FE = exp.flight.FE;
FE_ti = [FE.start_ts; FE.end_ts]';
sleep_ti = exp_get_sessions_ti(exp.details.exp_ID, 'Sleep1','Sleep2');
in_flight_IX = get_data_in_ti(exp.MUA.t, FE_ti);
in_rest_IX = get_data_in_ti(exp.MUA.t, exp.rest.ti);
in_sleep_IX = get_data_in_ti(exp.MUA.t,  sleep_ti);
epochs_IX = {in_flight_IX, in_rest_IX, in_sleep_IX };
epochs_name = {'flight','rest','sleep'};

%% plot FR for the entire session
figure
hold on
plot(exp.MUA.FR)
plot(exp.MUA.FR_filtered + exp.MUA.FR_MA)
% plot(exp.MUA.t, exp.MUA.FR)
% plot(exp.MUA.t, exp.MUA.FR_filtered + exp.MUA.FR_MA)
% plot(exp.MUA.t, exp.MUA.phase.*20)
% plot(ts,noise)
% plot(ts,noise_filtered)
ylabel('Firing rate (Hz)')
% ylabel('noise (uV)')
rescale_plot_data('x',[1/exp.MUA.fs 0]);
% rescale_plot_data('x',[1e-6/60 FE_ti(1)]);
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend('MUA','MUA filtered+slow moving average')

%% AC(MUA) + spectrum + xcorr(MUA,LFP) + phase-phase(LFP-MUA)
fig=figure
fig.WindowState = 'maximized';
ht = tiledlayout(length(epochs_name),4,'TileSpacing','compact');
for ii_epoch = 1:length(epochs_IX)
    IX = epochs_IX{ii_epoch};
    epoch_name = epochs_name{ii_epoch};

    FR = exp.MUA.FR(IX);
    Fs = exp.MUA.fs;
    
    nexttile
    [c, lags] = xcorr(FR);
    lags_s = lags ./ exp.MUA.fs;
    plot(lags_s,c);
    xlim([-1 1].*0.5)    
    hax=gca;
    hax.YLim(1) = 0;
    xlabel('Time lag (s)')
    ylabel('Counts')
    title('MUA firing rate autocorrelation')
    text(-.25,0.5,epoch_name,'Units','normalized','Rotation',90,'FontSize',16,'FontWeight','bold','HorizontalAlignment','center');

    nexttile
    n=15;
    pwelch(FR,2^n,2^(n-1),2^n,Fs);
    xlim([0 20])
    
    nexttile
    x = exp.MUA.FR;
    y = interp1(LFP_ts, noise_phase, exp.MUA.t, 'linear','extrap');
    IX = epochs_IX{ii_epoch};
    x = x(IX);
    y = y(IX);
    [c, lags] = xcorr(x,y,'normalized');
    plot(lags,c);
    rescale_plot_data('x', [1/exp.MUA.fs 0])
    xlim([-1 1].*1)
    xlabel('Time lag (s)')
    title('MUA-LFP cross-correlation')
    
    nexttile
    x = exp.MUA.phase;
    y = interp1(LFP_ts, noise_phase, exp.MUA.t, 'linear','extrap');
    IX = epochs_IX{ii_epoch};
    x = x(IX);
    y = y(IX);
    nbins = 200;
    [N,C]=hist3([x;y]',[nbins nbins]);
    N = N ./ mean(N,1);
%     contourf(C{1},C{2},N);
    imagesc('CData',N,'XData',C{1},'YData',C{2});
    axis equal
    axis tight
    title('MUA-LFP phase-phase density plot')
    xlabel('LFP Phase (rad)')
    ylabel('MUA Phase (rad)')
end
h=suptitle(exp.details.exp_ID);
h.Interpreter = 'none';
h.FontSize = 18;
h.FontWeight = 'bold';
h.Units='normalized';
h.Position(2) = -0.025;


%%









%%
