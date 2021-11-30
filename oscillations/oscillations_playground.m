%% load exp data
exp_ID = 'b0184_d191130';
exp = exp_load_data(exp_ID, 'details','MUA','flight','rest');

%% load LFP
[LFP, LFP_ts, LFP_fs, LFP_params, ch_valid] = LFP_load(exp_ID);

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

%%
figure
yyaxis left
plot(exp.MUA.t, exp.MUA.FR_filtered)
ylabel('Firing rate (Hz)')
yyaxis right
plot(ts,noise_filtered)
ylabel('noise (uV)')
xlabel('minutes')


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
    xlabel('Time lag (s)')
    ylabel('Counts')
    title('MUA firing rate autocorrelation')
    text(-.25,0.5,epoch_name,'Units','normalized','Rotation',90,'FontSize',16,'FontWeight','bold','HorizontalAlignment','center');

    nexttile
    n=15;
    pwelch(FR,2^n,2^(n-1),2^n,Fs);
    xlim([0 20])
    
    x = exp.MUA.phase;
    y = interp1(LFP_ts, noise_phase, exp.MUA.t, 'linear');
    IX = epochs_IX{ii_epoch};
    x = x(IX);
    y = y(IX);

    nexttile
    [c, lags] = xcorr(x,y,'normalized');
    plot(lags,c);
    rescale_plot_data('x', [1/exp.MUA.fs 0])
    xlim([-1 1].*1)
    xlabel('Time lag (s)')
    title('MUA-LFP cross-correlation')
    
    nexttile
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
IX = in_flight_IX;
% IX = in_rest_IX;
FR = exp.MUA.FR(IX);
Fs = exp.MUA.fs;
[c, lags] = xcorr(FR);
lags_s = lags ./ Fs;
figure
hold on
plot(lags_s,c);
xlim([-1 1].*0.5)

n=15;
figure
pwelch(FR,2^n,2^(n-1),2^n,Fs)
xlim([0 20])

%%
dir_OUT = 'D:\';
file_name = fullfile(dir_OUT, [exp_ID '_MUA_FR' '.ncs']);
nlx_csc_write(file_name, exp.MUA.FR, exp.MUA.t, exp.MUA.fs);
file_name = fullfile(dir_OUT, [exp_ID '_MUA_FR_filtered' '.ncs']);
nlx_csc_write(file_name, exp.MUA.FR_filtered, exp.MUA.t, exp.MUA.fs);


%%
cells_data = load('Z:\students_personal_backup_space\Shir\cells184_forYuval\cells_for_yuval_bat184.mat');
cells = cells_data.cells184;

%%
details = [cells.details];
cells(~contains({details.exp_ID},exp_ID)) = [];

%% phase-position plot
ii_cell = 2;
cell = cells(ii_cell);
figure
hold on
hax=gca;
hax.ColorOrder([1:4],:) = hax.ColorOrder([1 1 2 2],:);
for ii_dir = 1:2
%     subplot(2,1,ii_dir);
    spikes_pos = [cell.FE{ii_dir}.spikes_pos];
    spikes_ts  = [cell.FE{ii_dir}.spikes_ts];
    spikes_phase = interp1(exp.MUA.t, exp.MUA.phase, spikes_ts, 'linear');
    plot(spikes_pos, spikes_phase, '.')
    plot(spikes_pos, spikes_phase+2*pi, '.')
end
ylim([-pi 3*pi])
yline(pi)
title(cell.details.cell_ID)

%% phase AC
ii_cell = 9;
nBins = 30;
cell = cells(ii_cell);
all_fields = [cell.fields.all.light{:}];
spikes_ts  = [all_fields.spikes_ts];
spikes_phase = interp1(exp.MUA.t, unwrap(exp.MUA.phase), spikes_ts, 'linear');
spikes_phase = spikes_phase ./ (2*pi);
h=histogram(spikes_phase-spikes_phase','BinWidth',1/nBins, 'BinLimits',[-1 1].*4);
xlim([-1 1].*4)






%%




%%
