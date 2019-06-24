function wingbeat_detect(exp_ID)

%% load exp data
exp = exp_load_data(exp_ID,'details','path','flight','pos');
prm = PARAMS_GetAll();

%% read flight rhythm
files = dir(fullfile(exp.path.LFP, ['LFP_' exp_ID '*.ncs']));
file_IN = fullfile(files(1).folder, files(1).name);
[signal_raw, ts, fs] = Nlx_csc_read(file_IN, []);

%% filter the signal
passband = [3 10];
order = 300;
Wn   = passband / (fs/2);
[b,a] = fir1(order,Wn,'bandpass');
signal = filtfilt(b,a,signal_raw);

%% take signal only during flight
flight_IX = get_data_in_ti(ts,[[exp.flight.FE.start_ts]; [exp.flight.FE.end_ts]]');
signal(~ismember(1:length(ts), flight_IX)) = nan;

%% find artifact peaks = wingbeats (probably proxy by respiration)
m = std(signal(flight_IX));
[PKS,LOCS,width,prom] = findpeaks(signal...
    ,'MinPeakProminence', 0.01*m...
    ,'MinPeakDistance',0.080*fs...
    ,'WidthReference', 'halfprom'...
    ,'MinPeakWidth', 0.01*fs...
    ,'MaxPeakWidth', 0.5*fs...
    );

%% plot detection
figure('Units','normalized','Position',[0 0 1 1]);
hold on
behave_ts = exp_get_sessions_ti(exp_ID,'Behave');
t0 = behave_ts(1);
t = (ts-t0)*1e-6/60;
plot(t,signal,'k')
plot(t(LOCS),signal(LOCS),'r*')
xlabel('Time (minutes)', 'FontSize', 14)
ylabel('Voltage (uV)', 'FontSize', 14)
title('Wingbeat detection', 'FontSize', 14)
filename = fullfile('L:\Analysis\Results\exp\wingbeat',[exp_ID '_exp_wingbeat_detect']);
saveas(gcf, filename,'tif')
saveas(gcf, filename,'fig')
close(gcf)

%% create out data structure
wingbeat = struct();
wingbeat.ts = ts(LOCS);
wingbeat.pos = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos,       wingbeat.ts, 'linear');
wingbeat.vel = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.vel_csaps, wingbeat.ts, 'linear');

%% save data to file
filename = fullfile('L:\Analysis\Results\exp\wingbeat',[exp_ID '_exp_wingbeat']);
save(filename, 'wingbeat');

end

