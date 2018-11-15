function TTL_ts_usec = NLX_sync_TTL_extract(nlx_TTL_file)


%% read CSC file
[signal ts Fs] = Nlx_csc_read(nlx_TTL_file,[]);
Ts = 1/Fs;

%% notch filter (50 Hz)
use_50Hz_notch = 0;
if use_50Hz_notch
    Q = 5;
    Fo = 50;
    BW = (Fo/(Fs/2))/Q;
    N = round(Fs/Fo);
    [b,a] = iircomb(N,BW);  
    filtered_signal = filtfilt(b,a,signal);
else
    filtered_signal= signal;
end

%% detect
pulse_width_msec = 100;
pulse_width_in_samples = round(Fs*pulse_width_msec*1e-3);
pulse_rising_edge_pos = pulse_width_in_samples/2;
pulse_template = [zeros(1,round(pulse_width_in_samples/2)) ones(1,round(pulse_width_in_samples/2))];

template_match_score = conv(filtered_signal, fliplr(pulse_template), 'same');
template_match_score_th = mean(template_match_score) + 3*std(template_match_score);
[PKS,TTL_IX] = findpeaks(template_match_score, 'MINPEAKHEIGHT', template_match_score_th, 'MINPEAKDISTANCE', pulse_width_in_samples)

TTL_ts_usec = ts(TTL_IX);

%% DBG - plotting detection process
if 1
%%
figure
hold on
plot(ts, template_match_score )
plot(ts([1 end]), [template_match_score_th template_match_score_th ],'--m')
plot(ts(TTL_IX),PKS,'*g')

%%
figure
hold on
plot(ts, filtered_signal)
plot(ts(TTL_IX), filtered_signal(TTL_IX),'*g')

%% trigger signal around detected TTL rising edge to make sure it works fine
trigger_time_win_msec = 200;
trigger_time_win_in_samples = round(Fs*trigger_time_win_msec*1e-3);
signal_TTL_triggered = [];
for ii_TTL = 2:length(TTL_IX)-1
    
    IX_start = TTL_IX(ii_TTL) - trigger_time_win_in_samples;
    IX_end = TTL_IX(ii_TTL) + trigger_time_win_in_samples;
    IX = IX_start:IX_end;
    signal_TTL_triggered(ii_TTL,:) = filtered_signal(IX);
    
end
figure; hold on;
plot(signal_TTL_triggered')
plot( [trigger_time_win_in_samples trigger_time_win_in_samples], [min(signal_TTL_triggered(:)) max(signal_TTL_triggered(:))], '--g', 'LineWidth', 5)

end










