function PRE_sync_audio_to_nlg(audio_dir, nlg_dir, out_dir)
% TODO: make sure to plot figures showing good sync!

%% Audio TTL
audio_TTL_file_name = fullfile(audio_dir, 'EVENTS__Digital in.nev');
FieldSelection = [1 0 0 0 0];
ExtractHeader = 0;
ExtractMode = 1;
ModeArray = [];
audio_TTL_ts_usec = Nlx2MatEV( audio_TTL_file_name ,FieldSelection,ExtractHeader,ExtractMode,ModeArray);
audio_TTL_ts_msec = audio_TTL_ts_usec*1e-3;
audio_TTL_intervals = diff(audio_TTL_ts_msec);
audio_TTL_intervals_inc = diff(audio_TTL_intervals);

%% NLG TTL
nlg_TTL_file_name = fullfile(nlg_dir, 'EVENTS__Digital in.nev');
FieldSelection = [1 0 0 0 0];
ExtractHeader = 0;
ExtractMode = 1;
ModeArray = [];
nlg_TTL_ts_usec = Nlx2MatEV( nlg_TTL_file_name ,FieldSelection,ExtractHeader,ExtractMode,ModeArray);
nlg_TTL_ts_msec = nlg_TTL_ts_usec*1e-3;
nlg_TTL_intervals = diff(nlg_TTL_ts_msec);
nlg_TTL_intervals_inc = diff(nlg_TTL_intervals);

%% plot TTL timings
ax_h = [];
figure
ax_h(1) = subaxis(2,2,1);
plot(audio_TTL_ts_msec(2:end)*1e-3/60,audio_TTL_intervals*1e-3,'.')
% xlabel('Time (minutes)')
ylabel('inter-TTL-interval (sec)')
title('audio TTL timings')
ax_h(3) = subaxis(2,2,3);
plot(audio_TTL_ts_msec(3:end)*1e-3/60,audio_TTL_intervals_inc,'.')
xlabel('Time (minutes)')
ylabel('inter-TTL-interval increament (msec)')
ax_h(2) = subaxis(2,2,2);
plot(nlg_TTL_ts_msec(2:end)*1e-3/60,nlg_TTL_intervals*1e-3,'.')
title('nlg TTL timings')
% xlabel('Time (minutes)')
ylabel('inter-TTL-interval (sec)')
ax_h(4) = subaxis(2,2,4);
plot(nlg_TTL_ts_msec(3:end)*1e-3/60,nlg_TTL_intervals_inc,'.')
xlabel('Time (minutes)')
ylabel('inter-TTL-interval increament (msec)')
ylimits = get(ax_h([1 2]), 'ylim');
ylimits = [  min([ylimits{:}]) max([ylimits{:}])  ]
linkaxes(ax_h([1 2]),'y')
set(ax_h(1), 'ylim' , ylimits);
ylimits = get(ax_h([3 4]), 'ylim');
ylimits = [  min([ylimits{:}]) max([ylimits{:}])  ]
linkaxes(ax_h([3 4]),'y')
set(ax_h(3), 'ylim' , ylimits);
saveas(gcf, fullfile(out_dir, 'sync_audio_nlg__TTL_timinigs'), 'jpeg')
saveas(gcf, fullfile(out_dir, 'sync_audio_nlg__TTL_timinigs'), 'fig')

%% sync
time_conv_p_msec = sync_TTL_polyfit(round(audio_TTL_ts_msec), round(nlg_TTL_ts_msec), 2, 100);
saveas(gcf, fullfile(out_dir, 'sync_audio_nlg'), 'jpeg')

%% TODO: save polyfit values (if we want to convert again...)
save( fullfile(out_dir, 'sync_time_conv_p') , 'time_conv_p_msec')

%% convert timestamps
file_in = fullfile(audio_dir, 'audio.ncs');
file_out = fullfile(audio_dir, 'audio_nlg_ts.ncs');
[signal, ts_audio_usec, fs] = Nlx_csc_read(file_in, []);
ts_nlg_usec = polyval(time_conv_p_msec, ts_audio_usec*1e-3)*1e3;
nlx_csc_write(file_out, signal, ts_nlg_usec , fs);

%%






