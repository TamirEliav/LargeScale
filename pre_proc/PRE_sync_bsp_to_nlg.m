function PRE_sync_bsp_to_nlg(bsp_dir, nlg_dir, out_dir)
% TODO: make sure to plot figures showing good sync!

%%
% eval(expname); 

%%
% bsp_dir =  fullfile(main_dir, 'bsp');
% nlg_dir =  fullfile(main_dir, 'nlx');
% out_dir = fullfile(main_dir, 'sync');
% bsp_dir =  param.path.bsp;
% nlg_dir =  param.path.Nlx;
% out_dir = param.path.sync;
mkdir(out_dir)

%% BSP TTL
load( fullfile(bsp_dir, 'bsp_TTL.mat') );
bsp_TTL_ts_msec = round(1e-6.*bsp_TTL_ts_ns');
bsp_TTL_intervals = diff(bsp_TTL_ts_msec);
bsp_TTL_intervals_inc = diff(bsp_TTL_intervals);

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
ax_h(1) = subplot(2,2,1);
plot(bsp_TTL_ts_msec(2:end)*1e-3/60,bsp_TTL_intervals*1e-3,'.')
% xlabel('Time (minutes)')
ylabel('inter-TTL-interval (sec)')
title('bsp TTL timings')
ax_h(3) = subplot(2,2,3);
plot(bsp_TTL_ts_msec(3:end)*1e-3/60,bsp_TTL_intervals_inc,'.')
xlabel('Time (minutes)')
ylabel('inter-TTL-interval increament (msec)')
ax_h(2) = subplot(2,2,2);
plot(nlg_TTL_ts_msec(2:end)*1e-3/60,nlg_TTL_intervals*1e-3,'.')
title('nlg TTL timings')
% xlabel('Time (minutes)')
ylabel('inter-TTL-interval (sec)')
ax_h(4) = subplot(2,2,4);
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
saveas(gcf, fullfile(out_dir, 'sync_bsp_nlg__TTL_timinigs'), 'jpeg')
saveas(gcf, fullfile(out_dir, 'sync_bsp_nlg__TTL_timinigs'), 'fig')

%% sync
time_conv_p_msec = sync_TTL_polyfit(round(bsp_TTL_ts_msec), round(nlg_TTL_ts_msec), 2, 100);
saveas(gcf, fullfile(out_dir, 'sync_bsp_nlg'), 'jpeg')

%% TODO: save polyfit values (if we want to convert again...)
save( fullfile(out_dir, 'sync_time_conv_p') , 'time_conv_p_msec')

%% convert timestamps
files_to_convert = dir(fullfile(bsp_dir, '*pos*.mat'))
for ii_file = 1:length(files_to_convert)
    clear bsp_pos
    load( fullfile(bsp_dir, files_to_convert(ii_file).name) );
    bsp_pos.ts_nlg_usec = polyval(time_conv_p_msec,bsp_pos.ts_ns*1e-6)*1e3;
    save( fullfile(bsp_dir, files_to_convert(ii_file).name), 'bsp_pos');
end

%%






