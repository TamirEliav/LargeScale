function cell_calc_time_AC(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','FE','fields','spikes');
exp = exp_load_data(cell.details.exp_ID,'details');
prm = PARAMS_GetAll();

%%
calc_AC = @(ts)(cumpute_time_AC(ts*1e-6, ...
                         prm.oscillations.time_AC.bin_size,...
                         prm.oscillations.time_AC.win_size));

%% in-flight 
FE = [cell.FE{:}];
spikes_ts = [FE.spikes_ts];
time_AC.in_flight = calc_AC(spikes_ts);

%% by flight direction
for ii_dir = 1:2
    FE = cell.FE{ii_dir};
    spikes_ts = [FE.spikes_ts];
    time_AC.by_dir(ii_dir) = calc_AC(spikes_ts);
end

%% by fields
fields = [];
if ~isempty(cell.fields{1})
    fields = [fields cell.fields{1}];
end
if ~isempty(cell.fields{2})
    fields = [fields cell.fields{2}];
end
if length(fields)==0
    time_AC.in_field = [];
    time_AC.by_fields = [];
else
    %% in-field
    spikes_ts = [fields.spikes_ts];
    time_AC.in_field = calc_AC(spikes_ts);

    %% by fields
    for ii_field = 1:length(fields)
        spikes_ts = [fields(ii_field).spikes_ts];
        time_AC.by_fields(ii_field) = calc_AC(spikes_ts);
    end
end

%% sleep
sleep_ti = exp_get_sessions_ti(exp.details.exp_ID, 'Sleep1', 'Sleep2');
spikes_sleep_IX = get_data_in_ti(cell.spikes.ts, sleep_ti);
spikes_sleep_ts = cell.spikes.ts(spikes_sleep_IX);
time_AC.sleep = calc_AC(spikes_sleep_ts);

%% TODO: add resting on balls

%% all recorded spikes

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\time_AC',[cell_ID '_cell_time_AC']);
save(filename, 'time_AC');


end


function AC = cumpute_time_AC(ts,bin_size,win_size)

if isempty(ts)
    AC.lags = [];
    AC.c = [];
    AC.bin_size = bin_size;
    AC.win_size = win_size;
    return;
end

edges = ts(1) : bin_size : ts(end);
[N,edges] = histcounts(ts,edges);
maxlag = round(win_size/bin_size);
[c, lags] = xcorr(N, maxlag);
lags = lags .* bin_size;
t0_IX = lags==0;
c(t0_IX) = nan;
c(t0_IX) = max(c);

% c = smooth(c);

AC.lags = lags;
AC.c = c;
AC.bin_size = bin_size;
AC.win_size = win_size;

end












