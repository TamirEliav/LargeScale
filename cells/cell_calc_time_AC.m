function cell_calc_time_AC(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','FE','fields');
% exp = exp_load_data(cell.details.exp_ID,'details');
prm = PARAMS_GetAll();

%% in-flight 
FE = [cell.FE{:}];
spikes_ts = [FE.spikes_ts];
AC = cumpute_time_AC(spikes_ts*1e-6, ...
                     prm.oscillations.time_AC.bin_size,...
                     prm.oscillations.time_AC.win_size);
time_AC.in_flight = AC;

%% by flight direction
for ii_dir = 1:2
    FE = cell.FE{ii_dir};
    spikes_ts = [FE.spikes_ts];
    AC = cumpute_time_AC(spikes_ts*1e-6, ...
                         prm.oscillations.time_AC.bin_size,...
                         prm.oscillations.time_AC.win_size);
    time_AC.by_dir(ii_dir) = AC;
end

%% TODO: fix cases with fields
%% in-field
fields = [cell.fields{:}];
spikes_ts = [fields.spikes_ts];
AC = cumpute_time_AC(spikes_ts*1e-6, ...
                     prm.oscillations.time_AC.bin_size,...
                     prm.oscillations.time_AC.win_size);
time_AC.in_field = AC;

%% by fields
fields = [cell.fields{:}];
for ii_field = 1:length(fields)
    spikes_ts = [fields(ii_field).spikes_ts];
    AC = cumpute_time_AC(spikes_ts*1e-6, ...
                         prm.oscillations.time_AC.bin_size,...
                         prm.oscillations.time_AC.win_size);
    time_AC.by_fields(ii_field) = AC;
end

%% sleep

%% resting on balls

%% all recorded spikes

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\time_AC',[cell_ID '_cell_time_AC']);
save(filename, 'time_AC');


end


function AC = cumpute_time_AC(ts,bin_size,win_size)

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












