function cell_calc_cluster_quality2(cell_ID)

%%
% clear
% clc
% cell_ID = 'b0148_d170626_TT4_SS01';

%% load data
cell = cell_load_data(cell_ID,'details');
exp=exp_load_data(cell.details.exp_ID, 'path','details');
prm = PARAMS_GetAll();
time_segment_duration = 10*60*1e6; % un usec

%% load cell spikes data from NTT file of the relevant TT
NTT_file = fullfile(exp.path.spikes_sorting,...
                    ['spikes_' cell.details.TT_day_ID '.NTT']);
if isempty(cell.details.stable_ts)
    [Timestamps, CellNumbers, Samples, Header] = ...
         Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 1, [] );
else
    [Timestamps, CellNumbers, Samples, Header] = ...
         Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 4, cell.details.stable_ts );
end

%% parse header (convert bits to uVolts)
ADBitVolts = sscanf(Header{contains(Header,'ADBitVolts')},'-ADBitVolts %f %f %f %f');
Samples = Samples .* ADBitVolts' .* 1e6;

%% create time segments
if isempty(cell.details.stable_ts)
    time_segments = Timestamps([1 end]);
else
    time_segments = cell.details.stable_ts;
end
% time_segments = time_segments(1):time_segment_duration:time_segments(end);
time_segments = linspace(time_segments(1), time_segments(end), ceil(diff(time_segments)/time_segment_duration)+1);
time_segments = time_segments([1:end-1;2:end])';

%% calc cluster qualirt measures in time segments
L_Ratio=[];
Isolation_dis=[];
for ii_seg = 1:size(time_segments,1)
    % get all samples in the time segment
    Samples_seg_IX = Timestamps>time_segments(ii_seg,1) & Timestamps<time_segments(ii_seg,2);
    Samples_seg = Samples(:,:,Samples_seg_IX);
    CellNumbers_seg = CellNumbers(Samples_seg_IX);
    % get the IX of the time segment samples which are from our unit
    IX = find(  CellNumbers_seg == cell.details.unit );
    % calc quality measures
    if length(IX) > 8 % #samples > #features
        [L_Ratio(ii_seg), Isolation_dis(ii_seg), ~] = ...
            cell_calc_cluster_quality(Samples_seg, IX);
    else
        L_Ratio(ii_seg) = nan;
        Isolation_dis(ii_seg) = nan;
    end
end

%% add data to struct
cluster_quality = struct();
cluster_quality.Isolation_dis = Isolation_dis;
cluster_quality.L_Ratio = L_Ratio;
cluster_quality.Isolation_dis_mean = nanmean(Isolation_dis);
cluster_quality.Isolation_dis_median = nanmedian(Isolation_dis);
cluster_quality.L_Ratio_mean = nanmean(L_Ratio);
cluster_quality.L_Ratio_median = nanmedian(L_Ratio);
cluster_quality.time_segments = time_segments;
cluster_quality.n_time_segments = size(time_segments,1);
cluster_quality.time_segment_duration = time_segment_duration;

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\cluster_quality',[cell.details.cell_ID '_cell_cluster_quality']);
save(filename, 'cluster_quality');
    













end