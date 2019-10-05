function cell_calc_cluster_control(cell_ID)

%%
clear
clc
cell_ID = 'b0148_d170626_TT4_SS01';

%% load data
cell = cell_load_data(cell_ID,'details','fields','FE','spikes');
prm = PARAMS_GetAll();

%%
cluster_control = struct();
fields = [cell.fields{:}];
fields([fields.in_low_speed_area]) = [];
% if length(fields)<2
%     % not enough fields for comparison - TODO: do something...
% else
%     
% end

[~,maxch] = max(max(mean(cell.spikes.waveforms,3),[],1));
% add spikes waveforms for field struct
for ii_field = 1:length(fields)
    field = fields(ii_field);
    IX = ismember(cell.spikes.ts,field.spikes_ts);
    fields(ii_field).spikes_wvfs = cell.spikes.waveforms(:,:,IX);
    fields(ii_field).spikes_wvfs_avg = mean(cell.spikes.waveforms(:,:,IX),3);
    fields(ii_field).spikes_wvfs_avg_max = mean(cell.spikes.waveforms(:,maxch,IX),3);
end
% 
for ii_field = 1:length(fields)
    in_field_IX = ismember(1:length(fields),ii_field);
    in_field_wvfrm  = cat(3,fields( in_field_IX).spikes_wvfs);
    out_field_wvfrm = cat(3,fields(~in_field_IX).spikes_wvfs);
    % get only the maximal channel waveform
    in_field_wvfrm = squeeze(in_field_wvfrm(:,maxch,:));
    out_field_wvfrm = squeeze(out_field_wvfrm(:,maxch,:));
    sdf = pdist2(in_field_wvfrm', out_field_wvfrm');
end

%%
X = cat(3,fields.spikes_wvfs);
X = squeeze(X(:,maxch,:));
sdf = pdist(X');
nSpikes_per_field = cellfun(@length, {fields.spikes_ts});
delimiters = [0 cumsum(nSpikes_per_field)];
figure
hold on
imagesc(squareform(sdf));
xlimits = get(gca,'xlim');
ylimits = get(gca,'ylim');
for ii_line = 1:length(delimiters)
    plot(delimiters([ii_line ii_line]),ylimits,'r','LineWidth',.1);
    plot(xlimits,delimiters([ii_line ii_line]),'r','LineWidth',.1);
end
axis ij
colorbar

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\cluster_control',[cell_ID '_cell_cluster_control']);
save(filename, 'cluster_control');
    













end