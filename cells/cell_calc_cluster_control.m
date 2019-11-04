function cell_calc_cluster_control(cell_ID)

%%
% clear
% clc
% cell_ID = 'b0148_d170626_TT4_SS01';

%% load data
cell = cell_load_data(cell_ID,'details','fields','FE','spikes');
prm = PARAMS_GetAll();

%%
cluster_control = struct();
fields = [cell.fields{:}];
fields([fields.in_low_speed_area]) = [];
if length(fields)<2
    % not enough fields for comparison....
    cluster_control.pval_ranksum = nan;
    cluster_control.pval_ks = nan;
    cluster_control.field_size_vs_spike_size_r = nan;
    cluster_control.field_size_vs_spike_size_p = nan;

    % save data to file
    filename = fullfile('L:\Analysis\Results\cells\cluster_control',[cell_ID '_cell_cluster_control']);
    save(filename, 'cluster_control');
end

%% get spikes waveforms
[~,maxch] = max(max(mean(cell.spikes.waveforms,3),[],1));
% add spikes waveforms for field struct
for ii_field = 1:length(fields)
    field = fields(ii_field);
    IX = ismember(cell.spikes.ts,field.spikes_ts);
    fields(ii_field).spikes_wvfs = cell.spikes.waveforms(:,:,IX);
    fields(ii_field).spikes_wvfs_avg = mean(cell.spikes.waveforms(:,:,IX),3);
    fields(ii_field).spikes_wvfs_avg_max = mean(cell.spikes.waveforms(:,maxch,IX),3);
end

%% create in/out-field mask
nSpikes_per_field = cellfun(@length, {fields.spikes_ts});
mask_same_field = [];
for ii_field = 1:length(nSpikes_per_field)
    n = nSpikes_per_field(ii_field);
    field_mask = tril(ones(n),-1);
    mask_same_field = blkdiag(mask_same_field,field_mask);
end
mask_other_field = tril(ones(size(mask_same_field)),-1);
mask_other_field = mask_other_field - mask_same_field;
mask_same_field = boolean(mask_same_field);
mask_other_field = boolean(mask_other_field);

%% calc dist and stats
X = cat(3,fields.spikes_wvfs);
X = squeeze(X(:,maxch,:));
dist_metrics = {'euclidean';'correlation'};
spikes_dist_M = {};
for ii_dist = 1:length(dist_metrics)
    spikes_dist = pdist(X', dist_metrics{ii_dist});
    spikes_dist_M{ii_dist} = squareform(spikes_dist);
    % calc stats
    spikes_dist_same = spikes_dist_M{ii_dist}(mask_same_field);
    spikes_dist_other = spikes_dist_M{ii_dist}(mask_other_field);
    pval_ranksum(ii_dist) = ranksum(spikes_dist_same, spikes_dist_other);
    [~,pval_ks(ii_dist)] = kstest2(spikes_dist_same, spikes_dist_other);
end

%% plot a figure!
figure('Units','centimeters','Position',[5 5 30 15])
pnl = panel();
pnl.pack('h',[80 20]);
pnl(1).pack('h',2, 'v',[70 30]);
pnl(2).pack('v',{60 [] 45});
pnl.margin = [10 15 10 15];
pnl.de.margin = 10;
pnl(2).marginleft = 15;
delimiters = [0 cumsum(nSpikes_per_field)];
for ii_dist = 1:length(dist_metrics)
    
    pnl(1,ii_dist,1).select();
    hold on
    imagesc(spikes_dist_M{ii_dist});
    axis ij
    axis tight
    axis square
    xlimits = get(gca,'xlim');
    ylimits = get(gca,'ylim');
    for ii_line = 1:length(delimiters)
        plot(delimiters([ii_line ii_line]),ylimits,'r','LineWidth',.1);
        plot(xlimits,delimiters([ii_line ii_line]),'r','LineWidth',.1);
    end
    colorbar
    ha=gca;
    ha.XTick = [];
    ha.YTick = [];
    title(dist_metrics{ii_dist})

    pnl(1,ii_dist,2).select();
    hold on
    spikes_dist_same = spikes_dist_M{ii_dist}(mask_same_field);
    spikes_dist_other = spikes_dist_M{ii_dist}(mask_other_field);
    [y,x] = ksdensity(spikes_dist_same);
    plot(x,y,'k');
    [y,x] = ksdensity(spikes_dist_other);
    plot(x,y,'m');
    legend({'Same field';'Other field'},'Location','northeast')
    xlabel('Distance')
    ylabel('probability')
    text(0.3,0.98, sprintf('P(Ranksum)=%.2g',pval_ranksum(ii_dist)), 'Units','normalized');
    text(0.3,0.90, sprintf('P(KS)=%.2g',     pval_ks(ii_dist)),      'Units','normalized');
end

pnl(2,1).select();
plot([fields.spikes_wvfs_avg_max])
xlabel('Sample no.')
ylabel('Voltage (\muV)')

pnl(2,3).select();
x = [fields.width_prc];
y = range([fields.spikes_wvfs_avg_max]);
plot(x,y,'.','MarkerSize',15);
[r,p] = corr(x',y');
text(1,1,{sprintf('r=%.2g',r);sprintf('P=%.2g',p)},...
    'units','normalized','HorizontalAlignment','right','VerticalAlignment','top');
xlabel('Field size (m)')
ylabel('Spike Voltage (\muV)')

h=pnl.title(sprintf('%s (%d)',cell.details.cell_ID,cell.details.cell_num));
h.Interpreter='none';
h.FontSize = 13;
h.Position = [0.5 1.06];

% save figure
dir_out = 'L:\Analysis\Results\cells\cluster_control';
filename = sprintf('%s(%d)_cell_cluster_control', cell.details.cell_ID, cell.details.cell_num);
filename = fullfile(dir_out, filename);
saveas(gcf, filename, 'tif');

%% add data to struct
cluster_control.pval_ranksum = pval_ranksum;
cluster_control.pval_ks = pval_ks;
cluster_control.field_size_vs_spike_size_r = r;
cluster_control.field_size_vs_spike_size_p = p;

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\cluster_control',[cell_ID '_cell_cluster_control']);
save(filename, 'cluster_control');
    













end