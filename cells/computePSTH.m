function [PSTH,spike_density_filtered,time_spent_filtered] = computePSTH(pos,pos_fs,spikes_pos,bin_edges,min_time_spent,ker_SD)
% 1 dimensional inputs: bsp_pos   1xN
%                       ts        1xN
%                       spike_pos  1xM
%                       bin_edges 1xL

%% default params
if nargin == 4 %so no input of minimum time spent and kernel s.d.
    min_time_spent = default_min_time_spent;
    ker_SD = 0.5;
end

%% binning
bin_size = (bin_edges(end) - bin_edges(1))/(length(bin_edges)-1);
bin_centers = (bin_edges(1)+bin_size/2):bin_size:bin_edges(end);
[spike_density, ~] = density_mat_1D(spikes_pos, bin_centers);
[pos_counts,    ~] = density_mat_1D(pos,        bin_centers);
time_spent = pos_counts*(1/pos_fs);

%% apply minimum time spent
mask = time_spent>=min_time_spent;
if sum(mask) <= 1
    PSTH = nan(size(mask));
    spike_density_filtered = nan(size(mask));
    time_spent_filtered = nan(size(mask));
    return
end
% take only the main continuous component of valid bins
CC = bwconncomp(mask);
[~,max_IX] = max(cellfun(@length, CC.PixelIdxList));
mask(:) = 0;
mask(CC.PixelIdxList{max_IX}) = 1;

%% create filter kernel
hsize = 1 + (5*ker_SD/bin_size);% %ceil(5*(param.general.sigma_PF_smoothing*3));
hsize = round(hsize);
alpha = hsize*bin_size/(2*ker_SD); %the equation is:
% alpha = hsize/2*sigma  where sigma is in number of samples - so the sigma_PF_smoothing
% (m) is devided by bin size (m) in order to get the sigma in number of bins unit.
ker = gausswin(hsize,alpha)'./(sqrt(2*pi)*ker_SD);%fspecial('gaussian',hsize,param.general.sigma_PF_smoothing);

%% filter and then divide
time_spent_filtered     = nanfilt(time_spent,   mask,ker);
spike_density_filtered  = nanfilt(spike_density,mask,ker);
PSTH = spike_density_filtered ./ time_spent_filtered ;

end

%%
function xfilt = nanfilt(x,mask,ker)

IX = 1:numel(x);
x(~mask) = interp1(IX(mask), x(mask), IX(~mask), 'nearest', 'extrap');
xfilt = imfilter(x,ker,'same','conv','symmetric');
xfilt(~mask) = nan;

end

