function [PSTH,spike_density_filtered,time_spent_filtered] = computePSTH(pos,pos_fs,spikes_pos,bin_edges,min_time_spent)
% 1 dimensional inputs: bsp_pos   1xN
%                       ts        1xN
%                       spike_pos  1xM
%                       bin_edges 1xL

default_min_time_spent = 0.15; %sec
if nargin == 4 %so no input of minimum time spent
    min_time_spent = default_min_time_spent;
end

bin_size = (bin_edges(end) - bin_edges(1))/(length(bin_edges)-1);
bin_centers = (bin_edges(1)+bin_size/2):bin_size:bin_edges(end);

% ts(isnan(ts)) = []; bsp_pos(isnan(bsp_pos)) = [];
% 
% spike_pos = interp1(ts, bsp_pos, spike_ts,'linear');

[spike_density, spike_bin] = density_mat_1D(spikes_pos,bin_centers);
[bsp_samples_spent, ~] = density_mat_1D(pos,bin_centers);

time_spent = bsp_samples_spent*(1/pos_fs);
%mask = time_spent~=0;
mask = time_spent>=min_time_spent; mask(1) = false; mask(end) = false;
% Isolated bins:
check = sum(abs(diff(mask)));
if check>2 % check if there is more than one block. takes the largest continous block of behavior.
    ind_low = find(diff(mask)<0);
    ind_high = find(diff(mask)>0);
    blocks_size = ind_low - ind_high;
    [max_block_size,ind_block] = max(blocks_size);
    if max_block_size*bin_size>100 %m
        mask(1:ind_high(ind_block)) = false;
        mask(ind_low(ind_block)+1:end) = false;
    end
end

PSTH = spike_density ./ time_spent;

%%
sigma_PF_smoothing = 0.5;
hsize = 1 + (5*sigma_PF_smoothing/bin_size);% %ceil(5*(param.general.sigma_PF_smoothing*3));
hsize = round(hsize);
alpha = hsize*bin_size/(2*sigma_PF_smoothing); %the equation is:
% alpha = hsize/2*sigma  where sigma is in number of samples - so the sigma_PF_smoothing
% (m) is devided by bin size (m) in order to get the sigma in number of bins unit.

gaussian_kernel = gausswin(hsize,alpha)./(sqrt(2*pi)*sigma_PF_smoothing);%fspecial('gaussian',hsize,param.general.sigma_PF_smoothing);

timespent_with_NaN_not_filtered = time_spent;
PSTH_not_filtered = spike_density ./ timespent_with_NaN_not_filtered;

time_spent_filtered(~mask) = 0;
spike_density_filtered(~mask) = nan;
spike_density_filtered(mask) = conv(spike_density(mask),gaussian_kernel,'same');
time_spent_filtered(mask)    = conv(time_spent(mask),gaussian_kernel,'same');

PSTH = spike_density_filtered ./ time_spent_filtered ;
% PSTH(mask) = conv(PSTH(mask),gaussian_kernel,'same');

end



