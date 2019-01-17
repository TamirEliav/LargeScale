function [L_Ratio,Isolation_dis,features_space] = cell_calc_cluster_quality(SpikesWaveforms,ind)

%%
num_spikes_total = size(SpikesWaveforms,3);
num_spikes_cluster = length(ind);
num_spikes_noise = num_spikes_total - num_spikes_cluster;

%% get feature space
SpikesWaveforms = permute(SpikesWaveforms,[3,1,2]); %Mx32x4
[nSpikes, nSamp, nch] = size(SpikesWaveforms);

nrg = sqrt(sum(SpikesWaveforms.^2,2))/nSamp; %energy per spikes per channel
tmp = permute(nrg,[1,3,2]);
inactive_ch = find(sum(abs(tmp),1)==0);
nrg(:,:,inactive_ch) = [];
SpikesWaveforms(:,:,inactive_ch) = [];

Nrg_factor = repmat(nrg,1,32,1);
all_data_normalized = SpikesWaveforms./(Nrg_factor+eps); %eps is added in 
%order to deal with the case of a non existing channel

nrg = reshape(nrg,nSpikes,4-length(inactive_ch)); %energy per spikes per channel Mx4
for i = 1:(4-length(inactive_ch))
    [coeffC,scores,latent] = pca(all_data_normalized(:,:,i));
    pc1(:,i) = scores(:,1);
end

features_space = [nrg pc1]; %Mx8
noise_feature_vect = features_space; noise_feature_vect(ind,:) = [];
cluster_feature_vect = features_space(ind,:);

%% Mahalanobis distance:
D_2_noise = mahal(noise_feature_vect,cluster_feature_vect);

%% L-Ratio:
CDF_noiseSpikes = cdf('Chisquare',D_2_noise,8-2*(length(inactive_ch))); % 8 degrees of freedom
L_Ratio = sum(ones([num_spikes_noise,1]) - CDF_noiseSpikes)/num_spikes_cluster;

%% Isolation Distance
D_2_noise_sorted = sort(D_2_noise,'ascend');
if num_spikes_cluster > length(D_2_noise_sorted)
    Isolation_dis = D_2_noise_sorted(end);
else
    Isolation_dis = D_2_noise_sorted(num_spikes_cluster);
end



