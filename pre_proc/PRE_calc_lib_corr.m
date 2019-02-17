function rrr = PRE_calc_lib_corr(wvfrms, lib)
% wvfrms - 32x4xN  (N spikes)
% lib    - Mx32    (M patterns)

% For each event take the channel with the largest peak (8th point)
    [~,max_ch_IX] = max(squeeze(abs(wvfrms(8,:,:))),[],1);
    largest_waveforms = zeros(size(wvfrms,1),size(wvfrms,3));
    for ch=1:4
        IX = find(max_ch_IX == ch);
        largest_waveforms(:,IX) = squeeze(wvfrms(:,ch,IX));
    end
    % now, let's calc the max corr
    xxx_lags_shifts = [1:30; 2:31; 3:32];
    ccc = [];
    rrr = [];
    for ii_shift = 1:size(xxx_lags_shifts,1)
        xxx_lags = xxx_lags_shifts(ii_shift, :);
        ccc = corr(largest_waveforms(xxx_lags,:), lib(:,2:end-1)');
        rrr(ii_shift,:) = max(ccc,[],2);
    end
    rrr = max(rrr,[],1);

end