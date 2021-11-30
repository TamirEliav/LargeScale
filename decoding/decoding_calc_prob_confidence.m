function [KS, HPD, sparsity] = decoding_calc_prob_confidence(prob,HPD_prc)
arguments
    prob
    HPD_prc = 0.95
end
%%
HPD_prc = 0.95;
nX = size(prob,1); % num position bins
nT = size(prob,2); % num time bins
prob2 = sort(prob,1,'descend');
prob3 = cumsum(prob2,1);
HPD = zeros(1,nT);
KS = zeros(1,nT);
sparsity = zeros(1,nT);
test_cdf = cdf('Discrete Uniform',1:nX,nX);
for ii_ts = 1:nT
%     HPD(ii_ts) = find(prob3(:,ii_ts)>HPD_prc,1,'first'); % this approach
%     doesn't work when the total probability is less than the 95%, which
%     is very often in our case (because we look at the probability at the
%     movement state only and not summed over all states!
    HPD(ii_ts) = nX - sum(prob3(:,ii_ts)>HPD_prc);
%     KS(ii_ts) = max(abs(prob3(:,ii_ts)'-test_cdf)); % this can be 
    sparsity(ii_ts) = cumputeSparsity(prob3(:,ii_ts));
end
KS = max(abs(prob3'-test_cdf),[],2)'; % use this!

% normalize
HPD = interp1([1 nX], [1 0], HPD, 'linear');
KS_min = max(abs(test_cdf-test_cdf)); % simply zero...
KS_max = max(abs(ones(1,nX)-test_cdf));
KS = interp1(linspace(KS_min,KS_max,nX), linspace(0,1,nX), KS, 'linear');
sparsity = sparsity .* sum(prob,1); % weight according to the probability to be in the state

% % average across time bins
% HPD = mean(HPD);
% KS = mean(KS);
% sparsity = mean(sparsity);
 
% fig=figure;
% fig.WindowState = 'maximized';
% subplot(311)
% hold on
% plot(HPD)
% plot(KS)
% plot(sparsity)
% legend(sprintf('HPD@%d%%',HPD_prc*100),'KS','sparsity','Location','northoutside')
% xlabel('Time (bins)')
% ylabel('decoding confidence')
% subplot(312)
% plot(sum(prob))
% xlabel('Time (bins)')
% ylabel('State prob.')
% subplot(313)
% imagesc(prob)
% linkaxes(findall(gcf, 'type', 'axes'),'x')
% xlabel('Time (bins)')
% ylabel('Position (bins)')
% figname = sprintf('decoding_confidence_measures_event_example_HPD_%d',HPD_prc*100);
% h=suptitle(figname);
% h.Interpreter = 'none';
% % dir_OUT = 'Z:\USE_THIS_DIRECTORY_FOR_FILE_TRANSFER\for_Nachum_from_Tamir\20211124__replay_quantifying_seq';
% % saveas(fig,fullfile(dir_OUT,figname),'jpg');
% % subplot(311)
% % ylim([0.9 1])
% % saveas(fig,fullfile(dir_OUT,[figname '_ylim_zoom']),'jpg');

end