%% test transfter function of nlg

%%
clear 
clc

%%
channels = 0:15;
% channels = 0;
MPP_thr = 100;
num_freq_clusters = 30;
for ch = channels

    %% load data
    logger_SN = 5;
    dir_IN = 'L:\DATA\2289_Sami\test\20180613__test_transfer_function';
    file_IN = fullfile(dir_IN, ['logger_' num2str(logger_SN)' '\nlx\CSC' num2str(ch) '.ncs']);
    [signal, ts, fs] = Nlx_csc_read(file_IN, []);

    %% detect peaks and anlayze individual cycles
    
    [~,peaks_IX] = findpeaks(signal, 'MinPeakProminence',MPP_thr);

    freqs = 1e6./diff(ts(peaks_IX));
    subs = zeros(size(ts));
    subs(peaks_IX) = 1;
    subs = cumsum(subs);
    valid_IX = find(subs~=0);
    amps = accumarray(subs(valid_IX)', signal(valid_IX)', [], @peak2peak);
    amps = amps .* sqrt(2)/4; % convert P2P to rms (as in the minirator)
    amps(end) = [];

    %%
    % Y = discretize(freq,edges);
    [idx,C,sumd,D] = kmeans(freqs', num_freq_clusters);
    % f_amps = accumarray(idx,amps,[],@max);
    f_amps = accumarray(idx,amps,[],@(x)(prctile(x,95)));
    f = accumarray(idx,freqs,[],@median);
    
    %% plot all cycles
    figure
    plot(freqs,amps, '.')
%     dscatter(freqs',amps)
    xlabel('Freq (Hz)')
    ylabel('Amplitude (uV)')
    ylim([0 1200])

    %% plot by freq clusters (and save fig)
    figure
    stem(f, f_amps, 'o')
    set(gca,'xscale','log')
    xlabel('Freq (Hz)')
    ylabel('Amplitude (uV)')
    ylim([0 1200])
    title(['Transfter function - logger #' num2str(logger_SN) ' channel ' num2str(ch)]);
    file_out = fullfile(dir_IN, ['Transfter_function_logger_' num2str(logger_SN) '_channel_' num2str(ch)]);
    saveas(gcf, file_out , 'tif')

end % loop over channels