function [] = NLG_check_flash_artifacts(main_dir, time_start_end, logger_num)

%% 
% main_dir = 'D:\Tamir\PROJECTS\Neurologger\testing\data\test_2016_06_10__flash_write_artifcats_different_loggers';
% logger_num = 14;

dir_out = fullfile(main_dir, 'artifacts');
mkdir(dir_out);

%%
data = [];
num_channels = 16;
for ii_ch = 1:num_channels
    disp(['Reading ch ' num2str(ii_ch) ' out of ' num2str(num_channels)])
    file_in = fullfile(main_dir, 'nlx', ['CSC' num2str(ii_ch-1) '.ncs']);
    [signal, ts, fs] = Nlx_csc_read(file_in, time_start_end);
    last_IX = find(abs(signal)<8000,1, 'last');
    data(ii_ch, :) = signal(257:last_IX); % remove start/end rec noise/invlaid data
end
data_reshaped = reshape(data, 16, 256, []);
data_mean = squeeze(mean(data_reshaped,3));
data_mean2 = data_mean - repmat(mean(data_mean(:,150:end),2),1,256);
data_std = squeeze(std(data_reshaped,1,3));

%% plot
figure
imagesc(data_mean2)
caxis( [ min(data_mean2(:)) max(data_mean2(:)) ] )
colorbar
colormap jet
title(['noise mean - logger #' num2str(logger_num)])
xlabel('Sample#')
ylabel('Channel#')
file_out = fullfile(dir_out, ['noise_mean_colorplot_logger_' num2str(logger_num)]);
saveas(gcf, file_out,'fig')
saveas(gcf, file_out,'jpeg')

figure
imagesc(data_std)
caxis( [ min(data_std(:)) max(data_std(:)) ] )
colorbar
colormap jet
title(['noise std - logger #' num2str(logger_num)])
xlabel('Sample#')
ylabel('Channel#')
file_out = fullfile(dir_out, ['noise_std_colorplot_logger_' num2str(logger_num)]);
saveas(gcf, file_out,'fig')
saveas(gcf, file_out,'jpeg')

figure
ax_h = [];
for ii_ch=1:16
    ax_h(ii_ch) = subaxis(4,4,ii_ch);
    plot( data_mean2(ii_ch,:) ) 
    ylim([min(data_mean2(:)) max(data_mean2(:))])
    xlim([1 256])
    title(['ch ' num2str(ii_ch)])
end
xlabel('sample#')
ylabel('Voltage (uVolt)')
suptitle(['noise mean - logger #' num2str(logger_num)])
file_out = fullfile(dir_out, ['noise_mean_logger_' num2str(logger_num)]);
saveas(gcf, file_out,'fig')
saveas(gcf, file_out,'jpeg')

figure
ax_h = [];
for ii_ch=1:16
    ax_h(ii_ch) = subaxis(4,4,ii_ch);
    plot( data_std(ii_ch,:) ) 
    ylim([0 max(data_std(:))])
    xlim([1 256])
    title(['ch ' num2str(ii_ch)])
end
xlabel('sample#')
ylabel('Voltage (uVolt)')
suptitle(['noise std - logger #' num2str(logger_num)])
file_out = fullfile(dir_out, ['noise_std_logger_' num2str(logger_num)]);
saveas(gcf, file_out,'fig')
saveas(gcf, file_out,'jpeg')

%% power spectrum
figure
WINDOW = 2^12;
NOVERLAP = 0;
NFFT = 2^12;
for ii_ch=1:16
    ax_h(ii_ch) = subaxis(4,4,ii_ch);
    pwelch(data(ii_ch,:), WINDOW, NOVERLAP, NFFT, fs)
    title(['ch ' num2str(ii_ch)])
end
suptitle(['power spectrum - logger #' num2str(logger_num)])
file_out = fullfile(dir_out, ['power_specturm_' num2str(logger_num)]);
saveas(gcf, file_out,'fig')
saveas(gcf, file_out,'jpeg')

% noise std
figure
ch_std = [];
for ii_ch=1:16
    ch_std(ii_ch) = std(squeeze(data(ii_ch,:)));
end
bar(ch_std)
title(['noise levels - logger #' num2str(logger_num)])
file_out = fullfile(dir_out, ['noise_levels_' num2str(logger_num)]);
saveas(gcf, file_out,'fig')
saveas(gcf, file_out,'jpeg')




%%