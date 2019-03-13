function PRE_calc_write_artifacts(exp_ID)

%%
% clear
% clc
% exp_ID = 'b9861_d180607';

%% load exp details
exp=exp_load_data(exp_ID,'details','path');
main_dir = exp.path.spikes_raw;
switch exp.details.NeuralLoggerType
    case 'SpikeLog-16'
        block_size = 2^8;
    case 'MiniBatLog-16'
        block_size = 2^11;
end
dir_OUT = 'L:\Analysis\Results\pre_proc\SD_write_artifacts';
mkdir(dir_OUT)

%% run over TT/channels and calc mean/median/std/...
tic
for TT = 1:size(exp.details.activeChannels,1)
    for ch = 1:size(exp.details.activeChannels,1)
        if ~exp.details.activeChannels(TT,ch)
            continue
        end
        %% load raw data
        filename = fullfile(main_dir, sprintf('spikes_%s_TT%d_ch%d.ncs',exp_ID,TT,ch));
        [signal, ts, fs, params] = Nlx_csc_read(filename, []);
        %% arrange by write blocks and calc
        blocks = reshape(signal, block_size, []);
        blocks_mean{TT,ch} = mean(blocks,2);
        blocks_abs_mean{TT,ch} = mean(abs(blocks),2);
        blocks_median{TT,ch} = median(blocks,2);
        blocks_abs_median{TT,ch} = median(abs(blocks),2);
        blocks_std{TT,ch} = std(blocks,0,2);
    end
end
toc

%% plot figure - abs mean
figure('Units','normalized','Position',[0 0 1 1])
pnl = panel();
pnl.pack(size(exp.details.activeChannels,1), size(exp.details.activeChannels,2));
h=pnl.title({'calculating the mean of absolute raw data over SD writing blocks';exp_ID});
h.Interpreter = 'none';
h.Position = [0.5 0.98];
pnl.margin = [20 30 20 20];
pnl.de.margin = 10;
for TT = 1:size(exp.details.activeChannels,1)
    for ch = 1:size(exp.details.activeChannels,1)
        if ~exp.details.activeChannels(TT,ch)
            continue
        end
        pnl(TT,ch).select()
        x = blocks_abs_mean{TT,ch};
        t = linspace(0,length(x)/fs*1e3, length(x));
        plot(t,x);
        ylim([0 1.1*max(x)]);
        xlabel('Time (msec)')
        ylabel('voltage (uV)')
    end
end
file_OUT = fullfile(dir_OUT, sprintf('%s_SD_write_artifacts_abs_mean',exp_ID));
saveas(gcf, file_OUT, 'tif')

%% plot figure - abs median
figure('Units','normalized','Position',[0 0 1 1])
pnl = panel();
pnl.pack(size(exp.details.activeChannels,1), size(exp.details.activeChannels,2));
h=pnl.title({'calculating the median of absolute raw data over SD writing blocks';exp_ID});
h.Interpreter = 'none';
h.Position = [0.5 0.98];
pnl.margin = [20 30 20 20];
pnl.de.margin = 10;
for TT = 1:size(exp.details.activeChannels,1)
    for ch = 1:size(exp.details.activeChannels,1)
        if ~exp.details.activeChannels(TT,ch)
            continue
        end
        pnl(TT,ch).select()
        x = blocks_abs_median{TT,ch};
        t = linspace(0,length(x)/fs*1e3, length(x));
        plot(t,x);
        ylim([0 1.1*max(x)]);
        xlabel('Time (msec)')
        ylabel('voltage (uV)')
    end
end
file_OUT = fullfile(dir_OUT, sprintf('%s_SD_write_artifacts_abs_median',exp_ID));
saveas(gcf, file_OUT, 'tif')


%% artifact/noise ratio
artifact_noise_ratio_abs_mean = nan(size(exp.details.activeChannels));
artifact_noise_ratio_abs_median = nan(size(exp.details.activeChannels));
for TT = 1:size(exp.details.activeChannels,1)
    for ch = 1:size(exp.details.activeChannels,1)
        if ~exp.details.activeChannels(TT,ch)
            continue
        end
        artifact_noise_ratio_abs_mean(TT,ch) = max(blocks_abs_mean{TT,ch}) / min(blocks_abs_mean{TT,ch});
        artifact_noise_ratio_abs_median(TT,ch) = max(blocks_abs_median{TT,ch}) / min(blocks_abs_median{TT,ch});
    end
end

figure('Units','normalized','Position',[0 0 1 1])
pnl = panel();
pnl.pack('h',2);
pnl.margin = 50;
h=pnl.title({'SD writing artifact vs. noise level ration'; exp_ID});
h.Interpreter = 'none';
h.Position = [0.5 1.05];
h.FontSize = 14;

pnl(1).select();
title('abs mean');
M = artifact_noise_ratio_abs_mean;
h=imagesc(M);
axis ij
set(h,'AlphaData',~isnan(M));
ylabel('tetrode')
xlabel('channel')
set(gca,'XTick', 1:size(exp.details.activeChannels,1), 'YTick', 1:size(exp.details.activeChannels,2));
h=colorbar;
set(get(h,'title'),'String','ratio');
[x y] = meshgrid(1:size(exp.details.activeChannels,1), 1:size(exp.details.activeChannels,2));
M_str = arrayfun(@(x)(sprintf('%.1f',x)),M,'UniformOutput',0);
text(x(:),y(:),{M_str{:}},'color','r','FontSize',12,'FontWeight','Bold')

pnl(2).select()
title('abs median');
M = artifact_noise_ratio_abs_median;
h=imagesc(M);
axis ij
set(h,'AlphaData',~isnan(M));
ylabel('tetrode')
xlabel('channel')
set(gca,'XTick', 1:size(exp.details.activeChannels,1), 'YTick', 1:size(exp.details.activeChannels,2));
h=colorbar;
set(get(h,'title'),'String','ratio');
[x y] = meshgrid(1:size(exp.details.activeChannels,1), 1:size(exp.details.activeChannels,2));
M_str = arrayfun(@(x)(sprintf('%.1f',x)),M,'UniformOutput',0);
text(x(:),y(:),{M_str{:}},'color','r','FontSize',12,'FontWeight','Bold')

file_OUT = fullfile(dir_OUT, sprintf('%s_artifacts_noise_ratio',exp_ID));
saveas(gcf, file_OUT, 'jpg')

%%
close all






