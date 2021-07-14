function ripples_trigger_LFP(exp_ID)

%% get exp info
exp = exp_load_data(exp_ID,'details','path','ripples');
prm = PARAMS_GetAll();

%% load raw LFP data
[LFP, ts, fs, params, ch_valid] = LFP_load(exp_ID);
nTT = size(LFP,2);
nCh = size(LFP,3);

%% average LFP over channels per TT
LFP = squeeze(nanmean(LFP,3));

%% create ripples triggered LFP figure
win_sec = 0.5;
h=figure;
set(gcf,'DefaultAxesFontSize',6);
h.Units = 'centimeters';
h.Position = [3 3 20 20]; 
pnl=panel();
pnl.pack(nTT+1,nTT);
pnl.margin=[40 25 20 30];
pnl.de.margin=10;
pnl.fontsize=6;
h=pnl.xlabel('Triggered signal');
h.Position(2)=-0.12;
h.FontSize = 12;
h=pnl.ylabel('Detection by');
h.Position(1)=-0.18;
h.FontSize = 14;
h=pnl.title({exp_ID;'ripples triggered LFP'});
h.Interpreter = 'none';
h.FontSize = 16;
h.Position(2)=1.1;
for TT_detect = 1:nTT
    ripples = exp.ripples.by_TT{TT_detect};
    if isempty(ripples)
        continue;
    end
    for TT_trigger = 1:nTT
        %%
        traces = trigger_signal_by_IX(LFP(:,TT_trigger), [ripples.peak_IX], round(win_sec*fs));
        traces_mean = mean(traces,1);
        traces_SEM = std(traces,0,1)/sqrt(size(traces,1));
        traces_t = linspace(-win_sec,win_sec, size(traces,2)).*1e3;
        
        %%
        pnl(TT_detect,TT_trigger).select();
        shadedErrorBar(traces_t,traces_mean,traces_SEM)
        xline(0)
        text(1,1,"n="+length(ripples),'Units','normalized','FontSize',6,'HorizontalAlignment','right');
        if TT_trigger == 1
            ylabel('Voltage (uV)')
            text(-0.65,0.5,"TT "+TT_detect,'Units','normalized','FontSize',12,'HorizontalAlignment','center');
        end
    end
end

%% trigger using detection from all tetrodes (and plot...)
for TT_trigger = 1:nTT
    ripples = exp.ripples.all;
    traces = trigger_signal_by_IX(LFP(:,TT_trigger), [ripples.peak_IX], round(win_sec*fs));
    traces_mean = mean(traces,1);
    traces_SEM = std(traces,0,1)/sqrt(size(traces,1));
    traces_t = linspace(-win_sec,win_sec, size(traces,2)).*1e3;
    pnl(nTT+1,TT_trigger).select();
    shadedErrorBar(traces_t,traces_mean,traces_SEM)
    xlabel('Time (ms)')
    xline(0)
    text(1,1,"n="+length(ripples),'Units','normalized','FontSize',6,'HorizontalAlignment','right');
    if TT_trigger==1
        text(-0.5,0.5,"All TT ",'Units','normalized','FontSize',12,'HorizontalAlignment','center');
    end
    text(0.5,-0.5,"TT "+TT_trigger,'Units','normalized','FontSize',12,'HorizontalAlignment','center');
end

%% save figure
fig_filename = fullfile('L:\Analysis\Results\exp\ripples', [exp_ID '_ripples_triggered_LFP']);
saveas(gcf,fig_filename,'jpeg')

end

        %%
        % traces_power = ripples.peak_zpripple(ripples.valid_ripples);
        % [~,sort_IX] = sort(traces_power,'descend');
        % traces = traces(sort_IX,:);
        % traces_power = traces_power(sort_IX);
        % for ii = 1:size(traces,1)
        %     plot(traces_t,traces(ii,:));
        %     title("traces # "+ii+" z="+traces_power(ii));
        %     pause
        % end

        %%
        % figure
        % violinplot(ripples.peak_zpripple, ripples.valid_ripples)
        % violinplot(ripples.peak_ripple_gamma_ratio, ripples.valid_ripples)
        % violinplot(ripples.ripple_gamma_power_ratio, ripples.valid_ripples)
        % violinplot(ripples.durations, ripples.valid_ripples)


        %%
%         nfft=2^8;
%         win = [nfft nfft/2] ./ fs;
%         win = [0.01 0.001];
%         clear params
%         params.Fs = fs;
%         params.trialave = 1;
%         params.err = 0;
%         params.fpass = [50 250];
%         [S,t,f] = mtspecgramc( traces', win, params);
%         figure
%         imagesc(t-t(end)/2,f,log(S)');
%         axis xy
%         colormap jet
%         imagesc(t,f, log(squeeze(mean(S,3))));
%         [S,t,f,Serr] = mtspecgramc( traces', movingwin, params );

        %%  
%         nfft = 2^8;
%         [s,f,t] = stft(traces',fs,'Window',hamming(nfft,'periodic'),'OverlapLength',nfft-1,'FFTLength',nfft);
%         imagesc(t,f,log(squeeze(mean(abs(s),3))))



