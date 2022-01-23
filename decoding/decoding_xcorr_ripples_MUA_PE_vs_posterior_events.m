function decoding_xcorr_ripples_MUA_PE_vs_posterior_events(decode)
%%
exp_ID = decode.exp_ID;
epoch_type = decode.epoch_type;
params_opt = decode.params_opt;

%% IN/OUT folders
dir_OUT = 'F:\sequences\xcorr_ripples_vs_posterior_events';
mkdir(dir_OUT);

%% load data
exp = exp_load_data(exp_ID,'details','path','ripples','MUA','PE');
[events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
seqs = [events.seq_model];

%%
if length(events)<=1
    fprintf('No enough events for %s %s dec param %d\n', exp_ID, epoch_type, params_opt);
    fig=figure;
    text(.5,.5,'Not enough events!')
    sgtitle(sprintf('%s_%s_opt_%d',exp_ID, epoch_type, params_opt),'interpreter','none');
    filename = fullfile(dir_OUT, sprintf('%s_%s_opt_%d',exp_ID, epoch_type, params_opt));
    saveas(fig, filename, 'jpg');
    return;
end

%% plot figure
win = .5;
features = {'score','compression','confidence_sparsity','distance','duration','start_pos'};

fig=figure;
fig.WindowState = 'maximized';
tiledlayout(3,1+length(features),'TileSpacing','tight');

% arrange ripple/MUA in one struct
signals = struct();
signals(1).fs = exp.ripples.fs;
signals(1).t = exp.ripples.t;
signals(1).val = exp.ripples.zpripple;
signals(1).events = exp.ripples.events;
signals(1).label = 'Ripple power (z)';
signals(1).label2 = "ripples";
signals(2).fs = exp.MUA.fs;
signals(2).t = exp.MUA.t;
signals(2).val = exp.MUA.zFR;
signals(2).events = exp.MUA.events;
signals(2).label = 'MUA firing rate (z)';
signals(2).label2 = "MUA";
signals(3).events = exp.PE.thr;
signals(3).label2 = "PE";

for ii_signal = 1:length(signals)
    
    %% get signal
    signal = signals(ii_signal);
    if isempty(signal.t)
        continue;
    end
    
    %% trigger signals (ripples/MUA) around posterior events
    win_samples = round(win*signal.fs);
    trigger_IX = interp1(signal.t, 1:length(signal.t), [events.peak_ts], 'nearest');
    [trig_signal] = trigger_signal_by_IX(signal.val, trigger_IX, win_samples);
    t = linspace(-win, win, size(trig_signal,2));

    %% plot posterior-event-triggered signals
    nexttile
    hold on
    t = linspace(-win, win, size(trig_signal,2));
    shadedErrorBar(t, trig_signal, {@nanmean,@nansem});
    xlabel('Time from posterior event (s)')
    ylabel(signal.label)

    %%  plot scatters of ripples/MUA vs posterior event features
    for ii_feature = 1:length(features)
        feature_fn = features{ii_feature};

        %% group by nearby ripple/MUA/PE/other
        t0 = [events.peak_ts];
        t1 = [exp.ripples.events.peak_ts];
        t2 = [exp.MUA.events.peak_ts];
        t3 = [exp.PE.thr.peak_ts];
        IX1 = any(abs((t0-t1').*1e-6)<win);
        IX2 = any(abs((t0-t2').*1e-6)<win);
        IX3 = any(abs((t0-t3').*1e-6)<win);
        g = zeros(size(events));
        g(IX1) = 1;
        g(IX2) = 2;
        g(IX3) = 3;
        g = g+1;
        clrs = [0 0 0; 0 0 1; 0 1 0; 1 0 0];

        %% plot values at posterior event peak
%         nexttile
%         hold on
%         center_IX = win_samples+1;
%         x = [seqs.(feature_fn)];
%         y = trig_signal(:,center_IX);
%         scatter(x,y,5,clrs(g,:),'filled');
%         lsline
%         xlabel(feature_fn,'Interpreter','none');
%         ylabel(signal.label+ " @ time=0");

        %% plot values at posterior event window
        nexttile
        hold on
        x = [seqs.(feature_fn)];
        y = max(trig_signal,[],2);
        scatter(x,y,5,clrs(g,:),'filled');
        lsline
        xlabel(feature_fn,'Interpreter','none');
        ylabel(signal.label+ " @ max in win");

        %% create legend
        h=splitapply(@(c)(plot(nan,nan,'.','MarkerSize',15,'Color',c)),clrs',1:size(clrs,1));
        h=legend(h,"other","near ripple","near MUA","near PE",'Location','bestoutside');
        h.Position = [0.88 0.94 .05 .05];
    end
end

% plot events xcorr
for ii_signal = 1:length(signals)
    signal = signals(ii_signal);
    nexttile
    t1 = [events.peak_ts];
    t2 = [signal.events.peak_ts];
    tdiff = (t1-t2').*1e-6;
    histogram(tdiff(:), 'BinLimits',[-1 1].*win, 'BinWidth',.050);
    text(.95,.95,sprintf('%.2g%%',100*sum(any(abs(tdiff)<win))/size(tdiff,2)),'Units','normalized','HorizontalAlignment','right');
    title('posterior vs. ' + signal.label2);
    xlabel('Time from posterior event (s)');
    ylabel(signal.label2+' count')
end

% plot movemevnt probability triggered by ripples/MUA/PE
mvmnt_states_IX = contains(decode.state,'movement');
mvmnt_prob = decode.posterior_state(mvmnt_states_IX,:);
mvmnt_prob = max(mvmnt_prob,[],1);
% mvmnt_prob = sum(mvmnt_prob,1);
win = 2;
win_samples = round(win*decode.Fs);
for ii_signal = 1:length(signals)
    signal = signals(ii_signal);
    nexttile
    hold on
    ts = [signal.events.peak_ts];
    ib = find_nearest_point(ts,decode.time);
    valid = abs(ts - decode.time(ib)) < win*1e6/decode.Fs;
    ib(~valid)=[];
    [trig_signal] = trigger_signal_by_IX(mvmnt_prob, ib, win_samples);
    t = linspace(-win, win, size(trig_signal,2));
    shadedErrorBar(t, trig_signal, {@nanmean,@nansem});
    xlabel(sprintf('Time from %s (s)',signal.label2))
    ylabel('Movement prob.')
    text(.95,.95,"n="+length(ib),'Units','normalized','HorizontalAlignment','right');
end

sgtitle(sprintf('%s_%s_opt_%d',exp_ID, epoch_type, params_opt),'interpreter','none');
filename = fullfile(dir_OUT, sprintf('%s_%s_opt_%d',exp_ID, epoch_type, params_opt));
saveas(fig, filename, 'jpg');


%%


