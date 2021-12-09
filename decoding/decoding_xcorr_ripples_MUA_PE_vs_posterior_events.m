function decoding_xcorr_ripples_MUA_PE_vs_posterior_events(exp_ID, epoch_type, params_opt)
arguments
    %% 
    exp_ID = 'b9861_d180526'
    epoch_type {mustBeMember(epoch_type,{'sleep','rest','flight'})} = 'sleep'
    params_opt = 11;
%     event_type {mustBeMember(event_type,{'PE','posterior','ripples','MUA'})} = 'posterior'
end

%% IN/OUT folders
dir_OUT = 'F:\sequences\xcorr_ripples_vs_posterior_events';
mkdir(dir_OUT);

%% load data
exp = exp_load_data(exp_ID,'details','path','ripples','MUA','PE');
decode = decoding_load_data(exp_ID, epoch_type, params_opt);
[events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
seqs = [events.seq_model];

%%
if length(events)<=1
    fprintf('No enough events for %s %s dec param %d %s\n', exp_ID, epoch_type, params_opt);
    fig=figure;
    text(.5,.5,'Not enough events!')
    sgtitle(sprintf('%s_%s_opt_%d',exp_ID, epoch_type, params_opt),'interpreter','none');
    filename = fullfile(dir_OUT, sprintf('%s_%s_opt_%d',exp_ID, epoch_type, params_opt));
    saveas(fig, filename, 'jpg');
    return;
end

%%
win = .5;
sort_features = {'score','compression','confidence_sparsity'};

fig=figure;
fig.WindowState = 'maximized';
tiledlayout(3,1+2*length(sort_features));

% ripples
win_samples = round(win*exp.ripples.fs);
trigger_IX = interp1(exp.ripples.t, 1:length(exp.ripples.t), [events.peak_ts], 'nearest');
[trig_signal] = trigger_signal_by_IX(exp.ripples.zpripple_all, trigger_IX, win_samples);
nexttile
hold on
t = linspace(-win, win, size(trig_signal,2));
shadedErrorBar(t, trig_signal, {@nanmean,@nansem});
xlabel('Time from posterior event (s)')
ylabel('Ripples power (z)')
for ii_feature = 1:length(sort_features)
    nexttile
    feature_fn = sort_features{ii_feature};
    [~,sort_IX] = sort([seqs.(feature_fn)],'descend');
    imagesc(trig_signal(sort_IX,:),'XData',t)
    title(feature_fn,'Interpreter','none')
    ylabel(sprintf('events (sorted by %s)',feature_fn),'Interpreter','none')
    xlabel('Time from event (s)')
    nexttile
    center_IX = win_samples+1;
    x = [seqs.(feature_fn)];
    y = trig_signal(:,center_IX);
    plot(x,y,'.');
    lsline
    xlabel(feature_fn,'Interpreter','none')
    ylabel('Ripples power (z)')
end

% MUA
win_samples = round(win*exp.MUA.fs);
trigger_IX = interp1(exp.MUA.t, 1:length(exp.MUA.t), [events.peak_ts], 'nearest');
[trig_signal] = trigger_signal_by_IX(exp.MUA.zFR, trigger_IX, win_samples);
nexttile
hold on
t = linspace(-win, win, size(trig_signal,2));
shadedErrorBar(t, trig_signal, {@nanmean,@nansem});
xlabel('Time from posterior event (s)')
ylabel('MUA firing rate (z)')
for ii_feature = 1:length(sort_features)
    nexttile
    feature_fn = sort_features{ii_feature};
    [~,sort_IX] = sort([seqs.(feature_fn)],'descend');
    imagesc(trig_signal(sort_IX,:),'XData',t)
    title(feature_fn,'Interpreter','none')
    ylabel(sprintf('events (sorted by %s)',feature_fn),'Interpreter','none')
    xlabel('Time from event (s)')
    nexttile
    center_IX = win_samples+1;
    x = [seqs.(feature_fn)];
    y = trig_signal(:,center_IX);
    plot(x,y,'.');
    lsline
    xlabel(feature_fn,'Interpreter','none')
    ylabel('MUA firing rate (z)')
end

nexttile
x1 = [events.peak_ts];
x2 = [exp.ripples.all.peak_ts];
xdiff = (x1-x2').*1e-6;
histogram(xdiff(:), 'BinLimits',[-1 1].*win, 'BinWidth',.050);
text(.95,.95,sprintf('%.2g%%',100*sum(any(abs(xdiff)<win))/size(xdiff,2)),'Units','normalized','HorizontalAlignment','right');
title('posterior vs. ripples')
xlabel('Time lag (s)')
ylabel('Counts')

nexttile
x1 = [events.peak_ts];
x2 = [exp.MUA.events.peak_ts];
xdiff = (x1-x2').*1e-6;
histogram(xdiff(:), 'BinLimits',[-1 1].*win, 'BinWidth',.050);
text(.95,.95,sprintf('%.2g%%',100*sum(any(abs(xdiff)<win))/size(xdiff,2)),'Units','normalized','HorizontalAlignment','right');
title('posterior vs. MUA')
xlabel('Time lag (s)')
ylabel('Counts')

nexttile
x1 = [events.peak_ts];
x2 = [exp.PE.thr.peak_ts];
xdiff = (x1-x2').*1e-6;
histogram(xdiff(:), 'BinLimits',[-1 1].*win, 'BinWidth',.050);
text(.95,.95,sprintf('%.2g%%',100*sum(any(abs(xdiff)<win))/size(xdiff,2)),'Units','normalized','HorizontalAlignment','right');
title('posterior vs. PE')
xlabel('Time lag (s)')
ylabel('Counts')

mvmnt_states_IX = contains(decode.state,'movement');
mvmnt_prob = decode.posterior_state(mvmnt_states_IX,:);
mvmnt_prob = max(mvmnt_prob,[],1);
% mvmnt_prob = sum(mvmnt_prob,1);
win = 2;
win_samples = round(win*decode.Fs);

nexttile
hold on
ts = [exp.ripples.all.peak_ts];
ib = find_nearest_point(ts,decode.time);
valid = abs(ts - decode.time(ib)) < 2e6/decode.Fs;
ib(~valid)=[];
[trig_signal] = trigger_signal_by_IX(mvmnt_prob, ib, win_samples);
t = linspace(-win, win, size(trig_signal,2));
shadedErrorBar(t, trig_signal, {@nanmean,@nansem});
xlabel('Time from ripple (s)')
ylabel('Movement prob.')
text(.95,.95,"n="+length(ib),'Units','normalized','HorizontalAlignment','right');

nexttile
hold on
ts = [exp.MUA.events.peak_ts];
ib = find_nearest_point(ts,decode.time);
valid = abs(ts - decode.time(ib)) < 2e6/decode.Fs;
ib(~valid)=[];
[trig_signal] = trigger_signal_by_IX(mvmnt_prob, ib, win_samples);
t = linspace(-win, win, size(trig_signal,2));
shadedErrorBar(t, trig_signal, {@nanmean,@nansem});
xlabel('Time from MUA (s)')
ylabel('Movement prob.')
text(.95,.95,"n="+length(ib),'Units','normalized','HorizontalAlignment','right');

nexttile
hold on
ts = [exp.PE.thr.peak_ts];
ib = find_nearest_point(ts,decode.time);
valid = abs(ts - decode.time(ib)) < 2e6/decode.Fs;
ib(~valid)=[];
[trig_signal] = trigger_signal_by_IX(mvmnt_prob, ib, win_samples);
t = linspace(-win, win, size(trig_signal,2));
shadedErrorBar(t, trig_signal, {@nanmean,@nansem});
xlabel('Time from PE (s)')
ylabel('Movement prob.')
text(.95,.95,"n="+length(ib),'Units','normalized','HorizontalAlignment','right');

sgtitle(sprintf('%s_%s_opt_%d',exp_ID, epoch_type, params_opt),'interpreter','none');
filename = fullfile(dir_OUT, sprintf('%s_%s_opt_%d',exp_ID, epoch_type, params_opt));
saveas(fig, filename, 'jpg');
    
%%




%%


