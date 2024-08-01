%%
clear
clc

%%
dir_out = 'E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\internal_review_figures';

%%
epoch_type_clrs = [.6 .1 .8;.1 .8 .1];
IRI_bin_size = 0.1;
IRI_bin_limits = [0 20];
AC_bin_limits = [0 8];
AC_bin_width = 0.050;
AC_bin_edges = AC_bin_width/2*[-1 1]+AC_bin_limits;
AC_bin_edges = AC_bin_edges(1):AC_bin_width:AC_bin_edges(end);
AC_bin_centers = edges2centers(AC_bin_edges);

%% load all replays
epoch_types = {'sleep','rest'};
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
bats = unique(T.bat_num);
nExp = length(exp_list);
nEpochTypes = length(epoch_types);
events_all_per_session = cell(nEpochTypes, nExp);
IRI_all = cell(nEpochTypes, nExp);
replay_AC_all = zeros(nEpochTypes, nExp,length(AC_bin_centers));
for ii_exp = 1:nExp
    for ii_epoch_type = 1:nEpochTypes
        
        %% load events
        exp_ID = exp_list{ii_exp};
        exp = exp_load_data(exp_ID,'details');
        epoch_type = epoch_types{ii_epoch_type};
        fprintf('%d: %s (%s)\n',ii_exp,exp_ID, epoch_type);
        params_opt = 11;
        event_type = 'posterior';
        [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
        [~, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
        events(~TF) = [];
        if isempty(events)
            continue;
        end

        %% calc inter-replay-intervals
        replay_ts = mean([[events.start_ts];[events.end_ts]]);
        replay_ts = replay_ts .* 1e-6;
        g = findgroups([events.epoch_num]);
        IRI = splitapply(@(x){diff(x)},replay_ts,g);
        IRI = [IRI{:}];
%         IRI = IRI .* 1e-6; % convert to seconds
        IRI_all{ii_epoch_type,ii_exp} = IRI;
        replay_ts_diff = replay_ts-replay_ts';
        replay_AC = histcounts(replay_ts_diff(:),AC_bin_edges);
        replay_AC_all(ii_epoch_type,ii_exp,:) = replay_AC;

        %%
        events_all_per_session{ii_epoch_type,ii_exp} = events;
        
    end
end

%% save results to a mat file
filename = 'inter_replay_interval_all';
file_OUT = fullfile(dir_out,filename);
save(file_OUT,'T','exp_list','epoch_types','IRI_all');

%% plot ISI hist per session/bat/grand avg
close all
fig=figure
fig.WindowState = 'maximized';
tiledlayout("flow")
% grand avg
nexttile
hold on
hax=gca;
hax.ColorOrder = epoch_type_clrs;
for ii_epoch_type = 1:nEpochTypes
    x = [IRI_all{ii_epoch_type,:}];
    hh=histogram(x,'Normalization','pdf','DisplayStyle','Stairs','BinWidth',IRI_bin_size);
    hh.LineWidth = 2;
end
x = [IRI_all{:,:}];
hh=histogram(x,'Normalization','pdf','DisplayStyle','Stairs','BinWidth',IRI_bin_size,'EdgeColor','k');
hh.LineWidth = 2;
xlim(IRI_bin_limits)
legend({epoch_types{:},'pooled'} )
xlabel('Inter-replay-interval (s)')
ylabel('PDF')
title('Grand average')

% per session (separated sleep/rest)
nexttile
hold on
for ii_epoch_type = 1:nEpochTypes
    for ii_exp = 1:nExp
        x = [IRI_all{ii_epoch_type,ii_exp}];
        hh=histogram(x,'Normalization','count','DisplayStyle','Stairs','BinWidth',IRI_bin_size);
        hh.LineWidth = 1;
    end
end
xlim(IRI_bin_limits)
xlabel('Inter-replay-interval (s)')
ylabel('counts')
title('Per sessions (separated to sleep/rest)')

% per bat
for ii_bat = 1:length(bats)
    nexttile
    hold on
    hax=gca;
    hax.ColorOrder = epoch_type_clrs;
    bat_num = bats(ii_bat);
    bat_sessions_IX = T.bat_num==bat_num;
    for ii_epoch_type = 1:nEpochTypes
        x = [IRI_all{ii_epoch_type,bat_sessions_IX}];
        hh=histogram(x,'Normalization','pdf','DisplayStyle','Stairs','BinWidth',IRI_bin_size);
        hh.LineWidth = 2;
    end
    xlim(IRI_bin_limits)
    title("bat "+bat_num);
    xlabel('Inter-replay-interval (s)')
    ylabel('PDF')
end
sgtitle('Inter-replay-interval histograms','FontSize',20,'FontWeight','Bold')
filename = 'inter_replay_interval_hist';
file_OUT = fullfile(dir_out,filename);
saveas(fig,file_OUT,'jpg')

%% plot Autocorrelation
% nexttile
% hold on
% hax=gca;
% hax.ColorOrder = epoch_type_clrs;
% for ii_epoch_type = 1:nEpochTypes
%     x = squeeze(replay_AC_all(ii_epoch_type,:,:));
%     x = sum(x);
%     x(1)=nan;
%     plot(AC_bin_centers, x,'LineWidth',2);
% end
% legend(epoch_types )




%%