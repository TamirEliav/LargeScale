%%
clear 
clc

%% folders
dir_OUT = 'F:\sequences\figures\novelty_exposure';

%% params
epoch_types = {'sleep','rest'};
days_types = {'Early (days 1-5)','Late (days>5)'};
epoch_type_clrs = {[.6 .1 .8],[.1 .8 .1]};
line_styles = {'-','--'};

%% properties to plot
features_names = {'duration';'compression';'distance';'distance_norm';};
xlable_strs = {
    'Replay duration (s)';
    {'Compression ratio';'(replay speed / flight speed)'};
    'Replay distance (m)';
    {'Replay distance','(norm. to environment size)'};
    };
ft_xlimits = [
    -0.0950    1.9950;
    -1.9500   40.9500;
    0.4500   56.5500;
    -0.0210    0.4410
    ];

%% arrange sessions to load
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
novelty_exposure_bats = [184 2382];
T(~ismember(T.bat_num, novelty_exposure_bats),:)=[];
exp_list = T.exp_ID;
early_days_exp_list = {
    'b0184_d191129',
    'b0184_d191130',
    'b0184_d191201',
    'b2382_d190623',
    'b2382_d190624',
    'b2382_d190627',
    };
early_days_IX = find(ismember(T.exp_ID, early_days_exp_list));
late_days_IX = find(~ismember(T.exp_ID, early_days_exp_list));
early_late_IX = {early_days_IX,late_days_IX};

%% load data
events_all_per_session = {};
flight_speed_all = zeros(size(exp_list));
for ii_epoch_type = 1:length(epoch_types)
    for ii_exp = 1:length(exp_list)
        exp_ID = exp_list{ii_exp};
        exp = exp_load_data(exp_ID,'details','flight');
        flight_speed_all(ii_exp) = mean(abs([exp.flight.speed_traj.vel_median_median]));
        epoch_type = epoch_types{ii_epoch_type};
        params_opt = 11;
        event_type = 'posterior';
        [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
        [~, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
        events(~TF) = [];
        if isempty(events)
            continue;
        end
        [events.recordingArena] = deal(exp.details.recordingArena);
        event_struct = events(1,[]);
        [events(2:end).prev_event] = disperse(events(1:end-1));
        events(1).prev_event = event_struct;
        events_all_per_session{ii_epoch_type}{ii_exp} = events;
    end
end

%% histograms
figure(Units="centimeters",Position=[5 5 25 13]);
tiledlayout(2,4,'TileSpacing','compact');
for ii_EL = 1:length(early_late_IX)
    ii_EL
    for ii_fn = 1:length(features_names)
        fn = features_names{ii_fn};
        nexttile
        cla
        hold on
        lw = 1.1;
        for ii_epoch_type = 1:length(epoch_types)
            events = [events_all_per_session{ii_epoch_type}{early_late_IX{ii_EL}}];
            seqs = [events.seq_model];
            X = [seqs.(fn)];
            histogram(X,'DisplayStyle','stairs','Normalization','pdf','EdgeColor',epoch_type_clrs{ii_epoch_type},'LineWidth',lw);
            text(0.5,0.9-ii_epoch_type*0.1, "\itn_{" + epoch_types{ii_epoch_type} + "}="+length(X),'units','normalized')
        end
        xlim(ft_xlimits(ii_fn,:))
        title(days_types{ii_EL})
        xlabel(xlable_strs{ii_fn});
        if ii_fn == 1
            ylabel('Probability density');
        end
        hax=gca;
    %     hax.TickLength(1) = [0.025];
        hax.XRuler.TickLength(1) = 0.03;
        hax.XRuler.TickLabelGapOffset = -1.2;
        hax.YRuler.TickLabelGapOffset = 1;
        fprintf('%s: median %.2g, mean = %.2g\n',fn,median(X),mean(X))
    end
end
figname = 'early_vs_late_exposure_days_repaly_properties_hist';
filename = fullfile(dir_OUT,figname);
saveas(gcf,filename,'tif')

%% boxplots
figure(Units="centimeters",Position=[5 5 30 20]);
tiledlayout(2,4,'TileSpacing','loose');
for ii_epoch_type = 1:length(epoch_types)
    for ii_fn = 1:length(features_names)
        fn = features_names{ii_fn};
        nexttile
        cla
        hold on
        lw = 1.1;
        X={};
        for ii_EL = 1:length(early_late_IX)
            events = [events_all_per_session{ii_epoch_type}{early_late_IX{ii_EL}}];
            seqs = [events.seq_model];
            X{ii_EL} = [seqs.(fn)];
            G{ii_EL} = ones(size(seqs)).*ii_EL;
        end
        ranksum_pval = ranksum(X{1},X{2});
        [~,KS_pval] = kstest2(X{1},X{2});
        boxplot([X{:}],[G{:}],'labels',days_types)
        ylabel(xlable_strs{ii_fn});
        title({epoch_types{ii_epoch_type}, ...
            sprintf('ranksum p=%.2g',ranksum_pval), ...
            sprintf('KS p=%.2g',KS_pval)})
    end
end
figname = 'early_vs_late_exposure_days_repaly_properties_boxplot';
filename = fullfile(dir_OUT,figname);
saveas(gcf,filename,'tif')

%% boxplots
for ii_bat = 1:length(novelty_exposure_bats)
    figure(Units="centimeters",Position=[5 5 30 20]);
    tiledlayout(2,4,'TileSpacing','loose');
    bat_num = novelty_exposure_bats(ii_bat);
    for ii_epoch_type = 1:length(epoch_types)
        for ii_fn = 1:length(features_names)
            fn = features_names{ii_fn};
            nexttile
            cla
            hold on
            lw = 1.1;
            X={};
            for ii_EL = 1:length(early_late_IX)
                IX1 = early_late_IX{ii_EL};
                IX2 = find(T.bat_num == bat_num);
                IX = intersect(IX1,IX2);
                events = [events_all_per_session{ii_epoch_type}{IX}];
                seqs = [events.seq_model];
                X{ii_EL} = [seqs.(fn)];
                G{ii_EL} = ones(size(seqs)).*ii_EL;
            end
            ranksum_pval = ranksum(X{1},X{2});
            [~,KS_pval] = kstest2(X{1},X{2});
            boxplot([X{:}],[G{:}],'labels',days_types)
            ylabel(xlable_strs{ii_fn});
            title({epoch_types{ii_epoch_type}, ...
                sprintf('ranksum p=%.2g',ranksum_pval), ...
                sprintf('KS p=%.2g',KS_pval)})
        end
    end
    sgtitle("bat"+bat_num)
    figname = sprintf('early_vs_late_exposure_days_repaly_properties_boxplot_bat_%d',bat_num);
    filename = fullfile(dir_OUT,figname);
    saveas(gcf,filename,'tif')
end

%% flight speed early vs late days
rng(0)
figure
tiledlayout("flow",'TileSpacing','compact')
for ii_bat = 1:length(novelty_exposure_bats)
    bat_num = novelty_exposure_bats(ii_bat);
    early_IX = find(T.bat_num == bat_num & ismember(T.exp_ID, early_days_exp_list));
    late_IX = find(T.bat_num == bat_num & ~ismember(T.exp_ID, early_days_exp_list));
    nexttile
    hold on
    y = [flight_speed_all(early_IX); flight_speed_all(late_IX)];
    x = [1.*ones(length(early_IX),1);
         2.*ones(length(late_IX),1) ];
    swarmchart(x,y)
    xticks([1 2])
    xticklabels(days_types)
    title("bat"+bat_num)
    ylabel('Flight speed (m/s')
%     histogram(flight_speed_all(early_IX),'DisplayName',days_types{1})
%     histogram(flight_speed_all(late_IX),'DisplayName',days_types{2})
end
figname = 'early_vs_late_exposure_days_flight_speeds';
filename = fullfile(dir_OUT,figname);
saveas(gcf,filename,'tif')


%%
close all
directionality_contrast_index = [];
directionality_fraction = [];
directionality_binom_pval = [];
nSeqs = [];
for ii_epoch_type = 1:length(epoch_types)
    for ii_exp = 1:length(exp_list)
        events = events_all_per_session{ii_epoch_type}{ii_exp};
        seqs = [events.seq_model];
        n1 = sum([seqs.direction]==1);
        n2 = sum([seqs.direction]==-1);
        directionality_contrast_index(ii_epoch_type,ii_exp) = (n1-n2)/(n1+n2);
        directionality_fraction(ii_epoch_type,ii_exp) = (n1)/(n1+n2);
        directionality_binom_pval(ii_epoch_type,ii_exp) = myBinomTest(n1,n1+n2,0.5);
        nSeqs(ii_epoch_type,ii_exp) = length(seqs);
    end
end
vals = directionality_contrast_index;
% vals = directionality_fraction;
% vals = directionality_binom_pval;
vals(nSeqs<20)=nan;
[G,GID] = findgroups(T.bat_num);
figure
hold on
splitapply(@(x)plot(x,'-o'),vals',G);
% splitapply(@(x,n)bubblechart(1:length(x),x,n),vals',nSeqs',G);
ylim([-1 1])
% ylim([0 1])
xlabel('#session')
ylabel('direcitonality')
legend([string(epoch_types)+GID]');
figname = 'novelty_directionality';
filename = fullfile(dir_OUT,figname);
saveas(gcf,filename,'tif')







%%

%%

