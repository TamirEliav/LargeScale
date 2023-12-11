%%
clear
clc

%% graphic params
directions_clrs = {[0    0.3843    0.7451];[ 0.5216    0.2471         0]};
epoch_types_line_style = {'-','--',':','-.'};
epoch_type_clrs = {[.6 .1 .8],[.1 .8 .1]};

%% choose bats / sessions
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
clear exp_list
groupsummary(T,'bat_num')
if exist('bats_to_include','var')
    T = groupfilter(T,"bat_num",@(x)ismember(x,bats_to_include),'bat_num');
end
bats = unique(T.bat_num);

%% load and analyze data
nDirs = 2;
nbins = 50;
nSessions = height(T);
epoch_types = {'sleep','rest','sleep1','sleep2'};
nEpochTypes = length(epoch_types);
coverage_all = zeros(nSessions,nEpochTypes,nDirs,nbins); % [exp] X [nEpochTypes] X [2 directions] X [space]
ccc_all = zeros(nSessions,nDirs,nEpochTypes,nEpochTypes); % [exp] X [2 directions] X [epochType] X [nEpochTypes]
n_seqs_all = zeros(nSessions,nEpochTypes,nDirs); % [exp] X [nEpochTypes] X [2 directions]
sessions_duration = nan(nSessions,nEpochTypes);
for ii_exp = 1:nSessions
    % load exp data
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details','rest');
    params_opt = 11;
    [coverage_all(ii_exp,1,:,:), n_seqs_all(ii_exp,1,:)] = calc_coverage(exp_ID, 'sleep', params_opt,nbins);
    [coverage_all(ii_exp,2,:,:), n_seqs_all(ii_exp,2,:)] = calc_coverage(exp_ID, 'rest',  params_opt,nbins);
    [coverage_all(ii_exp,3,:,:), n_seqs_all(ii_exp,3,:)] = calc_coverage(exp_ID, 'sleep', params_opt,nbins,'epoch_num',1);
    [coverage_all(ii_exp,4,:,:), n_seqs_all(ii_exp,4,:)] = calc_coverage(exp_ID, 'sleep', params_opt,nbins,'epoch_num',2);
    ccc1 = corr(squeeze(coverage_all(ii_exp,:,1,:))');
    ccc2 = corr(squeeze(coverage_all(ii_exp,:,2,:))');
    ccc_all(ii_exp,1,:,:) = ccc1;
    ccc_all(ii_exp,2,:,:) = ccc2;

    sleep1_duration = diff(exp_get_sessions_ti(exp_ID,{'Sleep1'}))*1e-6;
    sleep2_duration = diff(exp_get_sessions_ti(exp_ID,{'Sleep2'}))*1e-6;
    sessions_duration(ii_exp,1) = sleep1_duration+sleep2_duration;
    sessions_duration(ii_exp,2) = sum(diff(exp.rest.ti,1,2).*1e-6);
    sessions_duration(ii_exp,3) = sleep1_duration;
    sessions_duration(ii_exp,4) = sleep2_duration;
end
replay_rate = sum(n_seqs_all,3) ./ sessions_duration;

%% calc shuffled corr
ccc_shuffled = {};
ccc_shuffled_sleep = {};
for ii_dir=1:2
    coverage_sleep = squeeze(coverage_all(:,1,ii_dir,:))';
    coverage_rest = squeeze(coverage_all(:,2,ii_dir,:))';
    coverage_sleep1 = squeeze(coverage_all(:,3,ii_dir,:))';
    coverage_sleep2 = squeeze(coverage_all(:,4,ii_dir,:))';
    
    ccc = corr(coverage_sleep,coverage_rest);
    mask = tril(true(size(ccc)),-1);
    sum(isnan(ccc(mask)),'all')
    ccc_shuffled{ii_dir} = ccc(mask);

    ccc = corr(coverage_sleep1,coverage_sleep2);
    mask = tril(true(size(ccc)),-1);
    sum(isnan(ccc(mask)),'all')
    ccc_shuffled_sleep{ii_dir} = ccc(mask);
end

%%
fig = figure;
fig.WindowState = 'maximized';
tiledlayout('flow','TileSpacing','compact')
nexttile
hold on
sdf = ccc_all(:,:,1,2);
histogram(sdf(:),'Normalization','pdf','FaceColor',[1 1 1].*0.5,'BinWidth',0.2)
histogram([ccc_shuffled{:}],'normalization','pdf','DisplayStyle','stairs','EdgeColor','k','lineWidth',2)
xlabel('corr')
ylabel('pdf')
legend("data","shuffle")
title('corr: Sleep vs. Rest')

nexttile
hold on
sdf = ccc_all(:,:,3,4);
histogram(sdf(:),'Normalization','pdf','FaceColor',[1 1 1].*0.5,'BinWidth',0.2)
histogram([ccc_shuffled_sleep{:}],'normalization','pdf','DisplayStyle','stairs','EdgeColor','k','lineWidth',2)
xlabel('corr')
ylabel('pdf')
legend("data","shuffle")
title('corr: Sleep1 vs. Sleep2')

nexttile
hold on
sdf1 = ccc_all(:,:,2,3);
sdf2 = ccc_all(:,:,2,4);
EDGES = linspace(-1,1,21);
histogram(sdf1(:),'Normalization','pdf','DisplayStyle','stairs','LineWidth', 2,'NumBins',21);
histogram(sdf2(:),'Normalization','pdf','DisplayStyle','stairs','LineWidth', 2,'NumBins',21);
[~,KS_pval] = kstest2(sdf1(:),sdf2(:));
RANKSUM_pval = ranksum(sdf1(:),sdf2(:));
SIGNRANK_pval = signrank(sdf1(:),sdf2(:));
if KS_pval>0.05
    KS_str = 'KS-test: n.s.';
end
if RANKSUM_pval>0.05
    RANKSUM_str = 'ranksum-test: n.s.';
end
if SIGNRANK_pval>0.05
    SIGNRANK_str = 'signrank-test: n.s.';
end
text(.4,.7,KS_str,'Units','normalized','FontSize',16)
text(.4,.8,RANKSUM_str,'Units','normalized','FontSize',16)
text(.4,.9,SIGNRANK_str,'Units','normalized','FontSize',16)
xlabel('corr')
ylabel('pdf')
legend("PRE-sleep vs rest","POST-sleep vs rest",'Location','best')

nexttile
plot(sessions_duration,'-o')
legend(epoch_types,'location','best')
xlabel('session no.')
ylabel('Session duration (s)')

nexttile
plot(replay_rate,'-o')
legend(epoch_types,'location','best')
xlabel('session no.')
ylabel('Replay rate (events/s)')
set(gca,'YScale','log')

nexttile
hold on
plot(replay_rate(:,3),replay_rate(:,4),'ok')
plot([1e-3 1e0],[1e-3 1e0],'Color',[1 1 1].*0.5)
pval = signrank(replay_rate(:,3),replay_rate(:,4));
text(0.1,0.9,sprintf('P signrank = %.g',pval),'units','normalized','FontSize',16)
axis square
axis equal
hax=gca;
hax.XLim(1) = 0;
hax.YLim(1) = 0;
hax.XScale = 'log';
hax.YScale = 'log';
refline(1,0)
xlabel('PRE-Sleep replay rate (events/s)')
ylabel('POST-Sleep replay rate (events/s)')

nexttile
hold on
plot(replay_rate(:,1),replay_rate(:,2),'ok')
plot([1e-3 1e0],[1e-3 1e0],'Color',[1 1 1].*0.5)
pval = signrank(replay_rate(:,1),replay_rate(:,2));
text(0.1,0.9,sprintf('P signrank = %.g',pval),'units','normalized','FontSize',16)
axis square
axis equal
hax=gca;
hax.XLim(1) = 0;
hax.YLim(1) = 0;
hax.XScale = 'log';
hax.YScale = 'log';
refline(1,0)
xlabel('Sleep replay rate (events/s)')
ylabel('Rest replay rate (events/s)')

nexttile
hold on
plot(replay_rate(:,3),replay_rate(:,2),'ok')
plot([1e-3 1e0],[1e-3 1e0],'Color',[1 1 1].*0.5)
pval = signrank(replay_rate(:,3),replay_rate(:,2));
text(0.1,0.9,sprintf('P signrank = %.g',pval),'units','normalized','FontSize',16)
axis square
axis equal
hax=gca;
hax.XLim(1) = 0;
hax.YLim(1) = 0;
hax.XScale = 'log';
hax.YScale = 'log';
refline(1,0)
xlabel('PRE-Sleep replay rate (events/s)')
ylabel('Rest replay rate (events/s)')

nexttile
hold on
plot(replay_rate(:,4),replay_rate(:,2),'ok')
plot([1e-3 1e0],[1e-3 1e0],'Color',[1 1 1].*0.5)
pval = signrank(replay_rate(:,4),replay_rate(:,2));
text(0.1,0.9,sprintf('P signrank = %.g',pval),'units','normalized','FontSize',16)
axis square
axis equal
hax=gca;
hax.XLim(1) = 0;
hax.YLim(1) = 0;
hax.XScale = 'log';
hax.YScale = 'log';
refline(1,0)
xlabel('POST-Sleep replay rate (events/s)')
ylabel('Rest replay rate (events/s)')

dir_OUT = 'F:\sequences\decoded_figs\replay_coverage';
filename = 'comparing_PRE_POST_sleep_replay_coverage_corr_and_rate';
file_OUT = fullfile(dir_OUT,filename);
exportgraphics(fig,[file_OUT '.pdf'])
% saveas(fig,file_OUT,'pdf')

%% save processed data
dir_OUT = 'E:\Tamir\work\PROJECTS\LargeScale\paper_replay\data_prepared_for_figures';
filename = 'replay_coverage';
file_OUT = fullfile(dir_OUT,filename);
save(file_OUT,'T','coverage_all','ccc_all','n_seqs_all','ccc_shuffled','epoch_types','sessions_duration','replay_rate');

%% sort exp_list by #replays
n_seq_pooled = sum(n_seqs_all(:,1:2,:),[2 3]);
[~,sorted_IX] = sort(n_seq_pooled,'descend');

%% plot coverage (fig per sessions)
xbins = linspace(0,1,nbins);
fig=figure;
fig.Units = 'centimeters';
fig.Position([3 4]) = [15 6];
axes('Units','normalized','Position',[0.1 0.2 0.85 0.65])
[~,sorted_IX] = sort(n_seq_pooled,'descend');
for ii = 1:length(sorted_IX)
    ii_exp = sorted_IX(ii);
    exp_ID = T.exp_ID{ii_exp};
    fprintf('%d\t%d\t%d\t%s\n',ii,n_seq_pooled(ii_exp),ii_exp,exp_ID);
    cla
    hold on
    lw = 2;
    for ii_epoch_type = 1:2
        for ii_dir = 1:2
            c = squeeze(coverage_all(ii_exp,ii_epoch_type,ii_dir,:));
            plot(xbins, c,'Color',directions_clrs{ii_dir},'LineStyle',epoch_types_line_style{ii_epoch_type},'LineWidth',lw);
        end
    end
%     plot(xbins, squeeze(c(1,1,:)),'-','Color',directions_clrs{1},'LineWidth',lw)
%     plot(xbins, squeeze(c(2,1,:)),'--','Color',directions_clrs{1},'LineWidth',lw)
%     plot(xbins, squeeze(c(1,2,:)),'-','Color',directions_clrs{2},'LineWidth',lw)
%     plot(xbins, squeeze(c(2,2,:)),'--','Color',directions_clrs{2},'LineWidth',lw)
    text(0.05,1.05,sprintf('#%d %s',ii,exp_ID),'Interpreter','none','HorizontalAlignment','left','Units','normalized');
    xlabel('Position (norm.)')
    ylabel('Replay coverage (counts)')
    hl=legend({'Sleep dir 1','Sleep dir 2','Awake dir 1','Awake dir 2'},'NumColumns',2);
    hl.Units = 'normalized';
    hl.Position([1 2]) = [0.6 0.85];
    dir_OUT = 'F:\sequences\decoded_figs\replay_coverage';
    filename = sprintf('ii_%.3d_%d_replay_coverage_%s',ii,ii_exp,exp_ID);
    file_OUT = fullfile(dir_OUT,filename);
    saveas(fig,file_OUT,'pdf')
end

%% plot full-session (fig per sessions)
fig=figure;
fig.Units = 'centimeters';
fig.Position = [5 2 15 15];
clear panels
panels(1) = axes('Units','normalized','Position',[0.2 0.2 0.25 0.7143]);
panels(2) = axes('Units','normalized','Position',[0.55 0.2 0.25 0.7143]);
% sorted_IX = 1:82;
[~,sorted_IX] = sort(n_seq_pooled,'descend');
replay_directionality = nan(height(T),2);
% seqs_all
for ii = 1:length(sorted_IX)
    ii_exp = sorted_IX(ii);
    exp_ID = T.exp_ID{ii_exp};
    fprintf('%d\t%d\t%d\t%s\n',ii,n_seq_pooled(ii_exp),ii_exp,exp_ID);
    cla
    hold on
    lw = 2;
    for ii_epoch_type = 1:2
        axes(panels(ii_epoch_type));
        hax=gca;
        cla 
        hold on
        epoch_type = epoch_types{ii_epoch_type};
        event_type = 'posterior';
        [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
        [~, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
        events(~TF) = [];
        if isempty(events)
            continue;
        end
        seqs= [events.seq_model];
        ndir1 = sum([seqs.direction]==1);
        ndir2 = sum([seqs.direction]==-1);
        replay_directionality(ii_exp,ii_epoch_type) = (ndir1-ndir2)/(ndir1+ndir2);
        
        seqs_edges = [seqs.start_pos_norm; seqs.end_pos_norm];
        seqs_IX = 1:length(seqs);
        h=plot(seqs_edges,[seqs_IX;seqs_IX],'-');
        dir_1_IX = [seqs.state_direction]==1;
        dir_2_IX = [seqs.state_direction]==-1;
        [h(dir_1_IX).Color] = disperse(repelem(directions_clrs(1),length(dir_1_IX)));
        [h(dir_2_IX).Color] = disperse(repelem(directions_clrs(2),length(dir_2_IX)));
        scatter(seqs_edges(1,:),seqs_IX, 3, [seqs.state_direction]==-1, "filled");
        hax.Colormap = cell2mat(directions_clrs);
        hax.XTick = [0 1];
        hax.YTick = [1 10*ceil(length(seqs)/10)];
        hax.XLim = [0 1];
        hax.XRuler.TickLabelGapOffset = -1;
        hax.YRuler.TickLabelGapOffset = 1;
        xlabel('Position (norm.)', 'Units','normalized', 'Position',[0.5 -0.017]);
        ylabel('Replay event no.', 'Units','normalized', 'Position',[-0.09 .5]);
        title(sprintf('Example %s session',epoch_type), 'Units','normalized', 'Position',[0.47 1],'FontWeight','normal');

        text(0.05,1.08,sprintf('#%d %s',ii,exp_ID),'Interpreter','none','HorizontalAlignment','left','Units','normalized');
    end
    dir_OUT = 'F:\sequences\decoded_figs\replay_coverage\full_sessions';
    mkdir(dir_OUT);
    filename = sprintf('ii_%.3d_%d_replay_coverage_%s',ii,ii_exp,exp_ID);
    file_OUT = fullfile(dir_OUT,filename);
%     saveas(fig,file_OUT,'pdf')
end
%% seems like in bat 184 the first days (this is a novelty bat) are VERY directional towards ball 1 ("home base")
figure
tiledlayout('flow')
nexttile
hold on
for ii_epoch_type = 1:2
    h=plot(sum(n_seqs_all(:,ii_epoch_type,:),3), replay_directionality(:,ii_epoch_type),'o','Color',epoch_type_clrs{ii_epoch_type});
    row = dataTipTextRow('Sexp_ID',T.exp_ID);
    h.DataTipTemplate.DataTipRows(end+1) = row;
    h.DataTipTemplate.Interpreter = 'none';
end
hax=gca;
hax.XScale = 'log'
nexttile
plot(replay_directionality(:,1),replay_directionality(:,2),'.')
axis equal
axis square
refline(1,0)

for ii_bat = 1:length(bats)
    nexttile
    hold on
    bat_num = bats(ii_bat);
    bat_str = sprintf('%.4d',bat_num)
    TF = contains(T.exp_ID,bat_str);
    for ii_epoch_type = 1:2
        bubblechart(1:sum(TF),replay_directionality(TF,1),sum(n_seqs_all(TF,ii_epoch_type,:),3),epoch_type_clrs{ii_epoch_type})
    end
    bubblesize([3 10])
    bubblelim([ min(sum(n_seqs_all(:,ii_epoch_type,:),3)) 
                max(sum(n_seqs_all(:,ii_epoch_type,:),3)) ]);
    ylim([-1 1]*1.1)
    title(bat_str)
end

%%
function [coverage, n_seqs] = calc_coverage(exp_ID, epoch_type, params_opt,nbins,opts)
arguments 
    exp_ID
    epoch_type
    params_opt
    nbins
    opts.epoch_num = []
end

coverage = zeros(2,nbins);
n_seqs = zeros(2,1);

[events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
if isempty(events)
    return;
end
seqs = [events.seq_model];
[seqs, TF] = decoding_apply_seq_inclusion_criteria(seqs);
events(~TF)=[];
if isempty(seqs)
    return;
end
if ~isempty(opts.epoch_num)
    TF = ismember([events.epoch_num], opts.epoch_num);
    events(~TF)=[];
    seqs(~TF)=[];
end
if isempty(seqs)
    return;
end

seqs_edges = [seqs.start_pos_norm; seqs.end_pos_norm];
seqs_edges = [min(seqs_edges ); max(seqs_edges )]';

xbins = linspace(0,1,nbins);
directions = [1 -1];
for ii_dir = 1:2
    direction = directions(ii_dir);
    coverage(ii_dir,:) = sum( xbins>seqs_edges([seqs.state_direction]==direction,1) & ...
                              xbins<seqs_edges([seqs.state_direction]==direction,2),1);
    n_seqs(ii_dir) = sum([seqs.state_direction]==direction);
end



end