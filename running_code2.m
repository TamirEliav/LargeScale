%%
exps_IDs = exp_t.exp_ID;

%%
exps = cellfun(@(exp_ID)( exp_load_data(exp_ID,'details','PE') ), exps_IDs);
PEs = {exps.PE};
PE_all = [exps.PE];
nPE = cellfun(@length, PEs);
mean_FR  = cellfun(@(PE)(mean([PE.peak_FR])), PEs);
mean_zFR  = cellfun(@(PE)(mean([PE.peak_zFR])), PEs);
order_feature = nPE .* mean_zFR;
[~,ranks]  = ismember(order_feature,flip(unique(order_feature)));
% [~,sort_IX] = sort(order_feature,'descend');
[exps.nPE] = disperse(nPE);
[exps.mean_zFR] = disperse(mean_zFR);
[exps.order_feature] = disperse(order_feature);
[exps.rank] = disperse(ranks);

%%
figure
tiledlayout('flow')
nexttile
h=scatter(nPE,mean_FR,20,order_feature,'filled');
h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('ranks',ranks);
h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('feature value',order_feature);
% lsline
xlabel('No. of events')
ylabel('Mean firing rate (Hz)')
colorbar

nexttile
h=scatter(nPE,mean_zFR,20,order_feature,'filled');
h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('ranks',ranks);
h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('feature value',order_feature);
% lsline
xlabel('No. of events')
ylabel('Mean firing rate (z)')
colorbar

%% sort exps
% exps=exps(sort_IX);

T1 = table( arrayfun(@(exp)(exp.details.exp_ID), exps,'UniformOutput',0),...
            arrayfun(@(exp)(exp.details.batNum), exps,'UniformOutput',1),...
            datetime(arrayfun(@(exp)(exp.details.date), exps,'UniformOutput',1),'Format','yyMMdd'),...
            arrayfun(@(exp)(exp.nPE), exps,'UniformOutput',1),...
            arrayfun(@(exp)(exp.mean_zFR), exps,'UniformOutput',1),...
            arrayfun(@(exp)(exp.order_feature), exps,'UniformOutput',1),...
            arrayfun(@(exp)(exp.rank), exps,'UniformOutput',1),...
            'VariableNames',{'exp_ID','batNum','date','nPE','mean_zFR','order_feature','rank'});
T1.Properties.RowNames = string(datestr(T1.date,'yymmdd'));

%%
cells_table_filename = "L:\Analysis\Code\inclusion_lists\cells_summary.xlsx";
cells_t = readtable(cells_table_filename);
T2 = groupcounts(cells_t,'date');
T2.Properties.VariableNames(2) = {'nCells'};
T2.date = datetime(T2.date,'Format','yyMMdd');
T2.Properties.RowNames = string(datestr(T2.date,'yymmdd'));

%%
T = outerjoin(T2,T1,'Keys','Row');
T = sortrows(T,'nCells','descend');
clc
T

%% creating replay examples
decoding_plot_event('b9861_d180527', 'sleep', 11, 824, 0.5)
decoding_plot_event('b9861_d180527', 'sleep', 11, 681, 0.5)
decoding_plot_event('b9861_d180527', 'sleep', 11, 503, 0.5)
decoding_plot_event('b9861_d180526', 'sleep', 11, 27,  0.5)
decoding_plot_event('b9861_d180526', 'sleep', 11, 42,  0.5)
decoding_plot_event('b9861_d180526', 'sleep', 11, 113, 0.5)
decoding_plot_event('b9861_d180526', 'sleep', 11, 113, 1)
decoding_plot_event('b9861_d180526', 'sleep', 11, 113, 1.5)
decoding_plot_event('b9861_d180526', 'sleep', 11, 113, 2)
decoding_plot_event('b9861_d180526', 'sleep', 13, 113, 1)
decoding_plot_event('b9861_d180526', 'sleep', 11, 132, 0.5)
decoding_plot_event('b9861_d180526', 'sleep', 11, 132, 1)
decoding_plot_event('b0184_d191201', 'sleep', 11, 136, 0.5)
decoding_plot_event('b0184_d191130', 'sleep', 11, 53,  0.5)
decoding_plot_event('b0184_d191130', 'sleep', 11, 333,  0.5)
decoding_plot_event('b0184_d191130', 'sleep', 11, 350,  0.5)

%%
decoding_plot_event('b9861_d180526', 'sleep', 14, 113, 1);
decoding_plot_event('b9861_d180527', 'sleep', 11, 824, 0.42);
decoding_plot_event('b0184_d191130', 'sleep', 11, 53,  [-.32 .11]);

%%
% decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',49);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',98);

%%
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',15);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',22);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',23);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',25);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',44);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',49);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',50);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',51);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',52);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',55);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',57);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',60);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',70);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',75);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',84);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',90);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',98);
decoding_seq_quantify_plot_event('b9861_d180526', 'sleep', 11, 'posterior',105);

%%
clc
X1 = [3 1 10 8 16];
u = [1 1 1 1 1];
% X2 = conv(X1,u,'full');
X2 = imfilter(X1,u,'symmetric','conv','same');
% X2 = imfilter(X1,u,'replicate','conv','same');
disp(X1)
disp(X2)


%%
clear; clc
exp_ID = 'b0184_d191130';
exp = exp_load_data(exp_ID, 'details', 'path');
dir_IN = exp.path.spikes_raw;
limits_ts = [];
TT = exp.details.refCh(1);
ch = exp.details.refCh(2);
csc_file = sprintf('spikes_%s_TT%d_ch%d.ncs',exp_ID,TT,ch);
csc_file = fullfile(dir_IN,csc_file);
csc_ref_ch = Nlx_csc_read(csc_file, limits_ts);
data = zeros([size(exp.details.activeChannels),length(csc_ref_ch)]);
%%
tic
for TT = 1:exp.details.numTT
    parfor ch = 1:4
        fprintf('TT%d_ch%d\n',TT,ch);
        % load ch data
        csc_file = sprintf('spikes_%s_TT%d_ch%d.ncs',exp_ID,TT,ch);
        csc_file = fullfile(dir_IN,csc_file);
        if ~exist(csc_file,'file')
            continue;
        end
        data(TT,ch,:) = Nlx_csc_read(csc_file,limits_ts);
    end
end
toc

%%
sdf=median(data,[1 2]);


%% select best TT for dislaying ripple traces
events_all = ripples.all;
IX = [events_all.peak_IX];
IX = IX+[-10:10]';
pripple_TT_at_peaks = pripple_TT(IX,:);
mean(pripple_TT_at_peaks)
n_contrib_events = zeros(1,nTT);
contrib_events = false(length(events_all),nTT);
for TT=1:nTT
    events_TT = ripples.by_TT{TT};
    if isempty(events_TT)
        continue;
    end
    t1 = [events_all.peak_ts];
    t2 = [events_TT.peak_ts];
    tdiff = abs(t1-t2');
    thr = 100; % usec
    TF = tdiff < thr;
    fprintf('TT %d: %d\n',TT,sum(any(TF)))
    contrib_events(:,TT) = any(TF);
    n_contrib_events(TT) = sum(any(TF));
end
figure
subplot(131);bar(1:nTT,n_contrib_events); ylabel('num events')
subplot(132);bar(1:nTT,mean(pripple_TT_at_peaks)); ylabel('pripple')
subplot(133);bar(1:nTT,n_contrib_events.*mean(pripple_TT_at_peaks)); ylabel('product')


%%
IX = [events_all.peak_IX];
IX = IX+[-200:200]';
% signal = LFP2;
% signal = sum(ripple,3);
signal = pripple_TT;
sdf = signal(IX,:);
sdf = reshape(sdf,size(IX,1),size(IX,2),nTT);

%%
close all
for TT=TT_to_use
    figure
    sdf2=squeeze(sdf(:,:,TT));
    sdf3=sdf2;
%     sdf3 = sdf3+[1:nEvents].*25;
%     plot(sdf3)
    imagesc(sdf3')
    title("TT"+TT)
end

%%
% sdf4= squeeze(sum(sdf,1));
sdf4 = squeeze(max(sdf,[],1));
% sdf4(isnan(sdf4))=0;
sdf4 = sdf4(:,TT_to_use);
[IDX C] = kmeans(sdf4,2);

[coeff,score,latent,tsquared,explained,mu] = pca(sdf4);
[~,sorted_IX] = sort(score(:,2));

% figure
% scatter(score(:,1),score(:,2),5,IDX,'filled')

figure
imagesc(coeff(:,[1 2]))
yticklabels(TT_to_use)
colormap bone
colorbar

figure
imagesc(contrib_events(sorted_IX,TT_to_use))

figure
subplot(121)
imagesc(squeeze(sdf(:,sorted_IX,1))')
% imagesc(squeeze(sdf(:,:,1))')
subplot(122)
imagesc(squeeze(sdf(:,sorted_IX,12))')
% imagesc(squeeze(sdf(:,:,12))')
















%%
