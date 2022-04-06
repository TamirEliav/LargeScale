%% PhD thesis figure 4.4 - replay - forward vs reverse

%%
close all
clear 
clc

%% options
% bats_to_include = [34 148 9861 2289];
% bats_to_include = [194 184 2382];
% bats_to_include = 34;
% bats_to_include = 148;
% bats_to_include = 9861;
% bats_to_include = 2289;
% bats_to_include = 194;
% bats_to_include = 184;
% bats_to_include = 2382;
%     bat_num    GroupCount
%     _______    __________
% 
%      2289           1    
%        34           2    
%      9861           4    
%       148          15    
%       184          17    
%       194          17    
%      2382          26    

%% define output files
res_dir = 'E:\Tamir\PhD\Thesis\resources\ch_4_seq';
mkdir(res_dir)
fig_name_str = 'Fig_4_5_replay_forward_vs_reverse';
fig_caption_str = ' ';
log_name_str = [fig_name_str '_log_file' '.txt'];
log_name_str = strrep(log_name_str , ':', '-');
log_name_str = strrep(log_name_str , ' ', '_');
log_name_out = fullfile(res_dir, log_name_str);

%% open log file
diary off
diary(log_name_out)
diary on
disp('Log file');
disp(['created: ', datestr(clock)]);
disp('======================================================');
disp([fig_name_str ':' fig_caption_str]);   
disp('======================================================');
disp('');

%% create figure
% figure_size_cm = [21.0 29.7]; % ~A4
figure_size_cm = [21.6 27.9]; % ~US letter
fig = figure;
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',7);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'DefaultAxesUnits','centimeters');
set(gcf,'PaperType','usletter')
% set(gcf,'PaperType','<custom>');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 figure_size_cm]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]); % position on screen...
set(gcf, 'Renderer', 'painters');
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

% create panels
panels_size = [4 4];
panels_A(1) = axes('position', [1.5 19 2 4]);
panels_B(1) = axes('position', [5 19 panels_size]);
panels_B(2) = axes('position', [10 19 panels_size]);
panels_B(3) = axes('position', [15 19 panels_size]);
panels_B(4) = axes('position', [5 13.5 panels_size]);
panels_B(5) = axes('position', [10 13.5 panels_size]);

%% choose bats / sessions
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
clear exp_list
groupsummary(T,'bat_num')
if exist('bats_to_include','var')
    T = groupfilter(T,"bat_num",@(x)ismember(x,bats_to_include),'bat_num');
end
bats = unique(T.bat_num)

%% load data
events = {};
for ii_exp = 1:height(T)
    % load exp data
    exp_ID = T.exp_ID{ii_exp};
%     exp = exp_load_data(exp_ID,'details','path','ripples','MUA','PE');
    epoch_type = 'sleep';
%     epoch_type = 'rest';
    params_opt = 11;
    [events_session, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    events{ii_exp} = events_session;
end
T.nEvents = cellfun(@length,events)'; % note this is without filtering sequence by features!
sortrows( groupsummary(T,'bat_num',["median","mean","max","sum"],"nEvents"),"sum_nEvents", 'descend')

%% pool data
events = [events{:}];

%% forward by score (cummulative thr)
thrs = linspace(0,1,25);
prc_forward = zeros(size(thrs));
for ii_thr = 1:length(thrs)
    seqs = [events.seq_model];
    seqs([seqs.score] < thrs(ii_thr)) = [];
    prc_forward(ii_thr) = mean([seqs.forward]);
end
figure
plot(thrs,prc_forward,'.-');
ylim([0 1]);
xlabel('Replay score thr')
ylabel('Fraction of forward replays');

%% forward by score (by bins)
seqs = [events.seq_model];
EDGES = linspace(min([seqs.score]), max([seqs.score]), 25);
CENTERS = edges2centers(EDGES);
BINS = discretize([seqs.score],EDGES);
fraction_forward = splitapply(@(x)(mean([x.forward])), seqs, BINS);
figure
bar(CENTERS, fraction_forward);
ylim([0 1])
xlabel('Replay score')
ylabel('Fraction of forward replays');

%% apply inclusion criteria 
seqs = [events.seq_model];
seqs = decoding_apply_seq_inclusion_criteria(seqs);
% remove high compression outliers (TODO: deal with those in a better way
% than just removing it)
[seqs([seqs.compression]>25).compression] = deal(nan);

%%
X = [];
X = [X; seqs.middle_pos_norm];
X = [X; seqs.distance];
X = [X; seqs.duration];
X = [X; seqs.compression];
X = [X; seqs.score];
X = [X; seqs.confidence_sparsity];
X = X';
labels = ["position","distance","duration","compression","score","confidence"];
g = categorical([seqs.forward],[true false],{'Forward','Reverse'});

%% bagplot (bi-variate boxplot)
% figure
% hold on
% res=bagplot(X(g=='Forward',[4 5]),'colorbag',0.8*[0 0 1],'colorfence',.95*[0 0 1],'databag',0,'datafence',0);
% xlim([0 30])
% ylim([0.5 1])
% figure
% res=bagplot(X(g=='Reverse',[4 5]),'colorbag',0.8*[1 0 0],'colorfence',.95*[1 0 0],'databag',0,'datafence',0);
% xlim([0 30])
% ylim([0.5 1])

%% options
cmap = brewermap(100,'RdBu');
nbins = [25 25];
sigma = 2;

%% score vs duration
axes(panels_B(1));
cla
hold on
x = [seqs.duration];
y = [seqs.score];
density_diff(x,y,g,nbins,sigma,cmap);
xlabel('Duration (s)')
ylabel('Score')
text(-0.25,1.1, 'B', 'Units','normalized','FontWeight','bold');

%% score vs distance
axes(panels_B(2));
cla
hold on
x = [seqs.distance];
y = [seqs.score];
density_diff(x,y,g,nbins,sigma,cmap);
xlabel('Distance (m)')
ylabel('Score')

%% score vs compression
axes(panels_B(3));
cla
hold on
x = [seqs.compression];
y = [seqs.score];
density_diff(x,y,g,nbins,sigma,cmap);
xlabel('Compression')
ylabel('Score')

%% Distance vs duration
axes(panels_B(4));
cla
hold on
x = [seqs.duration];
y = [seqs.distance];
density_diff(x,y,g,nbins,sigma,cmap);
xlabel('Duration (s)')
ylabel('Distance (m)')

%% Distance vs duration
axes(panels_B(5));
cla
hold on
x = [seqs.duration];
y = [seqs.compression];
density_diff(x,y,g,nbins,sigma,cmap);
xlabel('Duration (s)')
ylabel('Compression')
box on

%%
axes(panels_A(1));
cla
hold on
h=histogram(g,'Normalization','probability');
h.BarWidth = 0.7;
h.FaceColor = 0.5*[1 1 1];
ylabel('Fraction')
text(-0.4,1.1, 'A', 'Units','normalized','FontWeight','bold');
binom_pval = myBinomTest(sum(g=='Forward'),length(g),0.5,'two');
disp('Panel B:');
fprintf('Binomial test (two-sided), pal=%.2g\n',binom_pval);

%% correlation of distance vs duration
x = [seqs.duration]';
y = [seqs.distance]';
[ccc pval] = corr(x,y);
fprintf('\nDistance vs. duration correlations:\n');
fprintf('r = %.2g\n',ccc);
fprintf('p = %.2g\n',pval);

IX = g=='Forward';
[ccc pval] = corr(x(IX),y(IX));
fprintf('\nDistance vs. duration correlations (only for forward replay):\n');
fprintf('r = %.2g\n',ccc);
fprintf('p = %.2g\n',pval);

IX = g=='Reverse';
[ccc pval] = corr(x(IX),y(IX));
fprintf('\nDistance vs. duration correlations (only for reverse replay):\n');
fprintf('r = %.2g\n',ccc);
fprintf('p = %.2g\n',pval);

%% correlation of distance vs compression
x = [seqs.compression]';
y = [seqs.distance]';
[ccc pval] = corr(x,y);
fprintf('\nDistance vs. compression correlations:\n');
fprintf('r = %.2g\n',ccc);
fprintf('p = %.2g\n',pval);

IX = g=='Forward';
[ccc pval] = corr(x(IX),y(IX));
fprintf('\nDistance vs. compression correlations (only for forward replay):\n');
fprintf('r = %.2g\n',ccc);
fprintf('p = %.2g\n',pval);

IX = g=='Reverse';
[ccc pval] = corr(x(IX),y(IX));
fprintf('\nDistance vs. compression correlations (only for reverse replay):\n');
fprintf('r = %.2g\n',ccc);
fprintf('p = %.2g\n',pval);


%% gplotscatter
fig_gplotmat_scatter = figure;
% fig_gplotmat_scatter.WindowState = 'maximized';
fig_gplotmat_scatter.Units = 'centimeters';
fig_gplotmat_scatter.Position = [0 0 25 25];
h=gplotmatrix(X,[],g,["g","b"],[],[],true,[],labels)
sgtitle(['bats ' char(strjoin(""+bats,', '))],'Interpreter','none');

%% gplotdensity
sigma = 2;
nfeat = size(X,2);
fig_gplotmat_density = figure;
% fig_gplotmat_density.WindowState = 'maximized';
fig_gplotmat_density.Units = 'centimeters';
fig_gplotmat_density.Position = [0 0 25 25];
h=tiledlayout(nfeat,nfeat);
h.TileSpacing = 'tight';
cmap=brewermap(100,'RdBu');
for jj_feaure = 1:nfeat
    for ii_feaure = 1:nfeat
        nexttile
        if ii_feaure == jj_feaure
            %% violin
            h=violinplot(X(:,ii_feaure),g);
            [h.ViolinColor] = disperse(cmap([end 1],:)');
            ylabel(labels(ii_feaure));
        else
            %% density diff
            x = X(:,ii_feaure);
            y = X(:,jj_feaure);
            density_diff(x,y,g,nbins,sigma,cmap);
            xlabel(labels(ii_feaure));
            ylabel(labels(jj_feaure));
        end
    end 
end

%% save fig(s)
if exist('bats_to_include','var')
    bats_str = ['bats_' char(strjoin(""+bats,'_'))];
else
    bats_str = 'bats_all';
end

fig_name_out = fullfile(res_dir, sprintf('%s_%s_%s',fig_name_str, epoch_type, bats_str));
print(fig, fig_name_out, '-dpdf', '-cmyk', '-painters');

fig_name_out = fullfile(res_dir, sprintf('%s_%s_%s_%s',fig_name_str, epoch_type, bats_str,'gplotmat_scatter'));
saveas(fig_gplotmat_scatter, fig_name_out, 'tif');

fig_name_out = fullfile(res_dir, sprintf('%s_%s_%s_%s',fig_name_str, epoch_type, bats_str,'gplotmat_density'));
saveas(fig_gplotmat_density, fig_name_out, 'tif');

disp('figure was successfully saved to pdf/tiff/fig formats');
diary off

%%
figure
subplot(211)
hold on
x = [seqs.start_pos_norm];
histogram(x(g=='Forward'),linspace(0,1,100),'Normalization','probability');
histogram(x(g=='Reverse'),linspace(0,1,100),'Normalization','probability');
legend('Forward','Reverse')
xlabel('Start position (norm.)')
ylabel('Probability')
subplot(212)
hold on
x = [seqs.end_pos_norm];
histogram(x(g=='Forward'),linspace(0,1,100),'Normalization','probability');
histogram(x(g=='Reverse'),linspace(0,1,100),'Normalization','probability');
legend('Forward','Reverse')
xlabel('End position (norm.)')
ylabel('Probability')

%%
function density_diff(x,y,g,nbins,sigma,cmap)
    %% density diff
    xedges = linspace(min(x),max(x),nbins(1));
    yedges = linspace(min(y),max(y),nbins(2));
    N1 = histcounts2(x(g=='Forward'),y(g=='Forward'),xedges,yedges,'Normalization','pdf')';
    N2 = histcounts2(x(g=='Reverse'),y(g=='Reverse'),xedges,yedges,'Normalization','pdf')';
    N = N1-N2;
    N = imgaussfilt(N,sigma,'Padding','replicate');
    xcenters = edges2centers(xedges);
    ycenters = edges2centers(yedges);
    h=imagesc(xcenters,ycenters,N);
    axis xy
    axis tight
    colormap(cmap);
    hax=gca;
    m1 = min(N,[],'all');
    m2 = max(N,[],'all');
%     m = max(abs([m1 m2]));
    m = abs(m1);
%     m = abs(m2);
%     m = mean(abs([m1 m2]));
    hax.CLim = [-m m];
end







