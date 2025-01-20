%%
clear
clc
close all

%%
out_dir = 'E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\cell_revision';

%%
data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay pos between 25-115, dec acc_gt_70%.mat';
data = load(data_filename);

%% 
exp_list = unique(data.data_info.exp_ID(data.TF)');
X = [];
Y = [];
ID = [];
for ii_exp = 1:length(exp_list)
    %%
    exp_ID = exp_list{ii_exp};
    exp = exp_load_data(exp_ID, 'details', 'pos');
    co = exp.pos.proc_1D.co;
    x = co.pos(1:end-1);
    y = co.pos(2:end);
    id = repelem(ii_exp,length(x))';
    X = [X;x];
    Y = [Y;y];
    ID = [ID;id];
end

%%
fig1 = figure(Units="centimeters",Position=[5 5 30 30]);
tiledlayout(2,2,'TileSpacing','compact','Padding','loose');
nexttile
hold on
plot(X,Y,'o')
axis equal
h=lsline;h.Color = 'r';
[r,r_pval] = corr(X,Y,'rows','pairwise','type','Pearson');
[rho,rho_pval] = corr(X,Y,'rows','pairwise','type','Spearman');
xlabel('Crossover i position (m)')
ylabel('Crossover i+1 position (m)')
title({ ...
    'Temporal structure in cross-overs?'; ...
    sprintf('pearson, r=%.2g, P=%.2g',r,r_pval); ...
    sprintf('spearman, rho=%.2g, P=%.2g',rho,rho_pval); ...
    })
xlim([0 130])
ylim([0 130])
h=refline(1,0);h.Color = [1 1 1]*0.5;

nexttile
hold on
scatter(X,Y,30,ID,'filled','o')
hax=gca;
cmap = brewermap(length(unique(ID)),'Set1');
hax.Colormap = cmap;
axis equal
xlabel('Crossover i position (m)')
ylabel('Crossover i+1 position (m)')
title('colors = different sessions')
xlim([0 130])
ylim([0 130])
h=refline(1,0);h.Color = [1 1 1]*0.5;

nbins = 30;
ker_sd = 2;
N = histcounts2(X,Y,nbins);
N = N';
ker = fspecial('gaussian',round(5*ker_sd),ker_sd);
N2 = imfilter(N,ker);
nexttile
imagesc(N)
axis equal
axis off
axis xy
xlabel('Crossover i position (m)')
ylabel('Crossover i+1 position (m)')
title(sprintf('Density (nbins = %d)',nbins))

nexttile
imagesc(N2)
axis equal
axis off
axis xy
xlabel('Crossover i position (m)')
ylabel('Crossover i+1 position (m)')
title(sprintf('Density (smoothed, sd=%d bins)',ker_sd))

filename = fullfile(out_dir, '2bats_crossover_dynamics.pdf');
exportgraphics(fig1, filename)

%%
fig2 = figure(Units="centimeters",Position=[5 5 10 10]);
tiledlayout(1,1,'TileSpacing','compact','Padding','loose');
nexttile
hold on
scatter(X,Y,30,ID,'filled','o')
hax=gca;
cmap = brewermap(length(unique(ID)),'Set1');
hax.Colormap = cmap;
axis equal
xlabel('Crossover i position (m)')
ylabel('Crossover i+1 position (m)')
% title('colors = different sessions')
xlim([0 130])
ylim([0 130])
h=refline(1,0);h.Color = [1 1 1]*0.5;
filename = fullfile(out_dir, 'Figure_R3_crossover_temporal_structure.pdf');
exportgraphics(fig2, filename)

















%%