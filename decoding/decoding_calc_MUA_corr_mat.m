%%
clear
clc
% close all

%%
% exp_ID = 'b0184_d191124';
% exp_ID = 'b0184_d191125';
% exp_ID = 'b0184_d191126';
% exp_ID = 'b0184_d191127';
exp_ID = 'b0184_d191205';
% exp_ID = 'b0184_d191208';
dir_IN = 'F:\sequences\proc\';
dir_OUT = 'F:\sequences\figures\MUA_corr_mat';

%% load data
exp = exp_load_data(exp_ID,'details','path','LM');
load(fullfile(dir_IN,[exp_ID '.mat']));

%% params
pos_bin_size = 1;
% pos_std = 0.2; % does not apply here
mark_bin_size = 40;
mark_std = 20;
mark_limits = [0 750];

% pos_limits = [min(position) max(position)];
% pos_edges = pos_limits(1):pos_bin_size:pos_limits(end);
% [~,~,pos_BIN] = histcounts(position,pos_edges);
[N,pos_edges,pos_BIN] = histcounts(position,'BinWidth',pos_bin_size);
pos_centers = edges2centers(pos_edges);
% gDir = categorical(direction,[1 -1],["map1","map2"]);
direction12 = direction;
direction12(direction12==-1)=2;
is_FE = any(t>FE_ti(:,1)&t<FE_ti(:,2),1);

m1 = min(multiunits,[],"all");
m2 = max(multiunits,[],"all");
mark_edges = m1:mark_bin_size:m2;
mark_centers = edges2centers(mark_edges);
mark_centers = mark_limits(1):mark_bin_size:mark_limits(end);
xii = {};
for ch = 1:4
%     xii{ch} = min(x(:,ch)):mark_bin_size:max(x(:,ch));
    xii{ch} = mark_centers;
end
[x1,x2,x3,x4] = ndgrid(xii{:});
x1 = x1(:,:)';
x2 = x2(:,:)';
x3 = x3(:,:)';
x4 = x4(:,:)';
xi = [x1(:) x2(:) x3(:) x4(:)];

%%
nDir = 2;
nPos = length(pos_centers);
nMark = length(mark_centers);
nTT = size(multiunits,3);
nCh = size(multiunits,2);
fprintf('f size = %.2f MB\n', nMark^nCh*nDir*nPos*nTT*8*1e-6);

%%
f = nan(nMark^nCh,nDir,nPos,nTT);
tic
for ii_dir = 1:nDir
    TF_dir = is_FE & (direction12==ii_dir);
    for ii_pos_bin = 1:nPos
        TF = TF_dir & (pos_BIN==ii_pos_bin);
        parfor ii_TT = 1:nTT
            fprintf('dir: %d, pos_bin: %d, TT: %d\n',ii_dir,ii_pos_bin,ii_TT);
            x = multiunits(TF,:,ii_TT);
            x(any(isnan(x),2),:)=[];
            if isempty(x)
                continue;
            end
            f(:,ii_dir,ii_pos_bin,ii_TT) = mvksdensity(x,xi,'Bandwidth',mark_std);
        end
    end
end
toc

%%
f2=f;
f2 = permute(f2,[3 2 1 4]);
f2 = reshape(f2,nDir*nPos,[])';
ccc = corr(f2,'type','Pearson');
mask = eye(size(ccc),'logical');
% ccc(mask) = nan;
figure
imagesc(ccc)
colorbar
% colormap gray
hax=gca;
% hax.CLim = ([0.92 1]);
hax.CLim(1) = prctile(ccc(:),50);
sgtitle({"MUA map corr - " + exp_ID},'interpreter','none')
xlabel('Position x Direction')
ylabel('Position x Direction')
xticks([])
yticks([])
file_out = sprintf('%s_MUA_corr_mat',exp_ID);
file_out = fullfile(dir_OUT,file_out);
saveas(gcf,file_out,'jpg');

%%
return

%%
x = randn(100,2);
xi1 = linspace(-1,1,5);
xi2 = linspace(-2,2,10);
[XI1,XI2] = ndgrid(xi1,xi2);
xi = cat(3,XI1,XI2);
xi = reshape(xi,[],2);
% f = mvksdensity(x,xi,'Bandwidth',0.1);

%%
load hald
xi = ingredients(1:3,:);
f = mvksdensity(ingredients,xi,'Bandwidth',0.8);

%%
gridx1 = 0:2:22;
gridx2 = 20:5:80;
gridx3 = 0:2:24;
gridx4 = 5:5:65;
[x1,x2,x3,x4] = ndgrid(gridx1,gridx2,gridx3,gridx4);
